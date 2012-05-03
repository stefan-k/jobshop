################################################################################
## LIBRARY FOR EVOLUTIONARY ALGORITHMS                                        ##
################################################################################
#
# Authors: Joris Bayer
#          Stefan Kroboth
#
# TODO: 
#  * adopt design to better work with evolutionary strategies not just genetic 
#    algorithms.
#  * more printing functions (human readable)
#  * possibility to flag chromosomes as 'bad'
#  * test BitGene and adjust mutation/recombination for BitGene
#  * RFC: BitGene or BinaryGene
#  * different selection/competition models 
#    - replace parents
#    - in competition with parents
#    - decision wether a child becomes part of the population or not
#    - in short: all that (kappa, lambda) stuff (Papers?)
#  * design similar to the one shown in prototype.jl
#    - although in a more general way
#  * stopping criterion
#  * should Vector rather be AbstractArray{T,1} ?
#  * when adding Chromosome to Population, check that the number of
#    genes per chromosome are the same
#  * It could be that sort doesn't need to be implemented as long as
#    isless is defined for a certain type -- we should check that
#  * [PRIORITY] Remove the dirty hack in roulette!
#  * Implement adapting the standard deviation of genes
#  * [PRIORITY] Make gray2binary work for bitstrings
#  * discuss about how to split up this file
#  * [NOT IMPORTANT] try to get the type constraints to work
#  * add more to the todo list
#
# DONE:
#  * [DONE] Create random Gene with constraints
#  * [DONE] Create random Chromosome
#  * [DONE] properly implement sort!() for the population
#  * [DONE] handling of objective functions
#    - [DONE] new idea: get rid of the objective function within types entirely
#  * [DONE] roulette wheel selection 
#  * [DONE] crossover
#  * [DONE] mutation
#  * [DONE] Create type BitGene for "classical" genetic algorithms
#    - [DONE] bits() could be helpful for that
#    - [DONE] check if there is a way to convert to gray code

# I have no idea what I'm doing
abstract AbstractEvolutionary
abstract AbstractGene                 <: AbstractEvolutionary
abstract AbstractChromosome           <: AbstractEvolutionary
abstract AbstractPopulation           <: AbstractEvolutionary
abstract AbstractGenerations          <: AbstractEvolutionary
abstract AbstractGeneticProbabilities <: AbstractEvolutionary
abstract AbstractGeneticOperation     <: AbstractEvolutionary

abstract Recombination <: AbstractGeneticOperation
abstract Mutation      <: AbstractGeneticOperation
abstract Reproduction  <: AbstractGeneticOperation
abstract Immigration   <: AbstractGeneticOperation

################################################################################
## GENE TYPE                                                                  ##
################################################################################

type Gene <: AbstractGene
    gene::Number
    std::Float64
    upper_limit::Float64
    lower_limit::Float64
    # Point to discuss: Should 'factor' be a part of the Gene type?
    # Could be individually adapted, if necessary...
    function Gene(gene::Number) 
        new(gene, 0.5, NaN, NaN)
    end
    
    function Gene(gene::Number, std::Float64) 
        new(gene, std, NaN, NaN)
    end
    
    function Gene(gene::Number, std::Float64, upper_limit::Float64, 
                  lower_limit::Float64)
        if lower_limit >= upper_limit
            error("lower_limit must be less than upper_limit")
        end
        new(gene, std, upper_limit, lower_limit)
    end
end

# Copy-constructor
function copy(gene::Gene) 
    Gene(copy(gene.gene), 
         copy(gene.std),
         copy(gene.upper_limit),
         copy(gene.lower_limit))
end

# Utility functions
std(gene::Gene) = gene.std

function *(factor::Number, gene::Gene)
    g = copy(gene)
    g.gene *= factor
    return g
end

# increase standard deviation
function broaden_std(gene::Gene, factor::Float64)
    gene.std *= factor
end

# decrease standard deviation
function narrow_std(gene::Gene, factor::Float64)
    gene.std /= factor
end

# create random Gene
rand(T::Type{Gene}) = Gene(rand(), rand())
rand(T::Type{Gene}, std::Float64) = Gene(rand(), std)
rand(T::Type{Gene}, upper_limit::Float64, lower_limit::Float64) = Gene(rand()*(upper_limit-lower_limit)+lower_limit, rand(), upper_limit, lower_limit)
# create value in the given limit
rand(T::Type{Gene}, std::Float64, upper_limit::Float64, lower_limit::Float64) = Gene(rand()*(upper_limit-lower_limit)+lower_limit, std, upper_limit, lower_limit)


################################################################################
## BIT GENE TYPE                                                              ##
################################################################################

type BitGene <: AbstractGene
    gene::ASCIIString # in gray code
    gene_type::Type   # for converting back

    function BitGene(gene_type::Type, gene::ASCIIString)
        if 8*sizeof(gene_type) != length(gene)
            error("length of ASCII string does not agree with provided type")
        end
        new(gene, gene_type)
    end

    function BitGene{T<:Integer}(gene::T)
        new(int2gray(gene), T)
    end
end

rand(T::Type{BitGene}) = BitGene(randi(Int64))
rand(T::Type{BitGene}, U::Type) = BitGene(randi(U))

# Binary utility functions
int2gray(num::Integer) = bits(num $ (num >> 1))

function binary2int(bs::ASCIIString)
    tmp = copy(bs)
    s_len = length(tmp)
    out = 0
    for i = s_len-1:-1:0
        # this feels so dirty...
        out += 2^i*parse_int(ASCIIString([uint8(tmp[s_len-i])]))
    end
    return out
end

function gray2binary(num::Integer) 
    tmp = copy(num)
    shift = 1
    while shift < 8*sizeof(num)
        tmp = tmp $ (tmp >> shift)
        shift *= 2
    end
    return tmp
end

gray2binary(s::ASCIIString) = gray2binary(binary2int(s))

################################################################################
## CHROMOSOME TYPE                                                            ##
################################################################################

type Chromosome <: AbstractChromosome
    genes::Vector
    length::Int64
    fitness::Float64

    function Chromosome(genes::Vector{Gene}) 
        new(map(copy, genes), length(genes), Inf)
    end

    function Chromosome(genes::Vector{Gene}, fitness::Float64)
        new(copy(genes), length(genes), copy(fitness))
    end

    function Chromosome()
        new(Gene[], 0, Inf)
    end
end

# Copy-constructor
function copy(chr::Chromosome)
    Chromosome(copy(chr.genes), 
               copy(chr.fitness))
end

# referencing
ref(chromosome::Chromosome, ind...) = chromosome.genes[ind...]

# Utility functions
size(chromosome::Chromosome) = chromosome.length
length(chromosome::Chromosome) = chromosome.length
# print the stds of all genes of a chromosome indicated by a range/index
std(chromosome::Chromosome, ind...) = [std(gene) | gene = chromosome[ind...]]

# replace Gene of a Chromosome
function assign(chr::Chromosome, g::Gene, idx::Int64)
    chr.genes[idx] = copy(g)
end

# replace Genes of a Chromosome
function assign(chr::Chromosome, gs::Vector{Gene}, idx::Range1{Int64})
    chr.genes[idx] = map(copy, gs)
end

# replace Genes of a Chromosome
function assign(chr::Chromosome, gs::Vector{Gene}, idx::Vector{Int64})
    chr.genes[idx] = map(copy, gs)
end

# push gene onto the chromosome
function push(chr::Chromosome, g::Gene)
    push(chr.genes, copy(g))
    chr.length += 1
end

(+)(chr::Chromosome, g::Gene) = push(chr, g)

# apply an objective function to a chromosome
function apply_obj_func(chromosome::Chromosome, obj_func::Function)
    # untested
    chromosome.fitness = obj_func(chromosome)
end

# increase standard deviation of all genes in Chromosome
function broaden_std(chromosome::Chromosome, factor::Float64)
    for i = 1:length(chromosome)
        # this can probably be done with map
        broaden_std(chromosome[i], factor)
    end
end

# decrease standard deviation of all genes in Chromosome
function narrow_std(chromosome::Chromosome, factor::Float64)
    for i = 1:length(chromosome)
        narrow_std(chromosome[i], factor)
    end
end

rand(T::Type{Chromosome}, num::Int64, x...) = Chromosome([rand(Gene, x...) | i=1:num]) # neat

function rand(T::Type{Chromosome}, num::Int64, obj_func::Function, x...) 
    chr = rand(T, num, x...)
    chr.fitness = obj_func(chr)
    return chr
end

function print(chromosome::Chromosome)
    for i = 1:length(chromosome.genes)
        print("|")
        print(chromosome[i].gene)
    end
    print("|")
end

# Multiply all genes in the chromosome by a scalar
function *(factor::Number, chromosome::Chromosome)
    c = copy(chromosome)
    genes = Gene[]
    for g in c.genes
        push(genes,factor*g)
    end
    c.genes = genes

    return c
end

################################################################################
## POPULATION TYPE                                                            ##
################################################################################

type Population <: AbstractPopulation
    chromosomes::Vector
    pop_size::Int64

    # several chromosomes passed
    function Population(chromosomes::Vector) 
        new(chromosomes, length(chromosomes))
        # should we copy that too?
    end

    # one chromosome passed
    function Population(chromosome::Chromosome) 
        new([chromosome], 1)
        # should we copy that too?
    end

    # empty Population
    function Population()
        new(Chromosome[], 0)
    end
end

# Copy-constructor
function copy(pop::Population)
    Population(copy(pop.chromosomes))
end

# referencing
ref(population::Population, ind...) = population.chromosomes[ind...]

# Utility functions
size(population::Population) = population.pop_size
length(population::Population) = population.pop_size

# overload (+) for easier appending
(+)(pop::Population, chr::Chromosome) = push(pop, chr)

# sum of the fitness of all chromosomes of a population
function fitness_sum(population::Population) 
    fitness_sum = 0
    for i=1:length(population)
        fitness_sum += population[i].fitness
    end
    return fitness_sum
end

# sum of the inverse of fitness of all chromosomes of a population
# needed for roulette
function inv_fitness_sum(population::Population) 
    fitness_sum = 0
    for i=1:length(population)
        fitness_sum += 1/population[i].fitness
    end
    return fitness_sum
end

function print_population(population::Population)
    for i=1:length(population)
        tmp = population[i]
        print("$tmp \n")
    end
end

function print_population_evo(population::Population)
    for i=1:length(population)
        tmp = population[i]
        print("$tmp ")
    end
    print("\n")
end

# Modifiers
function push(population::Population, chromosome::Chromosome)
    push(population.chromosomes, copy(chromosome))
    population.pop_size += 1
end

function isless(chr1::Chromosome, chr2::Chromosome)
    return chr1.fitness < chr2.fitness
end

ismore(chr1::Chromosome, chr2::Chromosome) = !isless(chr1, chr2)

# in-place sort
function sort!(population::Population)
    sort!(isless, population.chromosomes)
end

# sort
function sort(population::Population)
    sort(isless, population.chromosomes)
end

# in-place reverse sort
function sortr!(population::Population)
    sort!(ismore, population.chromosomes)
end

# reverse sort
function sortr(population::Population)
    sort(ismore, population.chromosomes)
end

# increase standard deviation of all genes in Chromosome
function broaden_std(population::Population, factor::Float64)
    for i = 1:length(population)
        # this can probably be done with map
        broaden_std(population[i], factor)
    end
end

# decrease standard deviation of all genes in Chromosome
function narrow_std(population::Population, factor::Float64)
    for i = 1:length(population)
        narrow_std(population[i], factor)
    end
end

# well, that one was easy.
rand(T::Type{Population}, chr_num::Int64, gene_num::Int64, obj_func::Function, x...) = Population([rand(Chromosome, gene_num, obj_func, x...) | i = 1:chr_num])
rand(T::Type{Population}, chr_num::Int64, gene_num::Int64, x...) = Population([rand(Chromosome, gene_num, x...) | i = 1:chr_num])

################################################################################
## GENERATIONS TYPE                                                           ##
################################################################################

# keeps track of generations
type Generations <: AbstractGenerations
    populations::Array
    generations::Int64

    function Generations(population)
        new([copy(population)], 1)
    end

    function Generations()
        new(Population[], 0)
    end
end

# Copy-Constructor
# does this even make sense for Generations?
function copy(gen::Generations)
    Generations(copy(populations),
                copy(generations))
end

# referencing
ref(generations::Generations, ind...) = generations.populations[ind...]

# Utility functions
get_generations(generations::Generations) = generations.generations

# Modifiers
function push(generations::Generations, population::Population)
    push(generations.populations, copy(population))
    generations.generations += 1
end

(+)(generations::Generations, population::Population) = push(generations, population)

################################################################################
## DATASTRUCTURE FOR GENETIC ALGORITHM PROBABILITIES                          ##
################################################################################

type GeneticProbabilities <: AbstractGeneticProbabilities
    mutation::Float64
    recombination::Float64
    reproduction::Float64
    immigration::Float64

    function GeneticProbabilities(mutation::Float64, 
                                  recombination::Float64,
                                  reproduction::Float64,
                                  immigration::Float64)
        sum = mutation + recombination + reproduction + immigration
        new(mutation/sum, recombination/sum, reproduction/sum, immigration/sum)
    end
end

get_vector(gp::GeneticProbabilities) = [gp.mutation, gp.recombination, gp.reproduction, gp.immigration]


################################################################################
## FUNCTIONS                                                                  ##
################################################################################

# Roulette Wheel Selection on Population
function roulette(pop::Population)
    #sort!(pop) # I don't think sorting is necessary and might even lead to 
                # problems... gotta check that
    f_sum = inv_fitness_sum(pop)
    idx = rand()*f_sum
    x = 0
    elem = 1
    for i=1:length(pop)
        # DIRTY DIRTY HACK!
        # This is supposed to work without abs(), but this leads to problems when
        # negative fitness is allowed... we have to find a fix for this.
        x += abs(1/pop[i].fitness)
        if idx < x
            return elem
        end
        elem += 1
    end
    error("weird error that should not happen. You probably didn't define a fitness.")
end

# return several indices determined by roulette
roulette(pop::Population, num::Int64) = [ roulette(pop) | i = 1:num ]

# Roulette Wheel Selection on general Vectors (in case probabilies are passed)
function roulette(p::Vector{Float64})
    # Am I thinking wrong? Is sorting even necessary?
    #prop = sortr(p) # do not sort in-place!
    prop = copy(p)
    prop_sum = sum(prop) # In case it doesn't sum up to 1
    idx = rand()*prop_sum
    x = 0
    elem = 1
    for i=1:length(prop)
        x += prop[i]
        if idx < x
            return elem
        end
        elem += 1
    end
    error("dafuq?")
end

# roulette wheel selection on GeneticProbabilities
function roulette(gp::GeneticProbabilities) 
    idx = roulette(get_vector(gp))
    if idx == 1
        return Mutation
    elseif idx == 2
        return Recombination
    elseif idx == 3
        return Reproduction
    elseif idx == 4
        return Immigration
    else
        error("Error: Impossible Error.")
    end
end

# make sure the gene doesn't exceed it's limits
function assess_limits(g::Gene)
    # works with NaNs!
    if g.gene > g.upper_limit
        g.gene = g.upper_limit
    elseif g.gene < g.lower_limit
        g.gene = g.lower_limit
    end
end

# mutate a single gene
function mutate!(g::Gene)
    g.gene += g.std*randn()
    assess_limits(g)
end

function mutate(g::Gene)
    gn = copy(g)
    gn.gene += gn.std*randn()
    assess_limits(gn)
    return gn
end

# mutate a single gene with a predefined standard deviation
# -> ignores std settings in g.gene
function mutate!(g::Gene, std::Float64)
    g.gene += std*randn()
    assess_limits(g)
end

function mutate(g::Gene, std::Float64)
    gn = copy(g)
    gn.gene += std*randn()
    assess_limits(gn)
    return gn
end

# mutate a chromosome
function mutate!(chr::Chromosome)
    for i=1:length(chr)
        mutate!(chr[i])
    end
end

function mutate(chr::Chromosome)
    chrn = Chromosome()
    for i=1:length(chr)
        chrn + mutate(chr[i])
    end
    return chrn
end

# mutate a chromosome with a predefined standard deviation
# -> ignores std settings in g.gene
function mutate!(chr::Chromosome, std::Float64)
    for i=1:length(chr)
        mutate!(chr[i], std)
    end
end

function mutate(chr::Chromosome, std::Float64)
    chrn = Chromosome()
    for i=1:length(chr)
        chrn + mutate(chr[i], std)
    end
    return chrn
end

# mutate a whole population (not sure if anyone will ever need this)
function mutate!(pop::Population)
    for i=1:length(pop)
        mutate!(pop[i])
    end
end

function mutate(pop::Population)
    popn = Population()
    for i=1:length(pop)
        popn + mutate(pop[i])
    end
    return popn
end

# mutate a whole population with predefined std (not sure if anyone will ever need this)
function mutate!(pop::Population, std::Float64)
    for i=1:length(pop)
        mutate!(pop[i], std)
    end
end

function mutate(pop::Population, std::Float64)
    popn = Population()
    for i=1:length(pop)
        popn + mutate(pop[i], std)
    end
    return popn
end

# in-place crossover
function crossover!(chr1::Chromosome, chr2::Chromosome, slices::Int64)
    @assert length(chr1) == length(chr2)
    # weird, even works when rand produces 0
    idx = sort([1, int64(round(length(chr1) * rand(slices))), length(chr1)+1])
    tmp = copy(chr1)
    for i=1:length(idx)-1
        if i%2 == 0
            range = [idx[i]:idx[i+1]-1]
            chr1[range] = map(copy, chr2[range])
            chr2[range] = map(copy, tmp[range])
        end
    end
end

function crossover!(pop::Population, slices::Int64)
    idx = roulette(pop, 2)
    crossover!(pop[idx[1]], pop[idx[2]], slices)
end

# crossover
function crossover(chr1::Chromosome, chr2::Chromosome, slices::Int64)
    @assert length(chr1) == length(chr2)
    # weird, even works when rand produces 0
    idx = sort([1, int64(round((length(chr1)-1) * rand(slices)))+1, length(chr1)+1])
    #tmp = copy(chr1)
    chr1n = copy(chr1)
    chr2n = copy(chr2)
    for i=1:length(idx)-1
        if i%2 == 0
            range = [idx[i]:idx[i+1]-1]
            chr1n[range] = map(copy, chr2[range])
            chr2n[range] = map(copy, chr1[range])
        end
    end
    return chr1n, chr2n
end

function crossover(pop::Population, slices::Int64)
    idx = roulette(pop, 2)
    return crossover(pop[idx[1]], pop[idx[2]], slices)
end

crossover(chr1::Chromosome, chr2::Chromosome) = crossover(chr1, chr2, 2)
crossover(pop::Population) = crossover(pop, 2)

################################################################################
## GENETIC ALGORITHM                                                          ##
################################################################################

function genetic(pop::Population, probabilities::GeneticProbabilities, iter::Int64, obj_func::Function)
    # First version of a genetic algorithm - pretty basic, needs a lot more functionality and
    # probably even a better design and more flexibility.

    pop_o = copy(pop) # prevent in-place fiasco
    obj_func(pop_o) # make sure the the fitness for every chromosome is available
    
    # NOT keeping all generations around is a bit faster
    #gen = Generations() # Discussion: This might get pretty big... what to do?
    for j = 1:iter # max generations
        pop_n = Population()
        for i = 1:length(pop_o)
            operation = roulette(probabilities)
            if operation == Mutation
                chr = mutate(pop_o[roulette(pop_o)])
            elseif operation == Recombination
                chr = crossover(pop_o)
                chr = chr[1] # for future: take both offsprings!
            elseif operation == Reproduction
                chr = copy(pop_o[roulette(pop_o)])
            elseif operation == Immigration
                chr = rand(Chromosome, length(pop_o[1]), obj_func)
            end
            #[obj_func(chr[k]) | k=1:length(chr)] # for recombination?
            #obj_func(chr)
            pop_n + chr
        end
        obj_func(pop_n)
        sort!(pop_n)
        
        # Print best chromosome
        print("Best chromosome of generation $(dec(j,3)): ")
        print(pop_n[1])
        
        #gen + pop_n
        pop_o = copy(pop_n)
    end
end

################################################################################
## EVOLUTIONARY ALGORITHM                                                     ##
################################################################################

function evo_1plus1(pop::Population, epsilon::Number, factor::Float64, 
                    prop_pos::Float64, iter::Integer, inner_iter::Integer, 
                    obj_func::Function)
    # doesn't fit perfectly in the current design, therefore each chromosome
    # is assumed to only consist of one Gene.
    # We could also just let the User pass one Chromosome, which doesn't fit 
    # nicely in the theory which is defined as only to act on the Chromosome.

    pop_o = copy(pop)
    prev_fitness = obj_func(pop_o)


    while true # not really clear if this is necessary from the lecture notes
        pos_mut = 0;
        for i=1:iter
            # create new, mutated population
            pop_n = Population()
            for k=1:length(pop_o)
                pop_n + mutate(pop_o[k])
            end
            fitness = obj_func(pop_n)

            # abort mission
            if abs(prev_fitness - fitness) < epsilon
                print("Optimum: ")
                print_population_evo(pop_o)
                return pop_n
            elseif i%inner_iter == 0
                # adapt std after inner_iter iterations
                if pos_mut/i < prop_pos
                    broaden_std(pop_o, factor) # ooops, this might be implementend the wrong way, check
                else
                    narrow_std(pop_o, factor)
                end
            end

            # survival of the fittest
            if fitness < prev_fitness
                pop_o = copy(pop_n)
                print_population_evo(pop_o)
                pos_mut += 1
                prev_fitness = fitness
            end
        end
    end
end

function evo_slash(pop::Population, rho::Integer, lambda::Integer, iter::Integer,
                   factor::Float64, epsilon::Float64, obj_func::Function)
    # rho is the number of parents involved in one descendent. the design does
    # not account for that for now, therefore it is a useless parameter.
    pop_o = copy(pop)
    obj_func(pop_o)
    mu = length(pop_o)
    prev_fitness = Inf

    for i=1:iter
        pop_n = Population()
        while length(pop_n) < lambda
            chrs = crossover(pop_o)
            # randomly adjust std of Genes. Not sure if this is
            # supposed to be done for each new chromosome individually
            # or for the whole new population...
            for j=1:length(chrs)
                if randbit() == 1
                    broaden_std(chrs[j], factor)
                else
                    narrow_std(chrs[j], factor)
                end
                # push mutated thingy
                pop_n + mutate(chrs[j]) 
            end
        end

        obj_func(pop_n)
        sort!(pop_n)
        # slice new population
        pop_o = Population(pop_n[1:mu])

        # stopping criterion
        # not good at all....
        if abs(prev_fitness - pop_o[1].fitness) < epsilon
            print("Optimum: ")
            print(pop_o[1])
            println()
            return pop_o
        end

        print(pop_o[1])
        println()
        prev_fitness = copy(pop_o[1].fitness)
    end
    pop_o
end


