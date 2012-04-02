################################################################################
## LIBRARY FOR EVOLUTIONARY ALGORITHMS                                        ##
################################################################################
#
# Authors: Joris Bayer
#          Stefan Kroboth
#
# TODO: 
#  * [DONE] Create random Gene with constraints
#  * [DONE] Create random Chromosome
#  * [DONE] properly implement sort!() for the population
#  * try to get the type constraints to work
#  * more printing functions (human readable)
#  * possibility to flag chromosomes as 'bad'
#  * [PARTIAL?] handling of objective functions
#  * [DONE] roulette wheel selection 
#    - needs objective function and a working sort
#      + Update: why would it need the objective function? What was I thinking?
#  * [DONE] crossover
#  * [DONE] mutation
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
#  * add more to the todo list

# I have no idea what I'm doing
abstract AbstractEvolutionary
abstract AbstractGene        <: AbstractEvolutionary
abstract AbstractChromosome  <: AbstractEvolutionary
abstract AbstractPopulation  <: AbstractEvolutionary
abstract AbstractGenerations <: AbstractEvolutionary

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
    
    function Gene(gene::Number, std::Float64, upper_limit::Float64, lower_limit::Float64)
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
rand(T::Type{Gene}, upper_limit::Float64, lower_limit::Float64) = Gene(rand(), rand(), upper_limit, lower_limit)
# create value in the given limit
rand(T::Type{Gene}, std::Float64, upper_limit::Float64, lower_limit::Float64) = Gene(rand()*(upper_limit-lower_limit)+lower_limit, std, upper_limit, lower_limit)

################################################################################
## CHROMOSOME TYPE                                                            ##
################################################################################

type Chromosome <: AbstractChromosome
    genes::Vector
    length::Int64
    fitness::Float64
    obj_func::Function

    function Chromosome(genes::Vector{Gene}, obj_func::Function) 
        new(map(copy, genes), length(genes), Inf, obj_func)
        # map ensures that each element of the vector is copied,
        # not just the array itself. copy() would just copy the array
        # and leave all references in genes as the are, which is wrong.

        # the following might not work... 
        # maybe fitness shouldn't be calculated when created.
        # but I think it's a good idea to keep the objective function around.
        # Or maybe not (what if the objective function changes?). 
        # Comments requested ;)
        fitness = obj_func(genes)
    end

    function Chromosome(genes::Vector{Gene})
        print("Warning: No objective function passed!\n")
        new(copy(genes), length(genes), Inf, (x)->()) # hack, I don't know how to pass
                                                      # a 'None' function
    end

    function Chromosome(genes::Vector{Gene}, fitness::Float64)
        print("Warning: No objective function passed!\n")
        new(copy(genes), length(genes), copy(fitness), (x)->())
    end

    function Chromosome(genes::Vector{Gene}, obj_func::Function)
        new(copy(genes), length(genes), Inf, obj_func)
    end

    function Chromosome(genes::Vector{Gene}, fitness::Float64, obj_func::Function)
        new(copy(genes), length(genes), copy(fitness), obj_func)
    end

end

# Copy-constructor
function copy(chr::Chromosome)
    # not sure if function has to be copied
    Chromosome(copy(chr.genes), 
               copy(chr.fitness),
               copy(chr.obj_func))
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

# apply the objective function defined in the Chromosome type
function apply_obj_func(chromosome::Chromosome)
    chromosome.fitness = chromosome.obj_func(chromosome) # not sure if this works, gotta test
end

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

rand(T::Type{Chromosome}, num::Int64) = Chromosome([rand(Gene) | i=1:num]) # neat
rand(T::Type{Chromosome}, num::Int64, obj_func::Function) = Chromosome([rand(Gene) | i=1:num], obj_func) 
rand(T::Type{Chromosome}, num::Int64, std::Float64) = Chromosome([rand(Gene, std) | i=1:num])
rand(T::Type{Chromosome}, num::Int64, std::Float64, obj_func::Function) = Chromosome([rand(Gene, std) | i=1:num], obj_func)
rand(T::Type{Chromosome}, num::Int64, upper_limit::Float64, lower_limit::Float64) = Chromosome([rand(Gene, upper_limit, lower_limit) | i=1:num])
rand(T::Type{Chromosome}, num::Int64, upper_limit::Float64, lower_limit::Float64, obj_func::Function) = Chromosome([rand(Gene, upper_limit, lower_limit) | i=1:num], obj_func)
rand(T::Type{Chromosome}, num::Int64, std::Float64, upper_limit::Float64, lower_limit::Float64) = Chromosome([rand(Gene, std, upper_limit, lower_limit) | i=1:num])
rand(T::Type{Chromosome}, num::Int64, std::Float64, upper_limit::Float64, lower_limit::Float64, obj_func::Function) = Chromosome([rand(Gene, std, upper_limit, lower_limit) | i=1:num], obj_func)


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
end

# Copy-constructor
function copy(pop::Population)
    Population(copy(chromosomes),
               copy(pop_size))
end

# referencing
ref(population::Population, ind...) = population.chromosomes[ind...]

# Utility functions
size(population::Population) = population.pop_size
length(population::Population) = population.pop_size

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

################################################################################
## FUNCTIONS                                                                  ##
################################################################################

function roulette(pop::Population)
    sort!(pop)
    f_sum = inv_fitness_sum(pop)
    idx = rand()*f_sum
    x = 0
    elem = 1
    for i=1:length(pop)
        x += 1/pop[i].fitness
        if idx < x
            return elem
        end
        elem += 1
    end
    error("weird error that should not happen. You probably didn't define a fitness.")
end

function roulette(p::Vector{Float64})
    prop = sortr(p) # do not sort in-place!
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


# return several indices determined by roulette
roulette(pop::Population, num::Int64) = [ roulette(pop) | i = 1:num ]

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
function mutate(g::Gene)
    g.gene += g.std*randn()
    assess_limits(g)
end

# mutate a single gene with a predefined standard deviation
# -> ignores std settings in g.gene
function mutate(g::Gene, std::Float64)
    g.gene += std*randn()
    assess_limits(g)
end

# mutate a chromosome
function mutate(chr::Chromosome)
    for i=1:length(chr)
        mutate(chr[i])
    end
end

# mutate a chromosome with a predefined standard deviation
# -> ignores std settings in g.gene
function mutate(chr::Chromosome, std::Float64)
    for i=1:length(chr)
        mutate(chr[i], std)
    end
end

# mutate a whole population (not sure if anyone will ever need this)
function mutate(pop::Population)
    for i=1:length(pop)
        mutate(pop[i])
    end
end

# mutate a whole population with predefined std (not sure if anyone will ever need this)
function mutate(pop::Population, std::Float64)
    for i=1:length(pop)
        mutate(pop[i], std)
    end
end

function crossover(chr1::Chromosome, chr2::Chromosome, slices::Int64)
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

crossover(chr1::Chromosome, chr2::Chromosome) = crossover(chr1, chr2, 2)

################################################################################
## GENETIC ALGORITHM                                                          ##
################################################################################



