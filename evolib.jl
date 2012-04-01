################################################################################
## LIBRARY FOR EVOLUTIONARY ALGORITHMS                                        ##
################################################################################
#
# Authors: Joris Bayer
#          Stefan Kroboth
#
# TODO: 
#  * Create random Gene with constraints
#  * Create random Chromosome
#  * properly implement sort!() for the population
#  * try to get the type constraints to work
#  * more printing functions (human readable)
#  * possibility to flag chromosomes as 'bad'
#  * handling of objective functions
#  * roulette wheel selection 
#    - needs objective function and a working sort
#  * crossover
#  * mutation
#  * different selection/competition models 
#    - replace parents
#    - in competition with parents
#    - decision wether a child becomes part of the population or not
#    - in short: all that (kappa, lambda) stuff (Papers?)
#  * design similar to the one shown in prototype.jl
#    - although in a more general way
#  * stopping criterion
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
    # Point to discuss: Should 'factor' be a part of the Gene type?
    # Could be individually adapted, if necessary...
    function Gene(gene::Number) 
        new(gene, 0.5)
    end
    
    function Gene(gene::Number, std::Float64) 
        new(gene, std)
    end
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

################################################################################
## CHROMOSOME TYPE                                                            ##
################################################################################

type Chromosome <: AbstractChromosome
    genes::Vector
    length::Int64
    fitness::Float64
    obj_func::Function

    function Chromosome(genes::Vector{Gene}, obj_func::Function) 
        new(copy(genes), length(genes), Inf, obj_func)
        # one word to copy(): In my testing, I used just one gene several times
        # for the chromosome. If it is a reference, changing one gene changes 
        # (obviously) all of them. Therefore: copy()

        # the following might not work... 
        # maybe fitness shouldn't be calculated when created.
        # but I think it's a good idea to keep the objective function around.
        # Or maybe not (what if the objective function changes?). 
        # Comments requested ;)
        fitness = obj_func(genes)
    end
end

# referencing
ref(chromosome::Chromosome, ind...) = chromosome.genes[ind...]

# Utility functions
size(chromosome::Chromosome) = chromosome.length
length(chromosome::Chromosome) = chromosome.length
# print the stds of all genes of a chromosome indicated by a range/index
std(chromosome::Chromosome, ind...) = [std(gene) | gene = chromosome[ind...]]

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

################################################################################
## POPULATION TYPE                                                            ##
################################################################################

type Population <: AbstractPopulation
    chromosomes::Vector
    pop_size::Int64

    # several chromosomes passed
    function Population(chromosomes::Vector) 
        new(chromosomes, length(chromosomes))
    end

    # one chromosome passed
    function Population(chromosome::Chromosome) 
        new([chromosome], 1)
    end
end

# referencing
ref(population::Population, ind...) = population.chromosomes[ind...]

# Utility functions
size(population::Population) = population.pop_size
length(population::Population) = population.pop_size

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

function sort!(population::Population)
    # won't work this way! We will need the sort!(::Function,....) 
    # version which we have to pass a comparison function, so that
    # sorting is done based on fitness
    sort!(population.chromosomes)
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

# referencing
ref(generations::Generations, ind...) = generations.populations[ind...]

# Utility functions
get_generations(generations::Generations) = generations.generations

# Modifiers
function push(generations::Generations, population::Population)
    push(generations.populations, copy(population))
    generations.generations += 1
end
