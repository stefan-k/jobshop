abstract AbstractGene
abstract AbstractChromosome
abstract AbstractPopulation
abstract AbstractGenerations

################################################################################
## GENE TYPE                                                                  ##
################################################################################
type Gene{T<:Number} <: AbstractGene
    gene::T
    std::Float64
end

# Constructors
#Gene{T<:Number}(gene::T, std::Float64) = Gene(gene, std) # probably not necessary
Gene{T<:Number}(gene::T) = Gene(gene, 0.5)

# Utility functions
std(gene::Gene) = gene.std

################################################################################
## CHROMOSOME TYPE                                                            ##
################################################################################
type Chromosome <: AbstractChromosome
    genes::Vector
    length::Int64
    fitness::Float64
    obj_func::Function

    function Chromosome(genes::Gene, obj_func::Function) 
        new(genes, length(genes), Inf, obj_func)
        # might not work... 
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
std(chromosome::Chromosome, ind...) = [std(gene) | gene = chromosome[ind...]]

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
        new([population], 1)
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
