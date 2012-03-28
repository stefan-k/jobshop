abstract AbstractGene
abstract AbstractChromosome
abstract AbstractPopulation
abstract AbstractGenerations

## GENE TYPE ##
type Gene{T<:Number} <: AbstractGene
    gene::T
    std::Float64
end

# Constructors
#Gene{T<:Number}(gene::T, std::Float64) = Gene(gene, std) # probably not necessary
Gene{T<:Number}(gene::T) = Gene(gene, 0.5)

# Utility functions
std(gene::Gene) = gene.std

## CHROMOSOME TYPE ##
type Chromosome <: AbstractChromosome
    genes::Vector
    length::Int64
end

# Constructor
Chromosome(genes) = Chromosome(genes, length(genes))

# referencing
ref(chromosome::Chromosome, ind...) = chromosome.genes[ind...]

# Utility functions
size(chromosome::Chromosome) = chromosome.length
length(chromosome::Chromosome) = chromosome.length
std(chromosome::Chromosome, ind...) = [std(gene) | gene = chromosome[ind...]]

## POPULATION TYPE ##
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

## GENERATIONS TYPE ##
# keeps track of generations
type Generations <: AbstractGenerations
    populations::Array
    generations::Int64

    function Generations(population)
        new([population], 1)
    end
end

# Utility functions

# Modifiers
function push(generations::Generations, population::Population)
    push(generations.populations, copy(population))
    generations.generations += 1
end
