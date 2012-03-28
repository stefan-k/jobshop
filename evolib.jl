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
type Chromosome{T<:Gene} <: AbstractChromosome
    genes::Array{T}
    length::Int64
end

# Constructor
Chromosome{T<:Gene}(genes::Array{T}) = Chromosome(genes, length(genes))

# referencing
ref(chromosome::Chromosome, ind...) = chromosome.genes[ind...]

# Utility functions
size(chromosome::Chromosome) = chromosome.length
length(chromosome::Chromosome) = chromosome.length
std(chromosome::Chromosome, ind...) = [std(gene) | gene = chromosome[ind...]]

## POPULATION TYPE ##
type Population{T<:Chromosome} <: AbstractPopulation
#type Population{T} <: AbstractPopulation
    chromosomes::Array
    pop_size::Int64
    function Population(chromosomes) 
        new(chromosomes, length(chromosomes))
    end
end

# Constructor
# Doesn't work, complains about input being of Type Array{Any,2} which shouldn't be the case
#Population{T<:Chromosome}(chromosomes::Array{T}) = Chromosome(chromosomes, length(chromosomes))
#Population{T}(chromosomes::Array{Any,2}) = Chromosome(chromosomes, length(chromosomes))
#Population(chromosomes) = Chromosome(chromosomes, length(chromosomes))

# referencing
ref(population::Population, ind...) = population.chromosomes[ind...]

# Utility functions
size(population::Population) = population.size
length(population::Population) = population.size

## GENERATIONS TYPE ##
# keeps track of generations
#type Generations{T<:Population} <: AbstractGenerations
type Generations{T} <: AbstractGenerations
    populations::Array{T}
    generations::Int64

    function push(population::T)
        append(populations, population)
        generations+=1
    end
end

# Constructor
# Doesn't work either
#Generations() = Generations([], 0)
#Generations{T}(population::T) = Generations([population], 1)
