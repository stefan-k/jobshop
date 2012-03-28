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
Gene{T<:Number}(gene::T, std::Float64) = Gene(gene, std) # probably not necessary
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
    chromosomes::Array{T,1}
    size::Uint64
end

# Constructor
Population{T<:Chromosome}(chromosomes::Array{T,1}) = Chromosome(chromosomes, length(chromosomes))

# referencing
ref(population::Population, ind...) = population.chromosomes[ind...]

# Utility functions
size(population::Population) = population.size
length(population::Population) = population.size

## GENERATIONS TYPE ##
# keeps track of generations
type Generations{T<:Population} <: AbstractGenerations
    populations::Array{T}
    generations::Int64
end
