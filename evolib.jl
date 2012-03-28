abstract AbstractChromosome
abstract AbstractPopulation
abstract AbstractGene

type Gene{T <: Number} <: AbstractGene
    gene::T
    std::Float64
end

# Constructors
Gene{T <: Number}(gene::T, std::Float64) = Gene(gene, std) # probably not necessary
Gene{T <: Number}(gene::T) = Gene(gene, 0.5)

type Chromosome{T <: Gene} <: AbstractChromosome
    genes::Array{T,1}
    length::Uint64
end

# Constructor
Chromosome{T <: Gene}(genes::Array{T,1}) = Chromosome(genes, length(genes))

type Population{T <: Chromosome} <: AbstractPopulation
    chromosomes::Array{T,1}
    size::Uint64
end

# Constructor
Population{T <: Chromosome}(chromosomes::Array{T,1}) = Chromosome(chromosomes, length(chromosomes))


