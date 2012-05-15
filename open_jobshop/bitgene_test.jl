
load("permutation_genetic.jl")

srand(123456789)

problem = rand(OpenJobShopProblem, 5, 2)

genes = AbstractGene[]

# For every operation, add its unique ID to the chromosome as a gene
for job in problem.jobs
    for op in job.operations
        push(genes, BitGene(op.id))
    end
end

#new(genes) # I hope this does what i want: Return a new PermutationGeneticChromosome with 'genes' as genes

c1 = Chromosome(genes)
c2 = rand(Chromosome, BitGene, count_operations(problem))

println("Ordered chromosome:")
print(c1)
println()
println("Random chromosome:")
print(c2)
println()



#for i = 1:9 # Try some stuff out
#    chromosome = 10*rand(Chromosome,10)
#    chromosome = permutation_chromosome(chromosome)
#    schedule = schedule_from_permutation_chromosome(problem, chromosome)
#    println(compute_makespan(schedule))
#end
