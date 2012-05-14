
load("permutation_genetic.jl")


problem = rand(OpenJobShopProblem, 5, 2)

for i = 1:9 # Try some stuff out
    chromosome = 10*rand(Chromosome,10)
    chromosome = permutation_chromosome(chromosome)
    schedule = schedule_from_permutation_chromosome(problem, chromosome)
    println(compute_makespan(schedule))
end
