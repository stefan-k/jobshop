
load("hybrid_genetic.jl")

srand(123456789)
problem = rand(OpenJobShopProblem, 3, 3)

for i = 1:9 # Try some stuff out
    chromosome = 9*rand(Chromosome,9)
    print(chromosome)
    println()
    schedule = hybrid_schedule_builder(problem, chromosome)
    println(compute_makespan(schedule))
end
