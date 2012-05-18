
load("hybrid_genetic.jl")

srand(123456789)
problem = rand(OpenJobShopProblem, 5, 2)

for i = 1:9 # Try some stuff out
    chromosome = 10*rand(Chromosome,10)
    schedule = hybrid_schedule_builder(problem, chromosome)
    println(compute_makespan(schedule))
end
