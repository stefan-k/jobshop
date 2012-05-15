#
# Generate a random problem and solve it using the hybrid genetic algorithm 
#

load("hybrid_genetic.jl")


# Initialize
num_jobs = 5
num_machines = 9

srand(123) # always create the same test case, comment this out if you want a different test case in every run
problem = rand(OpenJobShopProblem, num_jobs, num_machines)

# Create initial schedule (just for comparison)
dumb_schedule = Schedule(problem)

# Solve
println()
println("Solving...")
population_size = 100
max_generations = 100
@time optimal_schedule = hybrid_genetic(problem, population_size, max_generations)


# Output
#println("Found schedule:")
print(optimal_schedule)
t1 = compute_makespan(dumb_schedule)
t2 = compute_makespan(optimal_schedule)
println()
print(num_jobs," jobs, ", num_machines, " machines, ", "Popsize=",population_size)
print(", Max generations=", max_generations)
print(", initial makespan: ", t1)
print(", optimal makespan: ", t2)
print(", reduced to: ")
printf("%.2f%%", (t2/t1)*100)
println()