#
# Generate a random problem and solve it using the permutation genetic algorithm 
#

load("permutation_genetic.jl")

# Initialize
num_jobs = 5
num_machines = 9
problem = rand(OpenJobShopProblem, num_jobs, num_machines)

# Create initial schedule (just for comparison)
initial_schedule = Schedule(problem)

# Solve
optimal_schedule = permutation_genetic(problem)


# Output

println("Found schedule:")
println(optimal_schedule)

t1 = compute_makespan(initial_schedule)
t2 = compute_makespan(optimal_schedule)
println()
println("Initial makespan:  ", t1)
println("Optimal makespan:  ", t2)
println("Reduced to:       ", (t2/t1)*100,"%" )
println()