#
# Generate a random problem and solve it using the permutation genetic algorithm 
#

load("permutation_genetic.jl")

# Initialize
num_jobs = 5
num_machines = 10
problem = rand(OpenJobShopProblem, num_jobs, num_machines)

# Solve
optimal_schedule = permutation_genetic(problem)

# Print
print(optimal_schedule)