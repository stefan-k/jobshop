# Test that creates the simplest of schedules
# for a random Open Job Shop probem

load("open_jobshop.jl")

num_jobs = 5
num_machines = 10

# Create problem to be solved:
problem = rand(OpenJobShopProblem, num_jobs, num_machines)

# Create initial schedule:
schedule = Schedule(problem)
print(schedule)
