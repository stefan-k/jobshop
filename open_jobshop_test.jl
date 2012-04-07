# Test that creates the simplest of schedules
# for a random Open Job Shop probem

load("open_jobshop.jl")


num_machines = 5
num_jobs = 3

# Generate jobs:
jobs = Job[]

for i in 1:num_jobs
 	operations = Operation[]
 	# Add a random number of operations to the job:
 	for j = 1:num_machines
 	 	push(operations,Operation(i, j, randi(20)))
 	end
 	push(jobs, Job(i,operations))
end

# Create problem to be solved:
problem = OpenJobShopProblem(num_machines, jobs)

# Create initial schedule:
schedule = Schedule(problem)
print(schedule)
