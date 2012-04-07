#
# Every Job has to pass from machine 1 to machine N
# Every job takes a certain amount of time on each machine
#


# One job consists of multiple operations!
type Operation
	machine::Int64
	duration::Int64 # Everything's discrete
end

function print(op::Operation)
	print("Op$(dec(op.machine,2))($(lpad(op.duration,2)))")	
end


type Job
	operations::Array{Operation}
end

function print(job::Job)
	for op in job.operations
		print(" > ")
		print(op)

	end
end


type OpenJobShopProblem
	jobs::Array{Job}
end

function print(problem::OpenJobShopProblem)
	println("--- Open Job Shop Problem ---")
	for i = 1:length(problem.jobs)
		print("Job $i: ")
		println(problem.jobs[i])
	end
	print("-----------------------------")
end

# Simple schedule representation
# Can be invalid!
type Schedule
	start_times::HashTable
end	

# Generate an inefficient but valid schedule
function initial_schedule(shop::OpenJobShopProblem)
	
end


# TODO move this test to own file

# Generate jobs:
jobs = Job[]

for i in 1:3
 	operations = Operation[]
 	# Add a random number of operations to the job:
 	for j = 1:5
 	 	push(operations,Operation(j, randi(20)))
 	end
 	push(jobs, Job(operations))
end

# Create problem to be solved:
problem = OpenJobShopProblem(jobs)
println(problem)

