#
# The permutation genetic algorithm from [Khuri] plus helper functions 
#

load("open_jobshop.jl")
load("../evolib.jl")

function permutation_genetic(problem::OpenJobShopProblem)

	chromosome = PermutationChromosome(problem)

	print("Created initial chromosome: ")
	println(chromosome)
	schedule = schedule_from_chromosome(problem, chromosome)

	# TODO: # Create random initial population

	return schedule

end

# TODO OpenJobShop namespace?
# type PermutationGeneticChromosome <: Chromosome # TODO This is more elegant. Why do i get a 

	
function PermutationChromosome(problem::OpenJobShopProblem)

	genes = Gene[]

	# For every operation, add its unique ID to the chromosome as a gene
	for job in problem.jobs
		for op in job.operations
			push(genes, Gene(op.id))
		end
	end

	#new(genes) # I hope this does what i want: Return a new PermutationGeneticChromosome with 'genes' as genes
	
	return Chromosome(genes)
end
	
# end


function generate_op_map(problem::OpenJobShopProblem)
	op_map = HashTable{Int64, Operation}()
	for job in problem.jobs
		for op in job.operations
			op_map[op.id] = op
		end
	end

	return op_map
end


#
#  Create valid schedule from a permutation
# TODO handle non-integer, non-unique chromosomes received from the evolib
#
function schedule_from_chromosome(problem::OpenJobShopProblem, chromosome)
	# The order of operations within a job is arbitrary!
	# An operation can be scheduled when the machine is ready *and* the previous op
	# from the same job has finished

	# Init
	op_map = generate_op_map(problem)
	time_tables = TimeTable[]
	machine_times = Int64[]
	for i = 1:problem.num_machines
		push(time_tables, TimeTable())
		push(machine_times, 1)
	end
	job_times = ones(length(problem.jobs))

	# Fill time tables:
	for i = 1:length(chromosome.genes)
		op_id = chromosome.genes[i].gene
		op = op_map[op_id]
		time_table = time_tables[op.machine]
		# Take first available time considering machine & job:
		start_time = max((machine_times[op.machine], job_times[op.job_index]))
		time_table[start_time] = op # Reference to operation!
		#println("Machine ", op.machine, ": op ", op.id, ": ",times[op.machine]," + ", op.duration, " = ", times[op.machine]+op.duration)
		
		# Update both times for the next op:
		machine_times[op.machine] = start_time + op.duration
		job_times[op.job_index]   = start_time + op.duration
	end
	

	return Schedule(time_tables)
end
