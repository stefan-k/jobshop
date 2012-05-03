#
# The permutation genetic algorithm from [Khuri] plus helper functions 
#

load("open_jobshop.jl")
load("../evolib.jl")

function permutation_genetic(problem::OpenJobShopProblem)

	chromosome = initial_chromosome(problem)

	print("Created initial chromosome: ")
	println(chromosome)
	schedule = schedule_from_permutation_chromosome(problem, chromosome)

	# TODO: # Create random initial population

	return schedule

end

function initial_chromosome(problem::OpenJobShopProblem)

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
# IDEA: The evolib creates real-valued genes, we need integers.
# Rounding is not enough though, because we need unique genes.
# This functions chooses the best matching gene for each operation and takes its
# index, so a valid permutation is created
#
# PROBLEM: It favors the left most genes over the right most genes
# IDEA: iterate over genes in sorted order
#
function permutation_chromosome(c::Chromosome)

	chromosome = copy(c)
	println("original chromo", chromosome)
	genes = map(x->x.gene, chromosome.genes)

	values = [1:1length(chromosome)]

	for i = 1:length(chromosome)
		gene = chromosome.genes[i].gene
		dists = map((x)->(abs(x - gene)), values)
		mini = find(dists == min(dists))
		index = mini[1]
		chromosome.genes[i].gene = index #ERROR
		values[index] = NaN # Make shure this one never gets chosen again
	end

	println("permutation chromo", chromosome)

	return chromosome
end


#
#  Create valid schedule from a permutation
# TODO handle non-integer, non-unique chromosomes received from the evolib
#
function schedule_from_permutation_chromosome(problem::OpenJobShopProblem, chromosome)
	# The order of operations within a job is arbitrary!
	# An operation can be scheduled when the machine is ready *and* the previous op
	# from the same job has finished

	# Initilize:
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
		
		# Update both times for the next op:
		machine_times[op.machine] = start_time + op.duration
		job_times[op.job_index]   = start_time + op.duration
	end
	
	# Create Schedule object:
	return Schedule(time_tables)
end
