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
# There is (or should be) a lot of intelligence in this function:
#  - Create valid schedule from a permutation
#  - Try to swap operations if it shortens the makespan
#
function schedule_from_chromosome(problem::OpenJobShopProblem, chromosome)

	# TODO handle non-integer, non-unique chromosomes received from the evolib

	op_map = generate_op_map(problem)
	
	# The order of operations within a job is arbitrary!
	
	# Naive (faulty) scheduling:
	# Put the tasks on the machines without checking if another job is running at the same time:

	# Initialize time tables:
	time_tables = TimeTable[]
	times = Int64[]
	for i = 1:problem.num_machines
		push(time_tables, TimeTable())
		push(times, 1)
	end

	# Fill time tables:
	for i = 1:length(chromosome.genes)
		op_id = chromosome.genes[i].gene
		op = op_map[op_id]
		time_table = time_tables[op.machine]
		time_table[times[op.machine]] = op # Reference to operation!
		#println("Machine ", op.machine, ": op ", op.id, ": ",times[op.machine]," + ", op.duration, " = ", times[op.machine]+op.duration)
		times[op.machine] += op.duration
	end
	

	return Schedule(time_tables)
end
