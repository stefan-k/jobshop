#
# The permutation genetic algorithm from [Khuri] plus helper functions 
#

load("open_jobshop.jl")
load("../evolib.jl")

#
# Prepare an OSSP problem for the genetic algorithm and run it
#
function permutation_genetic(problem::OpenJobShopProblem, probs::GeneticProbabilities, population_size, max_generations)

    objective_function = (x) -> ( makespan_objective_function(problem, x) )
    num_genes = count_operations(problem)

    population = rand(Population, population_size, num_genes, objective_function)
    result = genetic(population, probs, max_generations, objective_function)

    return schedule_from_chromosome(problem, result)

end

function permutation_genetic(problem::OpenJobShopProblem, population_size, max_generations)

    return schedule_from_chromosome(problem, GeneticProbabilities(1.0,1.0,1.0,1.0), population_size, max_generations)

end



#
# NOT USED
#
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
    op_map = Dict{Int64, Operation}(count_operations(problem))
    for job in problem.jobs
        for op in job.operations
            op_map[op.id] = op
        end
    end

    return op_map
end


#
# Compute the makespan for each chromosome and set it directly
#
function makespan_objective_function(problem::OpenJobShopProblem, population::Population)
    for i = 1:length(population)
        makespan_objective_function(problem, population.chromosomes[i])
    end
end

function makespan_objective_function(problem::OpenJobShopProblem, chromosome::Chromosome)
    chromosome.fitness = compute_makespan(schedule_from_chromosome(problem, chromosome))
end

#
# Function that combines the functions below
#
function schedule_from_chromosome(problem::OpenJobShopProblem, chromosome::Chromosome)
    return schedule_from_permutation_chromosome(problem, permutation_chromosome(chromosome))
end

# Much simpler and faster than the function below:
# just let the order of the genes decide the permutation,
# e.g. [2.5 7.3 2.2 1.6 9.2] -> [3 4 2 1 5]
function permutation_chromosome(c::Chromosome)    
    # Sort the genes:
    genes = map(x->x.gene, c.genes)
    (sorted, perm) = sortperm(genes)
    
    # Create the order from the sorting permutation
    num_genes = length(c)
    genes = Array(Gene, num_genes)
    for k=1:num_genes
        genes[perm[k]] = Gene(k)
    end

    return Chromosome(genes)
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
function OLD_permutation_chromosome(c::Chromosome)
    num_genes = length(c.genes)
    # Stretch chromosome over the range of values (creates a copy):
    chromosome = num_genes * c
    #genes = map(x->x.gene, chromosome.genes)
    values = [1:num_genes]

    for i = 1:num_genes
        gene = chromosome.genes[i].gene
        dists = map((x)->(abs(x - gene)), values)
        mini = find(dists == min(dists))
        index = mini[1]
        chromosome.genes[i].gene = index #ERROR
        values[index] = NaN # Make shure this one never gets chosen again
    end

    return chromosome
end

#
#  Create valid schedule from a permutation
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
