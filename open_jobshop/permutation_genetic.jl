#
# The permutation genetic algorithm from [Khuri] plus helper functions 
#

#load("open_jobshop.jl")
#load("../evolib.jl")
#load("gaps.jl")
    

#
# Prepare an OSSP problem for the genetic algorithm and run it
#
function permutation_genetic(problem::OpenJobShopProblem, 
                             probs::GeneticProbabilities, population_size::Int, 
                             max_generations::Int)

    objective_function = (x::Union(Population, Chromosome)) -> ( makespan_objective_function(problem, x) )
    num_genes = count_operations(problem)

    population = rand(Population, population_size, num_genes, objective_function)
    result, convergence = genetic(population, probs, max_generations, objective_function)

    return schedule_from_chromosome(problem, result), convergence

end

function permutation_genetic(problem::OpenJobShopProblem, population_size::Int, max_generations::Int)

    return permutation_genetic(problem, GeneticProbabilities(1.0,1.0,1.0,1.0), population_size, max_generations)

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
        #makespan_objective_function(problem, population.chromosomes[i])
        makespan_objective_function(problem, population[i])
    end
end

function makespan_objective_function(problem::OpenJobShopProblem, chromosome::Chromosome)
    chromosome.fitness = compute_makespan(schedule_from_chromosome(problem, chromosome))
end

#
# Function that combines the functions below
#
function schedule_from_chromosome(problem::OpenJobShopProblem, 
                                  chromosome::Chromosome)
    return permutation_schedule_builder(problem, 
                                  permutation_chromosome(chromosome))
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

    # Initialize:
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
        start_time = int(max((machine_times[op.machine], job_times[op.job_index])))
        time_table[start_time] = op # Reference to operation!
        
        # Update both times for the next op:
        machine_times[op.machine] = start_time + op.duration
        job_times[op.job_index]   = start_time + op.duration
    end
    
    # Create Schedule object:
    return Schedule(time_tables)
end

#
# Try the same as below with gap lists:
#
function permutation_schedule_builder(problem::OpenJobShopProblem, chromosome)
    
    # Initilize:
    op_map = generate_op_map(problem)
    #time_tables = TimeTable[]
    #for i = 1:problem.num_machines
    #    push(time_tables, TimeTable())
    #end

    scheduler = OpenJobShopScheduler(problem)
    
    # Fill time tables:
    for k = 1:length(chromosome.genes)
        op_id = chromosome.genes[k].gene
        op = op_map[op_id]

        schedule_operation(scheduler, op)
    
    end

    # Create Schedule object:
    return create_schedule(scheduler)

end

#
#  Create valid schedule from a permutation
#
function NEW_schedule_from_permutation_chromosome(problem::OpenJobShopProblem, chromosome)
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
        start_time = int(max((machine_times[op.machine], job_times[op.job_index])))
        no_space = true 
        for j in 1:start_time # loop over all elements of current timetable to see if there is some free space
            # IDEA save sorted time table keys in vriable until it changes
            for k in sort!(keys(time_table)) # loop over all already existing elements of the timetable and ...
                k_dur = time_table[k].duration
                if k > (j + op.duration) # stop looping if k is out of the range
                    break
                end
                if (k + k_dur) < j # next loop if out of range
                    continue
                end
                if !((k >= j && k <= (j + op.duration)) || ((k + k_dur) >= j && (k + k_dur) <= (j + op.duration))) &&
                   !((j >= k && j <= (k + k_dur)) || ((j + op.duration) >= k && (j + op.duration) <= (k + k_dur))) # ... check if fitting the operation at the current position would be a valid operation
                    mach_no_space = false  # there obviously is some space, we just have to check if this collides ...
                    
                    #M = op.machines
                    #for l = [1:M-1, M+1:problem.num_machines] #exclude current machine

                    for l = 1:problem.num_machines # ... with tasks of the same job on other machines
                        if l == op.machine # do not check the machine where we want to place the task
                            continue 
                        end
                        machine_table = time_tables[l] # get timetable of machine we want to check
                        for m in sort!(keys(machine_table)) # check all slots for collissions
                            if op.job_index == machine_table[m].job_index # only check if it's the same job
                                m_dur = machine_table[m].duration
                                if m > j + op.duration # stopp loop if we are above slots where it makes sense to check
                                    break
                                end
                                if m + m_dur < j
                                    continue
                                end
                                if (m >= j && m <= (j + op.duration)) || ((m + m_dur) >= j && (m + m_dur) <= (j + op.duration)) ||
                                   (j >= m && j <= (m + m_dur)) || ((j + op.duration) >= m && (j + op.duration) <= (m + m_dur)) # check for collissions
                                    mach_no_space = true
                                    break # break if there is a collission
                                end
                            end
                        end
                        if mach_no_space == true
                            break
                        end
                    end
                    if mach_no_space == false # add if there is space
                        no_space = false
                        time_table[j] = op # Reference to operation!
                        machine_times[op.machine] = max(machine_times[op.machine], j + op.duration)
                        job_times[op.job_index] = max(job_times[op.job_index], j + op.duration)
                        break
                    end
                end
            end
            if no_space == false
                break
            end
        end
        if no_space == true  # append if it doesn't fit somewhere in 
            time_table[start_time] = op # Reference to operation!
            machine_times[op.machine] = start_time + op.duration
            job_times[op.job_index]   = start_time + op.duration
        end

        #time_table[start_time] = op # Reference to operation!
        
        # Update both times for the next op:
        #machine_times[op.machine] = start_time + op.duration
        #job_times[op.job_index]   = start_time + op.duration
    end
    
    # Create Schedule object:
    return Schedule(time_tables)
end
