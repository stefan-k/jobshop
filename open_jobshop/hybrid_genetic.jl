#
# The hybrid genetic algorithm from [Khuri] plus helper functions 
#

load("open_jobshop.jl")
load("../evolib.jl")


#
# Prepare an OSSP problem for the hybrid genetic algorithm and run it
#
function hybrid_genetic(problem::OpenJobShopProblem, probs::GeneticProbabilities, population_size, max_generations)

    objective_function = (x) -> ( hybrid_makespan(problem, x) )
    num_genes = count_operations(problem)

    population = rand(Population, population_size, num_genes, objective_function)
    result = genetic(population, probs, max_generations, objective_function)

    return hybrid_schedule_builder(problem, result)

end

# Short constructor:
hybrid_genetic(problem::OpenJobShopProblem, population_size, max_generations) =
    hybrid_genetic(problem, GeneticProbabilities(1.0,1.0,1.0,1.0), population_size, max_generations)


function hybrid_makespan(problem::OpenJobShopProblem, population::Population)
    for i = 1:length(population)
        hybrid_makespan(problem, population.chromosomes[i])
    end
end

function hybrid_makespan(problem::OpenJobShopProblem, chromosome::Chromosome)
    chromosome.fitness = compute_makespan(hybrid_schedule_builder(problem, chromosome))
end

function hybrid_schedule_builder(problem::OpenJobShopProblem, chromosome::Chromosome)
   
    # Take real-valued genes and round them
    # (Better would be implementing BitGene)
    job_indices = [ convert(Int64, g.gene) for g in chromosome.genes ]
    
    num_jobs = count_jobs(problem)
    num_machines = count_machines(problem)
    
    # Work with sorted copies so the original arrays don't get messed with:
    unfinished_jobs = [ sortr(problem.jobs[i].operations) for i=1:num_jobs ]

    # Initialize timetables:
    job_times = ones(num_jobs)
    machine_times = ones(num_machines)
    time_tables = [ TimeTable(num_jobs) for i=1:num_machines ]

    for index in job_indices

        job_index = (index % length(unfinished_jobs)) + 1 # one-based indices
        unfinished_operations = unfinished_jobs[job_index]

        # print("job ",job_index," {")
        #     for o in unfinished_operations
        #         print(o)
        #         print(" ")
        #     end
        #     println("}")

        # Find longest op (operations are sorted descendantly):
        op = unfinished_operations[1]

        # Schedule longest:
        # Take first available time considering machine & job:
        time_table = time_tables[op.machine]
        start_time = max((machine_times[op.machine], job_times[op.job_index]))
        time_table[start_time] = op # Reference to operation!
        
        # Update both times for the next op:
        machine_times[op.machine] = start_time + op.duration
        job_times[op.job_index]   = start_time + op.duration

        # Remove operation from unfinished:
        del(unfinished_operations, 1)
        if isempty(unfinished_operations)
            del(unfinished_jobs, job_index)
        end
    end


   return Schedule(time_tables)
end