#
# The hybrid genetic algorithm from [Khuri] plus helper functions 
#

load("open_jobshop.jl")
load("../evolib.jl")


#
# Prepare an OSSP problem for the hybrid genetic algorithm and run it
#
function hybrid_genetic(problem::OpenJobShopProblem, population_size, max_generations)

    objective_function = (x) -> ( hybrid_makespan(problem, x) )
    num_genes = count_operations(problem)

    population = rand(Population, population_size, num_genes, objective_function)
    result = genetic(population, GeneticProbabilities(1.0,1.0,1.0,1.0), max_generations, objective_function)

    return hybrid_schedule_builder(problem, result)

end

function hybrid_makespan(problem::OpenJobShopProblem, population::Population)
    for i = 1:length(population)
        hybrid_makespan(problem, population.chromosomes[i])
    end
end

function hybrid_makespan(problem::OpenJobShopProblem, chromosome::Chromosome)
    chromosome.fitness = compute_makespan(schedule_builder(problem, chromosome))
end

function hybrid_schedule_builder(problem::OpenJobShopProblem, chromosome::Chromosome)
   
    # Take real-valued genes and round them
    # (Better would be implementing BitGene)
    order = [ round(g.gene) for g in chromosome.genes ]
    
    num_jobs = count_jobs(problem)
    num_machines = count_machines(problem)
    
    unfinished_jobs = 1:num_jobs #TODO better take actual objects?
    unfinished_operations = [ 1:num_machines for i=1:num_jobs ]

    #TODO for index in order
    #   i = index % length(unfinished_jobs)
    #   ops = unfinished_operations[i]
    #   #Schedule longest? unfinished op j in ops
    #   #remove ops[j]
    #   #if ops is empty, remove i from unfinished_jobs
    # end


   return Schedule(problem) # dummy solution
end