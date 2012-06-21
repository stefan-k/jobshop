#load("../evolib.jl") 
#load("open_jobshop.jl") 



################################################################################
## VIRTUAL POPULATION TYPE                                                    ##
################################################################################
#
# A virtual population holds probabilities of values for certain positions in
# the chromosome. The actual values do not matter, they are represented by
# their index.
#
# locus  .. the index where a value is stored
# allele .. value stored at a certain locus    
#
type VirtualPopulation
    
    probabilities :: Vector{Vector{Float64}} # *Not* a matrix

    # Constructor with given probabilities
    function VirtualPopulation(probabilities::Vector{Vector{Float64}})
        new(probabilities)
    end

    # Constructor with equiprobable alleles 
    function VirtualPopulation(sizes::Vector{Int64}) # Initialize with equal probabilities
        probabilities = [ [1.0/size for i=1:size] for size in sizes ]
        new(probabilities)
    end

end

# Check if the average maximum of the marginal probabilities is close to 1.0,
# Which means the virtual population has "decided" on a certain allele for each locus
#
# The threshold was chosen empirically.
#
steady_state(vp::VirtualPopulation, stop) = ( minmax(vp) >= stop)

# Helper function for steady state
minmax(vp::VirtualPopulation) = min([ max(locus) for locus in vp.probabilities ])


#
# Randomly choose chromosome from virtual population
#
function choose_chromosome(vp::VirtualPopulation)
    n = size(vp)
    genes = Array(Gene, n)
    for i=1:n
        randval = rand()
        limits = cumsum(vp.probabilities[i])
        for j= 1:length(limits)
            if randval <= limits[j]
                genes[i] = Gene(j)
                break
            end
        end
    end

    return Chromosome(genes)
end

# Helper functions:
function copy(vp::VirtualPopulation)
    probabilities = [ [p for p in locus] for locus in vp.probabilities ]
    return VirtualPopulation(probabilities)
end


size(vp::VirtualPopulation) = length(vp.probabilities)

function print(vp::VirtualPopulation)
    println("Virtual Population:")
    for locus in vp.probabilities
        for prob in locus
            printf("  %.2f ",prob)
        end
        printf("(%.2f total)\n", sum(locus))
    end
    println()
end


################################################################################
## SELFISH GENE ALGORITHM                                                     ##
################################################################################
#
# Let two individuals compete by comparing their fitness. The gene pool of the
# winner gains dominance, i.e. its probabilities increase in the virtual
# population.
# The loser's probabilities are decreased.
#


function selfish_gene(problem::OpenJobShopProblem, reward_step, stop, max_iter)

    num_jobs = count_jobs(problem)
    num_machines = count_machines(problem)

    sizes = Array(Int64, 2*num_jobs*num_machines)
    for i=1:num_jobs*num_machines
        sizes[2*(i-1)+1] = num_jobs
        sizes[2*(i-1)+2] = num_machines
    end
    
    vp = VirtualPopulation(sizes)
    best_sofar = choose_chromosome(vp)

    convergence = zeros(Float, max_iter)
    
    for i = 1:max_iter

        # Pick two random chromosomes form virtual population:
        c1 = choose_chromosome(vp)
        c2 = choose_chromosome(vp)

        # Function shortcut:
        fitness(c::Chromosome) = selfish_fitness(problem,c)
        
        if fitness(c1) < fitness(c2)
            winner = c1
            loser  = c2
        else
            winner = c2
            loser  = c1
        end

        # Rewrite probabilities:
        reward(vp, winner, 1+reward_step)
        punish(vp, loser , 1+reward_step)

        # Check if the virtual pop. is steady:
        #println(minmax(vp))
        if steady_state(vp, stop)
            #println()
            #println("  Population reached steady state at iteration ", i)
            break
        end

        # Compare the winner with previous best solution:      
        if fitness(winner) < fitness(best_sofar)
            #printf("  New best at iteration %10i: %i\n", i, fitness(winner))
            best_sofar = winner
        end

        convergence[i] = fitness(best_sofar)

    end #loop

    #println()

    return selfish_schedule_builder(problem, best_sofar), convergence

end

# Shortcut:
selfish_gene(problem::OpenJobShopProblem) = selfish_gene(problem, .04,.95, 10000)



function reward(vp::VirtualPopulation, chromosome::Chromosome, factor)
    values = [ convert(Int64, g.gene) for g in chromosome.genes ]
    @assert length(values) == size(vp)

    for i =1:length(values)
        value = values[i]
        vp.probabilities[i][value] *= factor # We use a *multiplicative* reward, it yields much better results!
        vp.probabilities[i] /= sum(vp.probabilities[i])
    end

end

punish(vp::VirtualPopulation, chromosome::Chromosome, factor) = reward(vp::VirtualPopulation, chromosome::Chromosome, 1.0/factor)


function selfish_fitness(problem, chromosome)   

    if chromosome.fitness == Inf
        #println("Compute schedule...")
        schedule = selfish_schedule_builder(problem, chromosome)
        chromosome.fitness = compute_makespan(schedule)
    end

    return chromosome.fitness

end


################################################################################
## SELFISH GENE SCHEDULE BUILDER                                              ##
################################################################################
#
# Interpret a chromosome from the virtual population as a schedule. The
# chromosome must be of size (2 * #jobs * #machines).
#
# Example:
# A chromosome {1,2, 3,4, ...} is interpreted as first scheduling the
# 2nd unfinished operation of the 1st unfinished job, then the 4th unfinished
# operation from the 3rd unfinished job, etc.
# These indices are circular, i.e. if there are only 2 unfinished jobs, the
# index 3 refers to the first unfinished job.
#
function selfish_schedule_builder(problem::OpenJobShopProblem, chromosome::Chromosome)
   
    num_jobs = count_jobs(problem)
    num_machines = count_machines(problem)
    
    indices = [ convert(Int64, g.gene) for g in chromosome.genes ]
    @assert ( length(indices) == 2*num_jobs*num_machines) # must be even
      
    # Work with copies so the original arrays don't get messed with:
    unfinished_jobs = [ copy(problem.jobs[i].operations) for i=1:num_jobs ]

    # Initialize timetables:
    scheduler = OpenJobShopScheduler(problem)

    for i =1:2:length(indices)

        job_index = ( (indices[i]  -1) % length(unfinished_jobs      )) + 1 # one-based indices
        unfinished_operations = unfinished_jobs[job_index]
        machine   = ( (indices[i+1]-1) % length(unfinished_operations)) + 1

        op = unfinished_operations[machine]

        schedule_operation(scheduler, op)
        
        # Remove operation from unfinished:
        del(unfinished_operations, machine)
        if isempty(unfinished_operations)
            del(unfinished_jobs, job_index)
        end
    end

   return create_schedule(scheduler)
end
