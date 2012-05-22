load("../evolib.jl") 
load("open_jobshop.jl") 



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
steady_state(vp::VirtualPopulation) = ( mean_max(vp) > .95)

# Helper function for steady state
mean_max(vp::VirtualPopulation) = mean([ max(locus) for locus in vp.probabilities ])


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
function selfish_gene(problem::OpenJobShopProblem)

    num_jobs = count_jobs(problem)
    num_machines = count_machines(problem)

    sizes = Array(Int64, 2*num_jobs*num_machines)
    for i=1:num_jobs*num_machines
        sizes[2*(i-1)+1] = num_jobs
        sizes[2*(i-1)+2] = num_machines
    end
    
    vp = VirtualPopulation(sizes)
    max_iter = 5000
    best_sofar = choose_chromosome(vp)
    
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
        factor = 1.1
        reward(vp, winner, factor)
        punish(vp, loser , factor)

        # Check if the virtual pop. is steady:
        if steady_state(vp)
            println()
            println("  Population reached steady state at iteration ", i)
            break
        end

        # Compare the winner with previous best solution:      
        if fitness(winner) < fitness(best_sofar)
            printf("  New best at iteration %10i: %i\n", i, fitness(winner))
            best_sofar = winner
        end

    end

    println()

    return selfish_schedule_builder(problem, best_sofar)

end


function reward(vp::VirtualPopulation, chromosome::Chromosome, factor)
    values = [ convert(Int64, g.gene) for g in chromosome.genes ]
    @assert length(values) == size(vp)

    for i =1:length(values)
        value = values[i]
        vp.probabilities[i][value] *= factor
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
    job_times = ones(num_jobs)
    machine_times = ones(num_machines)
    time_tables = [ TimeTable(num_jobs) for i=1:num_machines ]

    for i =1:2:length(indices)

        job_index = ( (indices[i]  -1) % length(unfinished_jobs      )) + 1 # one-based indices
        unfinished_operations = unfinished_jobs[job_index]
        machine   = ( (indices[i+1]-1) % length(unfinished_operations)) + 1

        op = unfinished_operations[machine]

        # Take first available time considering machine & job:
        time_table = time_tables[op.machine]
        start_time = max((machine_times[op.machine], job_times[op.job_index]))
        time_table[start_time] = op # Reference to operation!
        
        # Update both times for the next op:
        machine_times[op.machine] = start_time + op.duration
        job_times[op.job_index]   = start_time + op.duration

        # Remove operation from unfinished:
        del(unfinished_operations, machine)
        if isempty(unfinished_operations)
            del(unfinished_jobs, job_index)
        end
    end

   return Schedule(time_tables)
end


################################################################################
## TEST CASE                                                                  ##
################################################################################
#
# TODO: move to test file
#
function main()
    # Initialize
    num_jobs = 5
    num_machines = 9

    srand(123) # always create the same test case, comment this out if you want a different test case in every run
    problem = rand(OpenJobShopProblem, num_jobs, num_machines)

    # Create initial schedule (just for comparison)
    dumb_schedule = Schedule(problem)

    # Solve
    println()
    println("Solving...")
    println()
    @time optimal_schedule = selfish_gene(problem)


    # Output
    #println("Found schedule:")
    #print(optimal_schedule)
    t1 = compute_makespan(dumb_schedule)
    t2 = compute_makespan(optimal_schedule)
    println()
    print(num_jobs," jobs, ", num_machines, " machines")
    print(", initial makespan: ", t1)
    print(", best found makespan: ", t2)
    print(", reduced to: ")
    printf("%.2f%%\n", (t2/t1)*100)
    println()


end

main()