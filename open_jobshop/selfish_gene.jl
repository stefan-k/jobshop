load("../evolib.jl") 
load("open_jobshop.jl") 

#
# A virtual population holds probabilities of values for certain positions in the chromosome
# the actual values do not matter, they are represented by their index
#
type VirtualPopulation
    
    probabilities :: Array{Array{Float64}} # *Not* a matrix

    # locus  .. the index where a value is stored
    # allele .. value stored at a certain locus
    function VirtualPopulation(sizes::Array{Int64}) # Initialize with equal probabilities
        
        new([ [1.0/size for i=1:size] for size in sizes ])

        # probabilities = Array{Array{Float64}}(length(sizes))
        # for locus = 1:length(sizes)
        #     size = sizes[locus]
        #     p = 1.0 /size
        #     probabilities[locus] = [p for i=1:size]
        # end

        # new(probabilities)
    end

end

size(vp::VirtualPopulation) = length(vp.probabilities)

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


#
#
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

        #print(vp)

        #println("\nIteration ", dec(i,2),":")
        # Pick two random chromosomes form virtual population:
        c1 = choose_chromosome(vp)
        c2 = choose_chromosome(vp)
        #println("  C1: ", c1, "fitness: ", selfish_fitness(problem, c1))
        #println("  C2: ", c2, "fitness: ", selfish_fitness(problem, c2))
        
        #better_than = (a,b) -> (selfish_fitness(problem,a) < selfish_fitness(problem,b))
        #(s,p) = sort(better_than, [c1,c2])
        #println("  Order: ", p)

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
        reward(vp, loser , 1.0/factor)
        reward(vp, winner, factor)

        # Check if it's the best sofar        
        if fitness(winner) < fitness(best_sofar)
            printf("  New best at iteration %10i: %i\n", i, fitness(winner))
            best_sofar = winner
        end

    end

    return selfish_schedule_builder(problem, best_sofar)

end


function reward(vp::VirtualPopulation, chromosome::Chromosome, factor)
    values = [ convert(Int64, g.gene) for g in chromosome.genes ]
    @assert length(values) == length(vp.probabilities)

    for i =1:length(values)
        value = values[i]
        vp.probabilities[i][value] *= factor
        vp.probabilities[i] /= sum(vp.probabilities[i])
    end

end


function selfish_fitness(problem, chromosome)   

    if chromosome.fitness == Inf
        #println("Compute schedule...")
        schedule = selfish_schedule_builder(problem, chromosome)
        chromosome.fitness = compute_makespan(schedule)
    end

    return chromosome.fitness

end




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



# TODO: move to test file
function main()
    # Initialize
    num_jobs = 5
    num_machines = 9

    #srand(123) # always create the same test case, comment this out if you want a different test case in every run
    problem = rand(OpenJobShopProblem, num_jobs, num_machines)

    # Create initial schedule (just for comparison)
    dumb_schedule = Schedule(problem)

    # Solve
    println()
    println("Solving...")

    @time optimal_schedule = selfish_gene(problem)


    # Output
    #println("Found schedule:")
    #print(optimal_schedule)
    t1 = compute_makespan(dumb_schedule)
    t2 = compute_makespan(optimal_schedule)
    println()
    print(num_jobs," jobs, ", num_machines, " machines")
    print(", initial makespan: ", t1)
    print(", optimal makespan: ", t2)
    print(", reduced to: ")
    printf("%.2f%%", (t2/t1)*100)
    println()

end

main()