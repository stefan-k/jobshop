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
        println("limits", limits)
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
        println()
    end
    println()
end


#
#
#

function selfish_gene(problem::OpenJobShopProblem)

    schedule = Schedule(problem) # dummy solution

    num_jobs = count_jobs(problem)
    num_machines = count_machines(problem)

    sizes = Array(Int64, 2*num_machines)
    for i=1:num_machines
        sizes[2*(i-1)+1] = num_jobs
        sizes[2*(i-1)+2] = num_machines
    end
    vp = VirtualPopulation(sizes)
    print(vp)

    c = choose_chromosome(vp)
    print(c)

    return schedule

end





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