################################################################################
## TEST CASE                                                                  ##
################################################################################


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
    println()
    @time optimal_schedule = selfish_gene(problem)


    # Output
    #println("Found schedule:")
    print(optimal_schedule)
    t1 = compute_makespan(dumb_schedule)
    t2 = compute_makespan(optimal_schedule)
    println()
    print(num_jobs," jobs, ", num_machines, " machines")
    print(", lower bound: ", lower_bound(problem))
    #print(", initial makespan: ", t1)
    print(", best found makespan: ", t2)
    #print(", reduced to: ")
    #printf("%.2f%%\n", (t2/t1)*100)
    println()
    println()

end

main()
