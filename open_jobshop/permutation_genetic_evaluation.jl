#
# Create 10 random problems, test them for 10, 100 and 1000 generations with
# different genetic probabilities
#

load("permutation_genetic.jl")

function permutation_genetic_evaluation()

    println()
    println("Average makespan reduction over 10 random problems...")
    println()
    println("| p(mut) | p(com) | p(rep) | p(imi) |    1 gen. |   10 gen. |  100 gen. | 1000 gen. |")
    println("| ------ | ------ | ------ | ------ | --------- | --------- |---------- | --------- |")


    # Create 10 random problems of size 5x9:
    num_jobs = 5
    num_machines = 9
    srand(231305786158076) # always create the same test cases, comment this out if you want a different test case in every run
    problems = [rand(OpenJobShopProblem, num_jobs, num_machines) for i=1:10]

    # Iterate over probabilities:
    range = 0.0:1.0:1.0
    population_size = 100

    for p_mut=range, p_com=range, p_rep=range, p_imi=range

        if p_mut == p_com == p_rep == p_imi == 0.0
            continue
        end
        
        printf("|    %.1f |    %.1f |    %.1f |    %.1f |",p_mut,p_com,p_rep,p_imi)
        probs = GeneticProbabilities(p_mut,p_com,p_rep,p_imi)

        # Iterate over number of generations:
        for num_generations in (1,10,100,1000)
            reduction = Number[]
            for problem in problems
                worst = Schedule(problem)
                best  = permutation_genetic(problem, probs, population_size, num_generations)
                push(reduction, compute_makespan(best)/ compute_makespan(worst))
            end
            printf(" %9.2f |", 100*mean(reduction))
        end

        println()

    end

end

@time permutation_genetic_evaluation()