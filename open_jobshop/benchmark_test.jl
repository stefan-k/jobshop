load("benchmark_generator.jl")
load("permutation_genetic.jl")
load("hybrid_genetic.jl")
load("selfish_gene.jl")


#
# This evaluation script attempts to display the same test cases as in the [Khuri] paper
#


# TODO Parameters:
# * Selfish: mut 0.1, reward 0.4, 10000 computations or 95% on every locus
# * Run 100 times -> mean!

# TEST CASE:
function benchmark_test()

    # Evaluation parameters:
    num_runs = 100     # = 100 in paper
    p_mutation  = 0.1 # = 0.1 in paper
    p_crossover = 0.6 # = 0.6 in paper

    
    # Parameters for the genetic algorithms (PGA):
    population_size = 200 # = 200 in paper
    num_generations = 500 # = 200 in paper


    selfish_iterations = 10000 # = 10 000 in paper
    selfish_reward = 0.04      # = 0.04 in paper
    selfish_stop = .95         # = .95 in paper

    benchmarks = [
                #max
                #op
        #size   #duration  #Time seed  #Machine seed
         4  4  99        1166510396  164000672
         4  4  99        1624514147  1076870026 
         4  4  99        1116611914  1729673136
         4  4  99        410579806   1453014524
         4  4  99        1036100146  375655500
         4  4  99        597897640   322140729
         4  4  99        1268670769  556009645
         4  4  99        307928077   421384574
         4  4  99        667545295   485515899
         4  4  99        35780816    492238933
        #15 15  99         840612802  398197754
        #15 15  99        1314640371  386720536
        #15 15  99        1227221349  316176388

    ]

    (num_problems, ~) = size(benchmarks)

    println()
    println("### Evaluation results for ", num_runs, " runs ###")
    println()

    if num_runs > 50 
        println("WARNING: This can take hours...")
        println()
    end

    println("| Problem |   Size | lower bound | Permutation GA  | Hybrid GA     | Selfish Gene   |") # TODO mean, other algorithms
    println("|         |        |             |   best |   mean | best |   mean |  best |   mean |")
    println("| ------- | ------ | ----------- | --------------- | ------------- | -------------- |")


    probs = GeneticProbabilities(p_mutation,p_crossover,0.0,0.0)

    for i = 1:num_problems
        b = benchmarks[i,:]
        problem = generate_problem(b[1],b[2],b[3],b[4],b[5])
        
        printf("|    %2i%2i | %2i x%2i | %11i |", b[1], i, b[1],b[2], lower_bound(problem))

        #1) Permutation Genetic:
        makespans = zeros(Int64, num_runs,1) # TODO: Is there a cheaper way to allocate an array?
        for j=1:num_runs
            schedule= permutation_genetic(problem, probs, population_size, num_generations)
            makespans[j] = compute_makespan(schedule)
        end
        printf(" %6i | %6.1f |", min(makespans), mean(makespans))


        # 2) Hybrid Genetic:
        makespans = zeros(Int64, num_runs,1) # TODO: Is there a cheaper way to allocate an array?
        for j=1:num_runs
            #probs = GeneticProbabilities(p_mutation,p_crossover,0.0,0.0)
            schedule= hybrid_genetic(problem, probs, population_size, num_generations)
            makespans[j] = compute_makespan(schedule)
        end
        printf(" %4i | %6.1f |", min(makespans), mean(makespans))


         # 3) Selfish Gene:
        makespans = zeros(Int64, num_runs,1) # TODO: Is there a cheaper way to allocate an array?
        for j=1:num_runs
            #probs = GeneticProbabilities(p_mutation,p_crossover,0.0,0.0)
            schedule = selfish_gene(problem, selfish_reward, selfish_stop, selfish_iterations) # TODO more parameters!
            makespans[j] = compute_makespan(schedule)
        end
        printf(" %5i | %6.1f |\n", min(makespans), mean(makespans))
    
    end


end

# Start test case:
@time benchmark_test()