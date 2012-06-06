load("open_jobshop.jl")
load("../evolib.jl") 
load("benchmark_generator.jl")
load("hybrid_genetic.jl")
load("selfish_gene.jl")
load("permutation_genetic.jl")


#
# This evaluation script attempts to display the same test cases as in the [Khuri] paper
#


# TODO Parameters:
# * Selfish: mut 0.1, reward 0.4, 10000 computations or 95% on every locus
# * Run 100 times -> mean!

# TEST CASE:
function benchmark_test()

    # Evaluation parameters:
    #num_runs = 100     # = 100 in paper
    num_runs = 1     # = 100 in paper
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
4  4  99  1166510396  164000672
4  4  99  1624514147  1076870026 
4  4  99  1116611914  1729673136
4  4  99  410579806   1453014524
4  4  99  1036100146  375655500
4  4  99  597897640   322140729
4  4  99  1268670769  556009645
4  4  99  307928077   421384574
4  4  99  667545295   485515899
4  4  99  35780816    492238933
       
5  5  99  527556884       1343124817
5  5  99  1046824493      1973406531
5  5  99  1165033492      86711717
5  5  99  476292817      24463110 
5  5  99  1181363416     606981348
5  5  99  897739730      513119113
5  5  99  577107303      2046387124 
5  5  99  1714191910     1928475945 
5  5  99  1813128617     2091141708 
5  5  99  808919936     183753764

7 7 99  1840686215      1827454623
7 7 99  1026771938     1312166461
7 7 99  609471574     670843185
7 7 99  1022295947        398226875 
7 7 99  1513073047       1250759651 
7 7 99  1612211197      95606345
7 7 99  435024109      1118234860
7 7 99  1760865440        1099909092
7 7 99  122574075        10979313
7 7 99  248031774       1685251301 

15 15 99  1561423441     1787167667 
15 15 99  204120997     213027331 
15 15 99  801158374     1812110433
15 15 99  1502847623        1527847153 
15 15 99  282791231     1855451778 
15 15 99  1130361878        849417380
15 15 99  379464508        944419714 
15 15 99  1760142791       1955448160
15 15 99  1993140927      179408412
15 15 99  1678386613     1567160817

20 20 99  957638     9237185
20 20 99  162587311     1489531109 
20 20 99  965299017     1054695706
20 20 99  1158457671     1499999517
20 20 99  1191143707        1530757746
20 20 99  1826671743       901609771
20 20 99  1591533998      1146547719
20 20 99  937297777      92726463
20 20 99  687896268     1731298717 
20 20 99  687034842     684013066

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

    i2 = 0
    current = 0

    for i = 1:num_problems
        b = benchmarks[i,:]

        # size counter:
        if b[1] != current
            i2 = 0
            current = b[1]
        end
        i2 = i2 + 1

        problem = generate_problem(b[1],b[2],b[3],b[4],b[5])
        
        printf("|    %2i%2i | %2i x%2i | %11i |", b[1], i2, b[1],b[2], lower_bound(problem))

        #1) Permutation Genetic:
        makespans = zeros(Int64, num_runs,1) # TODO: Is there a cheaper way to allocate an array?
        tic()
        for j=1:num_runs
            schedule = permutation_genetic(problem, probs, population_size, num_generations)
            makespans[j] = compute_makespan(schedule)
        end
        toc()
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
