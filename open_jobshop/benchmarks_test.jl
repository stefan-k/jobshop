load("benchmark_generator.jl")

# TEST CASE:
function benchmark_test()

    problem = generate_problem(15,15, 99, 840612802, 398197754)
    println("Lower bound: ", lower_bound(problem))

    problem = generate_problem(15,15, 99, 1314640371, 386720536)
    println("Lower bound: ", lower_bound(problem))

    problem = generate_problem(15,15, 99, 1227221349, 316176388)
    println("Lower bound: ", lower_bound(problem))


    println("4x4 problems:")

    problem = generate_problem(4,4, 99, 1166510396,164000672)
    println("Lower bound: ", lower_bound(problem))

end

# Start test case:
benchmark_test()