# Testing the 1+1 evolutionary algorithm on the modified rosenbrock function
# from Magele's "Numerical Optimization" lecture.
# The global minimum is at [2.38, 5.32]

load("evolib.jl")

# objective function acting on the population
function rosenbrock(pop::Population)
    # thats a bit hacky for 1+1
    x1 = pop[1][1].gene
    x2 = pop[2][1].gene
    fitness = 100*(x2-x1^2)^2+(1-x1)^2-50*((x1+1)^2+(x2-1)^2)-10*((x1-1.5)^2+(x2-2.5)^2)+exp(2*x2-5)
end

# randomly initialize a population with 100 chromosomes
popul = rand(Population, 2, 1, 8.0, 6.0, -6.0 )

# launch genetic algorithm
@time evo_1plus1(popul, 0.000000000001, 0.85, 5.0, 1000, 100, rosenbrock)
