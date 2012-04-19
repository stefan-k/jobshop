# Testing the (mu/rho, lambda) evolutionary  algorithm on the modified rosenbrock function
# from Magele's "Numerical Optimization" lecture.
# The global minimum is at [2.38, 5.32]

load("evolib.jl")

# define the objective function (in-place!)
function rosenbrock(chr::Chromosome)
    x1 = chr[1].gene
    x2 = chr[2].gene
    chr.fitness = 100*(x2-x1^2)^2+(1-x1)^2-50*((x1+1)^2+(x2-1)^2)-10*((x1-1.5)^2+(x2-2.5)^2)+exp(2*x2-5)
end

# in case someone wants to calculate the objective function for a population
function rosenbrock(pop::Population)
    for i=1:length(pop)
        rosenbrock(pop[i])
    end
end

# randomly initialize a population with 100 chromosomes
popul = rand(Population, 100, 2, rosenbrock, 1.0, 6.0, -6.0 )

# launch genetic algorithm
@time evo_slash(popul, 2, 150, 1000, 0.85, 0.0000001, rosenbrock)
