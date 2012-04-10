# Testing the genetic algorithm on the modified rosenbrock function
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
    @parallel for i=1:length(pop)
        rosenbrock(pop[i])
    end
end

# define the probability of each genetic operation
probs = GeneticProbabilities(0.799, 0.1, 0.1, 0.001)

# randomly initialize a population with 100 chromosomes
popul = rand(Population, 100, 2, rosenbrock, 1.0, 6.0, -6.0 )

# launch genetic algorithm
@time genetic(popul, probs, 1000, rosenbrock)
