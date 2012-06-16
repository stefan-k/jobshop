#
# Generate a random problem and solve it using the permutation genetic algorithm 
#

load("open_jobshop.jl")
load("../evolib.jl")
load("permutation_genetic.jl")


# Initialize
num_jobs = 4
num_machines = 4
max_duration = 10

srand(123) # always create the same test case, comment this out if you want a different test case in every run
problem = rand(OpenJobShopProblem, num_jobs, num_machines, max_duration)

# Create initial schedule (just for comparison)
dumb_schedule = Schedule(problem)
initial_schedule = OLD_schedule_from_permutation_chromosome(problem, initial_chromosome(problem))

# Solve
println()
println("Solving...")
population_size = 200
max_generations = 200
@time optimal_schedule = permutation_genetic(problem, population_size, max_generations)


# Output
#println("Found schedule:")
print(optimal_schedule)
#t1 = compute_makespan(dumb_schedule)
t2 = compute_makespan(optimal_schedule)
println()
print(num_jobs," jobs, ", num_machines, " machines, ", "Popsize=",population_size)
print(", Num generations=", max_generations)
print(", lower bound: ", lower_bound(problem))
print(", best found makespan: ", t2)
#print(", reduced to: ")
#printf("%.2f%%", (t2/t1)*100)
println()
println()


# println("Compare with pure random number generator:")
# num_ops = num_jobs*num_machines
# function find_best_random(n)
#    chromosomes = [ rand(Chromosome, num_ops) for i=1:n ]
#    lessthan = (x,y)->(compute_makespan(schedule_from_chromosome(problem,x)) < compute_makespan(schedule_from_chromosome(problem,y)))
#    chromosomes = sort(lessthan, chromosomes)
#    return schedule_from_chromosome(problem, chromosomes[1])
# end
# printf("Creating max_generations x population_size (%i x %i) random chromosomes...\n", max_generations, population_size)
# @time  optimal_schedule = find_best_random(max_generations * population_size)
# t2 = compute_makespan(optimal_schedule)
# println()
# print(num_jobs," jobs, ", num_machines, " machines, ", "Popsize=",population_size)
# print(", Max generations=", max_generations)
# print(", initial makespan: ", t1)
# print(", optimal makespan: ", t2)
# print(", reduced to: ")
# printf("%.2f%%", (t2/t1)*100)
# println()
