load("../evolib.jl") 

# reimplementation of Chromosome because of the needed swappiness
type PermutationChromosome <: AbstractChromosome
    genes::Vector{Integer}
    swappiness::Integer # defines the number of Genes to swap
    length::Int64
    fitness::Float64

    function PermutationChromosome(jobs::Integer, machines::Integer, 
                                   swappiness::Integer)
        g = Integer[1:jobs*machines]
        swap_mutate!(g, jobs*machines) # scramble
        new(g, swappiness, jobs*machines, Inf)
    end

    function PermutationChromosome(genes::Vector{Integer}, swappiness::Integer,
                                   length::Int64, fitness::Float64)
        new(genes, swappiness, length, fitness)
    end
end

ref(chr::PermutationChromosome, ind...) = copy(chr.genes[ind...])
copy(chr::PermutationChromosome) = PermutationChromosome(copy(chr.genes),
                                                         copy(chr.swappiness),
                                                         copy(chr.length),
                                                         copy(chr.fitness))

function assign(chr::PermutationChromosome, gene::Integer, idx::Integer)
    chr.genes[idx] = copy(gene)
end

length(chr::PermutationChromosome) = chr.length

function isless(chr1::PermutationChromosome, chr2::PermutationChromosome)
    return chr1.fitness < chr2.fitness
end

#function sort!(pop::Vector{PermutationChromosome})
    #sort!(isless, pop)
#end

function swap_mutate!(genes::Vector{Integer}, swappiness::Integer)
    l = length(genes)
    for i in 1:swappiness
        idx = randi(l, 2)
        tmp = copy(genes[idx[1]])
        genes[idx[1]] = copy(genes[idx[2]])
        genes[idx[2]] = tmp
    end
    return genes
end

function swap_mutate!(chr::PermutationChromosome) 
    chr.genes = swap_mutate!(chr.genes, chr.swappiness)
    return chr
end

# TODO: Crossover!

function inv_fitness_sum(pop::Vector{PermutationChromosome}) 
    fitness_sum = 0
    for i=1:length(pop)
        fitness_sum += 1.0/pop[i].fitness
    end
    return fitness_sum
end

function roulette(pop::Vector{PermutationChromosome})
    f_sum = inv_fitness_sum(pop)
    idx = rand()*f_sum
    x = 0
    elem = 1
    for i=1:length(pop)
        x += abs(1.0/pop[i].fitness)
        if idx < x
            return elem
        end
        elem += 1
    end
    error("weird error that should not happen. You probably didn't define a fitness.")
end

type OSSP <: AbstractEvolutionary
    jobs::Array{Integer}
    machines::Array{Integer}
    duration::Array{Integer}
    num_jobs::Integer
    num_machines::Integer

    function OSSP(num_jobs::Integer, num_machines::Integer, 
                  duration::Vector{Int})
        jobs = Int[]
        machines = Int[]
        for i in 1:num_jobs
            jobs = vcat(jobs, repmat([i], num_machines, 1))
        end
        for i in 1:num_machines
            machines = vcat(machines, [1:num_machines])
        end
        new(jobs, machines, duration, num_jobs, num_machines)
    end
end

# lolwut
# TODO: Check if a operation can be scheduled inbetween two other operations
function makespan(chr::PermutationChromosome, p::OSSP)
    machine_count = zeros(Int, p.num_machines)
    job_count = zeros(Int, p.num_jobs)
    for i in 1:length(chr)
        idx = chr.genes[i]
        midx = p.machines[idx]
        jidx = p.jobs[idx]
        dur = p.duration[idx]
        job_count[jidx] = max(job_count[jidx], machine_count[midx]) + dur
        machine_count[midx] = max(job_count[jidx], machine_count[midx]) + dur
    end
    return max(machine_count)
end

# i <3 multiple dispatch
permutation_obj(chr::PermutationChromosome, p::OSSP) = chr.fitness = makespan(chr, p)

function permutation_obj(pop::Vector{PermutationChromosome}, p::OSSP)
    for i in 1:length(pop)
        pop[i].fitness = makespan(pop[i], p)
    end
end

function genetic(pop::Vector{PermutationChromosome}, p::OSSP, 
                 probabilities::GeneticProbabilities, iter::Int64, 
                 obj_func::Function)
    pop_o = copy(pop) # prevent in-place fiasco
    obj_func(pop_o, p) # make sure the the fitness for every chromosome is available
    
    best = copy(pop_o[1])
    best_generation = 0
    for j = 1:iter # max generations
        print("Generation ",dec(j,2), "... ")
        pop_n = PermutationChromosome[]
        for i = 1:length(pop_o)
            operation = roulette(probabilities)
            if operation == Mutation
                chr = swap_mutate!(pop_o[roulette(pop_o)])
            #elseif operation == Recombination
                #chr = crossover(pop_o)
                #chr = chr[1] # for future: take both offsprings!
            elseif operation == Reproduction
                chr = copy(pop_o[roulette(pop_o)])
            #elseif operation == Immigration
                #chr = rand(Chromosome, length(pop_o[1]), obj_func)
            end
            push(pop_n, chr)
        end
        obj_func(pop_n, p)
        sort!(pop_n)
        
        # Store best chromosome (don't know if the original algo does this as well)
        if pop_n[1].fitness < best.fitness
            print("New current best: $(pop_n[1].fitness)")
            best = copy(pop_n[1])
            best_generation = j
        end
        println()
        pop_o = copy(pop_n)
    end

    #println("Genetic algorithm found best solution w/ fitness $(best.fitness) in generation $best_generation.")

    return best
end

function main()
    jobs = 5
    machines = 9
    srand(123) # always create the same test case, comment this out if you want a different test case in every run
    #p = OSSP(jobs, machines, int([1:jobs*machines]))
    p = OSSP(jobs, machines, randi(20, jobs*machines)) # random times
    probs = GeneticProbabilities(0.30, 0.0, 0.70, 0.0)

    popu = [PermutationChromosome(jobs, machines, 2) for i = 1:1000]

    num_generations = 100

    best = genetic(popu, p, probs, num_generations, permutation_obj)
    println(best.genes)
    println(best.fitness)
end

@time main()
