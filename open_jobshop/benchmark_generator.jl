load("open_jobshop.jl")

type RandomNumberGenerator
    seed::Int
end

function rand(generator::RandomNumberGenerator)

    seed = generator.seed


    a::Int32 = 16807
    b::Int32 = 127773
    c::Int32 = 2836
    m::Int32 = 2^31 - 1

    #println("type of seed/b ",typeof(seed/b), " ", seed/b)
    k = floor(seed/b)
    seed = a*(seed % b) - k*c
    if seed < 0
        println("seed $seed")
        seed += m
    end

    generator.seed = seed
    
    #println("$seed, $(seed/m)")

    return (seed/m)
end

# Get a random number from an interval:
rand(generator::RandomNumberGenerator, left, right) = convert(Int64, floor(left+rand(generator)*(right-left+1)))

function benchmark_generator(sizes, amounts)

    generator = RandomNumberGenerator(123) 
    rand(generator)



    n = length(sizes)
    @assert n == length(amounts)

    problems = OpenJobShopProblem[]
    for i=1:n
        for j=1:amounts[i]

            push(problems, rand(OpenJobShopProblem, sizes[i], sizes[i], 99))
        end
    end

    return problems
       

end


function generate_problem(num_jobs, num_machines, max_duration, time_seed, machine_seed)

    time_rand = RandomNumberGenerator(time_seed)
    machine_rand = RandomNumberGenerator(machine_seed)

    # Generate jobs:
    jobs = Job[]
    operation_id = 1
    for i in 1:num_jobs
        
        durations = [ rand(time_rand, 1, max_duration) for j=1:num_machines ]
        
        # Swapping:
        for j=1:num_machines
            temp = durations[j]
            j2 = rand(machine_rand, j, num_machines)
            durations[j ] = durations[j2]
            durations[j2] = temp
        end
        
        # Create Operation objects:
        operations = Operation[]
        for j=1:num_machines
            push(operations,Operation(i, j, durations[j], operation_id))
            operation_id += 1
        end

        push(jobs, Job(i,operations))
    end

    return  OpenJobShopProblem(num_machines, jobs)

end
