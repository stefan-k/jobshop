#load("open_jobshop.jl")

type RandomNumberGenerator
    seed::Int
end

function rand(generator::RandomNumberGenerator)

    seed = generator.seed

    a::Int32 = 16807
    b::Int32 = 127773
    c::Int32 = 2836
    m::Int32 = 2^31 - 1 # highest positive integer for Int32

    k = floor(seed/b)
    seed = a*(seed % b) - k*c
    if seed < 0
        seed += m
    end

    generator.seed = seed

    return (seed/m) # Integer/Integer division returns Float in Julia!
end

# Get a random number from an interval:
rand(generator::RandomNumberGenerator, left, right) = convert(Int64, floor(left+rand(generator)*(right-left+1))) #TODO


function generate_problem(num_jobs, num_machines, max_duration, time_seed, machine_seed)

    time_rand = RandomNumberGenerator(time_seed)
    machine_rand = RandomNumberGenerator(machine_seed)

    # Generate jobs:
    jobs = Job[]
    operation_id = 1
    for i in 1:num_jobs
        
        durations = [ rand(time_rand, 1, max_duration) for j=1:num_machines ]
        
        # Create indices (in non-open job shop scheduling, this would represent the machine order)
        index = [ idx for idx = 1:num_machines ]

        # Swapping:
        for j=1:num_machines
            temp = index[j]
            j2 = rand(machine_rand, j, num_machines)
            index[j ] = index[j2]
            index[j2] = temp
        end

        # Now, set durations[index[j]] to durations[j]:
        durations2 = zeros(Int64,1,num_machines)
        for j=1:num_machines
            durations2[index[j]] = durations[j]
        end
        
        # Create Operation objects:
        operations = Operation[]
        for j=1:num_machines
            push(operations,Operation(i, j, durations2[j], operation_id))
            operation_id += 1
        end

        push(jobs, Job(i,operations))
    end

    return  OpenJobShopProblem(num_machines, jobs)

end
