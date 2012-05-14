#
# Every Job has to pass from machine 1 to machine N
# Every job takes a certain amount of time on each machine
#


################################################################################
## HELPER FUNCTIONS
################################################################################

#TODO move to own file?

#TODO built-in function exists?
# NOT USED
function minkey(table::Dict)
    mini = Inf
    for (k,v) in table
        if k < mini
            mini = k
        end
    end
    return mini
end

# Return entry with the smallest key which is greater or equal
# to mini
function mingeq(table::Dict, bottom::Int64)
    mini = typemax(Int64) # Inf doesnt work well with Int64 type
    for (k,v) in table
        
        if k < mini && k>= bottom
            mini = k
        end
        
    end
    
    return mini
end



################################################################################
## OPERATION TYPE                                                            ##
################################################################################

# One job consists of multiple operations!
type Operation
    job_index::Int64 # Pointer to job would make it circular
    machine::Int64
    duration::Int64 # Everything's discrete
    id::Int64 # Unique ID for each operation
end

function print(op::Operation)
    print("#$(dec(op.job_index,2))#$(dec(op.machine,2)):$(lpad(op.duration,2))")    
end


################################################################################
## JOB TYPE                                                            ##
################################################################################

type Job
    index::Int64
    operations::Array{Operation}
end

function print(job::Job)
    for op in job.operations
        print(" > ")
        print(op)

    end
end


################################################################################
## OPEN JOB SHOP TYPE                                                            ##
################################################################################

type OpenJobShopProblem
    num_machines::Int64
    jobs::Array{Job}
end

function print(problem::OpenJobShopProblem)
    println("--- Open Job Shop Problem ---")
    for i = 1:length(problem.jobs)
        print("Job $i: ")
        println(problem.jobs[i])
    end
    print("-----------------------------")
end

function rand(T::Type{OpenJobShopProblem}, num_jobs, num_machines)
    # Generate jobs:
    jobs = Job[]
    operation_id = 1
    for i in 1:num_jobs
        operations = Operation[]
        # Add a random number of operations to the job:
        for j = 1:num_machines
            push(operations,Operation(i, j, randi(20), operation_id))
            operation_id += 1
        end
        push(jobs, Job(i,operations))
    end

    return  OpenJobShopProblem(num_machines, jobs)
end

# Getter functions:
# TODO replace length(problem.jobs) by num_jobs(problem) everywhere
count_operations(problem::OpenJobShopProblem) = length(problem.jobs) * num_machines
count_jobs(problem::OpenJobShopProblem) = length(problem.jobs)

################################################################################
## SCHEDULE TYPE                                                            ##
################################################################################

typealias TimeTable Dict{Int64, Operation}

# Simple schedule representation
# Can be invalid!
type Schedule
    
    # A time_table for each machine
    time_tables::Array{TimeTable}

    function Schedule(time_tables::Array{TimeTable}) # Why is this necessary?
        new(time_tables)
    end

    # Create an initial valid schedule from an open job shop problem
    # by simply executing all jobs one after another
    function Schedule(problem::OpenJobShopProblem)
        
        # Initialize time tables:
        time_tables = TimeTable[]
        for i = 1:problem.num_machines
            push(time_tables, TimeTable(count_jobs(problem)))
        end

        # Fill time tables:
        time = 1
        for job in problem.jobs
            @assert length(job.operations) == problem.num_machines
            for op in job.operations
                time_table = time_tables[op.machine]
                time_table[time] = op # Reference to operation!
                time += op.duration
            end
        end
        new(time_tables)
    end
end 

# Return references to the operations that are next for a given time.
# e.g. if time==3 and two operations start at 5, both ops are returned.
function next_operations(schedule::Schedule, time::Int64)
    next_per_machine = map(mingeq, schedule.time_tables, time)
    next = find(next_per_machine == min(next_per_machine))
    
    return (next, next_per_machine[next][1])
end

function compute_makespan(schedule::Schedule)
    max_ops = map(max, schedule.time_tables)
    max_end_times = map((x)->(x[1]+x[2].duration-1), max_ops)
    
    return max(max_end_times)
end


# TODO: fix bug in makespan
function print(schedule::Schedule)
    num_machines = length(schedule.time_tables)

    # Print table header
    println()
    print("| TIME |")
    for i = 1:num_machines
        print(" Machine $(dec(i,2)) |")
    end
    println()
    
    time = 1
    current_start_times = (map((x)->(0), schedule.time_tables))
    #current_start_times = current_operations

    makespan = compute_makespan(schedule)
    
    while time <= makespan
        # TODO only get next_operations if necessary (not in every timestep)
        (machines, next_time) = next_operations(schedule, time)
        
        # Print timestep
        print("|", lpad(time,6), "|")
        
        # Print table cell for each machine:
        for i = 1:num_machines

            # Check if there's a new operation on this machine:
            if (time == next_time) && (contains(machines, i))
                current_start_times[i] = time
            end

            time_table = schedule.time_tables[i]
            current_start_time = current_start_times[i]
            
            # Check if the current operation has run out:
            if current_start_time > 0
                op = time_table[current_start_time]
                if (current_start_time + op.duration) <= time
                    current_start_time = 0
                end
            end
            
            # Finally, check if operation exists and print it:
            if current_start_time > 0
                print("  ", op, " |")
            else
                print("            |")
            end

        end
        println()
        
        time += 1
    end

    println()
    println("Makespan (total time): ", makespan, "s" )
    println()
end

