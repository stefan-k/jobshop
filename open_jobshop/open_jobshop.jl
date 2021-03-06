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

# Useful for functions like max, sort:
isless(op1::Operation, op2::Operation) = isless(op1.duration, op2.duration)


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

rand(T::Type{OpenJobShopProblem}, num_jobs, num_machines) = rand(T::Type{OpenJobShopProblem}, num_jobs, num_machines, 20)

function rand(T::Type{OpenJobShopProblem}, num_jobs, num_machines, max_duration)
    # Generate jobs:
    jobs = Job[]
    operation_id = 1
    for i in 1:num_jobs
        operations = Operation[]
        # Add a random number of operations to the job:
        for j = 1:num_machines
            push(operations,Operation(i, j, randi(max_duration), operation_id))
            operation_id += 1
        end
        push(jobs, Job(i,operations))
    end

    return  OpenJobShopProblem(num_machines, jobs)
end

# Lower bound estimation from [Taillard'89]:
# "The lower bound of the makespans corresponds to the maximum amount of time that a job or a machine requires"
function lower_bound(problem::OpenJobShopProblem)
    
    job_max = max([ sum([op.duration for op in job.operations]) for job in problem.jobs ])
    
    machine_total_times = zeros(Int64, count_machines(problem))#[ 0 for i=1:count_machines(problem) ]
    for job in problem.jobs
        for j=1:count_machines(problem)# job.operations
            machine_total_times[j] += job.operations[j].duration
        end
    end
    machine_max = max(machine_total_times)
    
    return max( [job_max, machine_max] )
end

# Getter functions:
# TODO replace length(problem.jobs) by count_jobs(problem) everywhere
count_operations(problem::OpenJobShopProblem) = length(problem.jobs) * problem.num_machines
count_jobs(problem::OpenJobShopProblem) = length(problem.jobs)
count_machines(problem::OpenJobShopProblem) = problem.num_machines


################################################################################
## SCHEDULE TYPE                                                            ##
################################################################################

typealias TimeTable Dict{Int64, Operation}

# Simple schedule representation
# Can be invalid!
type Schedule
    
    # A time_table for each machine
    time_tables::Array{TimeTable}
    makespan::Number

    function Schedule(time_tables::Array{TimeTable}) # Why is this necessary?
        new(time_tables, Inf)
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
        new(time_tables, Inf)
    end
end 

# Return references to the operations that are next for a given time.
# e.g. if time==3 and two operations start at 5, both ops are returned.
function next_operations(schedule::Schedule, time::Int64)
    next_per_machine = map(mingeq, schedule.time_tables, time)
    next = find(next_per_machine .== min(next_per_machine))
    
    return (next, next_per_machine[next][1])
end

function compute_makespan(schedule::Schedule)

    if schedule.makespan != Inf
        return schedule.makespan
    end

    max_ops = map(max, schedule.time_tables)
    max_end_times = map((x)->(x[1]+x[2].duration-1), max_ops)
    schedule.makespan = max(max_end_times)

    return schedule.makespan
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
                print("  ")
                print(op)
                print(" |")
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


################################################################################
## SCHEDULER TYPE                                                             ##
################################################################################
# The scheduler keeps track of open gaps where operations fit, considering both
# machine and job restrictions.

load("gaps.jl")

type OpenJobShopScheduler # Name 'Scheduler' was taken

    num_jobs        ::  Int64
    num_machines    ::  Int64
    gap_lists       ::  Array{GapList,2}
    time_tables     ::  Vector{TimeTable}

    function OpenJobShopScheduler(problem::OpenJobShopProblem)

        num_jobs = count_jobs(problem)
        num_machines = count_machines(problem)

        time_tables = [ TimeTable() for i=1:problem.num_machines ]
        pseudo_inf = 999999999999
        gap_lists = [ [Gap(1,pseudo_inf)] for i=1:num_jobs,j=1:num_machines ]
    
        new(num_jobs, num_machines, gap_lists, time_tables)
    end

end

function schedule_operation(scheduler::OpenJobShopScheduler, op::Operation)
    time_table = scheduler.time_tables[op.machine]

    gap_list = scheduler.gap_lists[op.job_index, op.machine]
    start_time = find_gap(gap_list, op.duration)
    #println(start_time)

    @assert start_time > 0
    time_table[start_time] = op # Reference to operation!

    gap = Gap(start_time, op.duration)
    # The time slots this gap occupies have to be subtracted from the open gaps of job & machine:
    # 1) for the job:
    for i = 1:scheduler.num_jobs
        scheduler.gap_lists[i, op.machine] = subtract(scheduler.gap_lists[i, op.machine], gap)
    end

    for j = 1:scheduler.num_machines
        scheduler.gap_lists[op.job_index, j] = subtract(scheduler.gap_lists[op.job_index, j], gap)
    end
end

    create_schedule(scheduler::OpenJobShopScheduler) = Schedule(scheduler.time_tables)