#--- Abstract evolutionary framework:
# TODO mutation
# TODO keep only children better than parents
load("evolib.jl")


# Simple evolutionary algorithm without mutation:
function evolution(initial_generation::Function,
                   replicator_function::Function,
                   fitness_function::Function,  
  	               stopping_criterion::Function)

	current_generation = copy(initial_generation) 
	generation_counter = 0

	print_generation(0,current_generation)
	
	while( ! stopping_criterion(current_generation, generation_counter))

		# Pair up randomly:
		shuffle!(current_generation) # in-place shuffle
		
		# Replicate:
		n = length(current_generation)
		for i = 1:n
			parent1 = current_generation[i]
			j = (i % n) + 1 # Can Julia also do arrays starting at 0? That would be nice!
			parent2 = current_generation[j]
			children = replicator_function(parent1, parent2)
			#println("Pair $(dec(parent1,10)) with $(dec(parent2,10)) -> $(map(dec,children,10))")
			current_generation = append(current_generation, children)
		end

		sort!(current_generation) #todo sort with fitness function
		current_generation = current_generation[1:n]
		generation_counter = generation_counter + 1
		print_generation(generation_counter, current_generation)
	end

	return current_generation[1]

end


#--- Test case: numbers

# Simplest fitness function: return argument itself
#function fitness(x::Number)
	#return x
#end
# easier way:
fitness(x::Number) = x

function print_generation(i, specimens)
	print("Generation $(dec(i,3)): ")
	println(map(hex,specimens,8))
end

# Simple replicator function:
function combine_integers(parent1::Uint32,parent2::Uint32)
	low_bits = randi(Uint32)
	high_bits = ~low_bits
	
	mutation_bit = 2 ^ randi((1,32))-1
	child1 = (parent1 & high_bits) | (parent2 & low_bits)
	child1 = child1 $ mutation_bit
	child1 = convert(Uint32, child1)

	mutation_bit = 2 ^ randi((1,32))-1
	child2 = (parent1 & low_bits ) | (parent2 & high_bits)
	child2 = child2 $ mutation_bit
	child2 = convert(Uint32, child2)

	return child1, child2 # return tuple instead of dense array
end

# Simple stopping criterion:
function stop_after_10(generation, counter)
	if counter > 100
		return true
	end

	return generation[1] < 5
end


#--- START

println()
println("Very simple example: shrink numbers until they approach zero:")
println()

initial_generation = [randi(Uint32) | i=1:5]
winner = evolution(initial_generation, combine_integers, fitness, stop_after_10)

println()
println("And the winner is: $winner")
println()
