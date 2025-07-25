using LinearAlgebra 

#calculate the distance between 2 atoms
function distance_two_atoms(atom1, atom2)
    return norm(atom2 - atom1)
end

#for rounding numbers 
function round_num(num::Number, digits)

    places_moved = 10^digits
    num *= places_moved
    num = round(num) #round it
    num/= 10^digits #move it however many decimal places back so u can get the original version 

    return num 
end

function round_num(vector::AbstractVector, digits)
    for i in eachindex(vector) 
        vector[i] = round_num(vector[i], digits)
    end
    return vector
end 

#lennard jones function 
function lennard_jones(distance, epsilon = 1, sigma = 1) #defaults, if u change them in calls parameters will just be overridden
    if distance < 0.8
        distance = 0.8
    end
      return 4 * epsilon * ((sigma / distance)^ 12 - (sigma / distance)^ 6)
end

#--
#sum of potential energy between multiple pairs by adding up the lj vals of each distance 
function sum_lj_energy(distance) #takes in a list of distances 
    value  = 0
    for i in distance 
        value+=lennard_jones(i)
    end

    return value 
end 

#sum of potential energy given one distance val and one theta val 
function sum_lj_energy(distance, theta)
    if abs(cos(theta)) < 1e-6
    error("theta too close to ±90 degrees, causing division by zero.")
    end

    #for plotting purposes but not nessesary
    atom1 = [0,0]
    atom2 = [distance, 0]

    midpoint_distance = distance/2
    midpoint = [midpoint_distance, 0]

    #
   r13 = distance/(2*cos(theta)) #half of the distance aka half of r12 divided by cos theta

   r23 = r13 

   return lennard_jones(distance) + lennard_jones(r13) + lennard_jones(r23) 
end 

function sum_lj_energy_isosceles(distance, height)
    atom1 = [0, 0]
    atom2 = [distance, 0]
    midpoint = [distance/2, 0]

    atom3 = [distance/2, height]  # directly above midpoint, so isosceles

    r12 = distance
    r13 = sqrt((atom3[1] - atom1[1])^2 + (atom3[2] - atom1[2])^2)
    r23 = r13  # equal sides 

    return lennard_jones(r12) + lennard_jones(r13) + lennard_jones(r23)
end

function sum_lj_energy(distance, x3, y3) #by using coordinates or fixing one distance and positions, you automatically respect geometric constraints and ensure that the 3 distances in the situation end up forming a triangle
    #assumption that atom 1 is at (0,0)
    atom1 = [0,0]
    atom2 = [distance, 0.0]
    atom3 = [x3, y3]

    r13 = norm(atom3 .- atom1) #the dot just emphasizes that it subtracts each element in the list
    r23 = norm(atom3 .- atom2)

    return lennard_jones(distance) + lennard_jones(r13) + lennard_jones(r23)
end 

function sum_lj_energy(positions) #a positions matrix where columns = 3 and rows = however many atoms there are 
    #positions = reshape(positions, :, 3) # now N rows, 3 columns
    lj_vals = []
    total_pot = 0
    #println(typeof(positions)) #make sure that positions is a matrix 
    for i in 1:size(positions,1)-1 #(positions, 1) gives the amt of rows in the matrix (positions, 2 gives amt of columns)
        for j in (i+1):size(positions, 1)
            distance = distance_two_atoms(positions[i, :], positions[j, :]) #positions
            lj = lennard_jones(distance)
            push!(lj_vals, lj)
            total_pot += lj
        end
    end
    return total_pot
end

#--- END OF TOTAL POT FUNCTIONS

#turn the total potential energy calculation into a function so you can use the result for optimization 
function opt(x) #where x is a FLAT vector --> [1 2 3 4] rather than [1, 2, 3, 4]
    x = reshape(x, :, 3) #reshape x into a MATRIX so that it can be processed by sum_lj_energy-> can’t declare x as a matrix in opt’s argument because Optim.optimize always passes a flat vector regardless of how u try and set the parameter up
    energy_val = sum_lj_energy(x)
    return energy_val 
end

#pick a random index number in a list that excludes the midpoint index in that set 
function rand_exclude_midpoint(list_length, excluded_val)
        randomNum = rand(1:list_length)
        while(randomNum == excluded_val)
            randomNum = rand(1:list_length)
        end
        return randomNum 
    end 

function generateRandomNums(n, range) #generates a list of random numbers where none of them are the same
    #n = how many random numbers u need #range = the range of ur random numbers
    randomNumList = []

for i in 1:n
    randNum = rand(1:range)
    for i in randomNumList #runs through random num list 
        while randNum == i #while the random number generated is equal to 1 or more nums already in the list
            randNum = rand(1:range) #regenerate it
        end
    end
    println(randNum)
    push!(randomNumList, randNum) #push the number into the list 
end

end 

function generateRandomNumsExcludeMidpoint(n, range, excluded) #generates a list of random numbers where none of them are the same AND it excludes a certain index
    #n = how many random numbers u need #range = the range of ur random numbers
    randomNumList = []

for _ in 1:n #dont need 'i' 
    randNum = rand(1:range)
        while (randNum in randomNumList || randNum == excluded) #while the random number generated is equal to 1 or more nums already in the list
            randNum = rand(1:range) #regenerate it
        end
    println(randNum)
    push!(randomNumList, randNum) #push the number into the list 
end

end 
