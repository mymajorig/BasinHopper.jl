using LinearAlgebra 

#lennard jones function 
function lennard_jones(distance, epsilon = 1, sigma = 1) #defaults, if u change them in calls parameters will just be overridden
      return 4 * epsilon * ((sigma / distance)^ 12 - (sigma / distance)^ 6)
end

#sum of potential energy between multiple pairs by adding up the lj vals of each distance 
function d_sum_lj_energy(distance) #takes in a list of distances 
    value  = 0
    for i in distance 
        value+=lennard_jones(i)
    end

    return value 
end 

#sum of potential energy given one distance val and one theta val 
function rtheta_sum_lj_energy(distance, theta)
    if abs(cos(theta)) < 1e-6
    error("theta too close to ±90 degrees, causing division by zero.")
    end

    atom1 = [0,0]
    atom2 = [r12, 0]

    midpoint_distance = r12/2
    midpoint = [midpoint_distance, 0]

   r13 = r12/(2*cos(theta))

   r23 = r13 

   return lennard_jones(r12) + lennard_jones(r13) + lennard_jones(r23)

end 

#calculate the distance between 2 atoms
function distance_two_atoms(atom1, atom2)
    return norm(atom2 - atom1)
end






