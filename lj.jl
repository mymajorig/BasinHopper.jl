

function lennard_jones(distance, epsilon = 1, sigma = 1) #defaults, if u change them in calls parameters will just be overridden
      return 4 * epsilon * ((sigma / distance)^ 12 - (sigma / distance)^ 6)
end


function sum_lj_energy(distance)
    value  = 0
    for i in distance 
        value+=lennard_jones(i)
    end

    return value 
end 