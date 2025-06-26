 
include("potential_energy_functions.jl")



#TESTING LENNARD JONES
lennard_jones(1.1)
lennard_jones(2.4)

#TESTING sum_lj_energy that takes in a list of distances

distance_list = [10, 3, 4, 2]
d_sum_lj_energy(distance_list)

#TESTING distance between 2 atoms
a1 = [2, 5]
a2 = [1,4]
a3 = [3,9]
a4 = [15, 90]

distance_two_atoms(a1, a2)
distance_two_atoms(a3, a4)