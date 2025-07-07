 
include("potential_energy_functions.jl")

#TESTING LENNARD JONES
lennard_jones(1.1)
lennard_jones(2.4)

#TESTING sum_lj_energy that takes in a list of distances

distance_list = [10, 3, 4, 2]
sum_lj_energy(distance_list)

#testing sum_lj_energy with a distance and an angle
sum_lj_energy(2, 45)

#TESTING distance between 2 atoms
a1 = [2, 5]
a2 = [1,4]
a3 = [3,9]
a4 = [15, 90]

distance_two_atoms(a1, a2)
distance_two_atoms(a3, a4)

#testing rounding decimals 
f = 6.2385279
round_num(f, 5)
round_num(f,4)


