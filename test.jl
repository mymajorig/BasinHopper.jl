 
include("potential_energy_functions.jl")

#TESTING LENNARD JONES
lennard_jones(1.1)
lennard_jones(2.4)

#TESTING sum_lj_energy that takes in a list of distances

distance_list = [10, 3, 4, 2]
sum_lj_energy(distance_list)

#TESTING sum_lj_energy with a distance and an angle
sum_lj_energy(2, 45)

#TESTING lj energy but with a matrix that takes in coordinates of multiple atoms 
a1= [2,3,4]
a2 = [4,6,2]
a3 = [8,9,5]
a4 = [0,3,8]

sample = [a1'; a2'; a3'; a4'] # the semicolons make the computer interpret each a# as a ROW -> 4 rows and 1 column of atomic coordinates 
#sample2 = vcat(a1', a2', a3', a4') - alternate way of doing this 

#TESTING distance between 2 atoms
at1 = [2, 5]
at2 = [1,4]
at3 = [3,9]
at4 = [15, 90]

distance_two_atoms(at1, at2)
distance_two_atoms(at3, at4)

#testing rounding decimals 
f = 6.2385279
round_num(f, 5)
round_num(f,4)

#testing sum_lj_energy where the input is a matrix
atom1 = [4, 5, 30]
atom2 = [3, 12, 1]
atom3 = [13, 2, 18]
atom4 = [14, 5, 27]
positions = [atom1'; atom2'; atom3'; atom4']

