using Plots 
using Optim


function plot_matrix_optimization(points::Matrix{Float64})

    points = collect(transpose(points)) #makes a matrix where all x vals are in first row, y vals are in second row, z vals in third row --> collect() turns the ReshapedArray into a Matrix{Float64}
    points_flattened = vec(points)

    x_vals = points[1, :] #going over all collumns in the first row 
    y_vals = points[2, :]#going over all rows in the second column --> gathering all y's
    z_vals = points[3, :] #going over all rows in the third column --> gathering all z's

    energy_list = []
    iteration_list = []

    x0 = points_flattened #flatten the matrix by turning it into a vector so it can be aptly used in the optimization function

    #adding constraints by making a vector the same length as the #of coordinates so that each coordinate can be assigned a constraint
    lower_bound = fill(-10.0, length(x0)) 
    upper_bound = fill(10.0, length(x0))

    res = optimize(opt, lower_bound, upper_bound, x0, Fminbox(BFGS()), Optim.Options(store_trace = true, extended_trace = true, g_tol = 1e-4,    # looser gradient tolerance
    x_abstol = 1e-4,))

    #filling the 2 lists with content so its plottable 
    trace_steps = Optim.trace(res) #getting access to all of the step objects (iterations and energy vals) 
    for (i, step) in enumerate(trace_steps)
        push!(energy_list, step.value)
    end

    iteration_list = 1:length(trace_steps)
    println("Iteration Data: ", iteration_list)
    println("Energy Data: ", energy_list)
    println("Local Min energy value: ", res.minimum)
    println("Local Min atomic configuration: ", res.minimizer)
    
for i in 2:length(trace_steps)
    x1 = trace_steps[i-1].metadata["x"]
    x2 = trace_steps[i].metadata["x"]
    if x1 == x2
        println("Iteration $i is identical to Iteration $(i-1)")
    end
end

#-----------------

#PLOTTING PLOTTING 
    #Plotting the graph
    theme(:default)
    graph = plot(iteration_list, energy_list, color = "royalblue3", seriestype = :line, xlabel = "Iteration", ylabel = "Energy", background = "oldlace", label = "trendline")
    
    
    #Plotting the chosen points for visualization 
    midpoint = floor(Int, length(iteration_list)/2)
    random_nums = generateRandomNumsExcludeMidpoint(3, length(iteration_list)-1, midpoint) #will generate 3 random indeces for atomic coordinate representation (excluding the midpoint bc i want to make it guaranteed)
    #making the range length(iteration_list) -1 bc dont want last index to be generated
    rand1, rand2 = random_nums[1], random_nums[2] #getting the 2 random numbers from the list
    println("rand 1 is iteration #", rand1)
    println("rand 2 is iteration #", rand2)
    println("midpoint is #", midpoint)
    
    scatter!(
    [iteration_list[1], iteration_list[midpoint], iteration_list[rand1], iteration_list[rand2], iteration_list[length(iteration_list)]], #if you want to plot multiple scatter points, you have to put the values in nested lists 
    [energy_list[1], energy_list[midpoint], energy_list[rand1], energy_list[rand2], energy_list[length(iteration_list)]], 
    label = "chosen points", marker = (:star5, 6, :magenta3)) 
    
    display(graph) #displays the full graph with scatter points 
    
    #Seperate plots
    chosen_steps = [1, midpoint, rand1, rand2, length(iteration_list)]

    for c in chosen_steps #running through each step to get the coordinates for all atomic positions at that step and plot said coordinates 
        coordinate = trace_steps[c].metadata["x"] #gives the full vector of all atomix coordinates [x1, y1, z1, x2, y2, z2...] --> cant rly access x, y,z with this so need to reshape it!
        coordinate_reshaped = reshape(coordinate, 3, :) #turns into a 3 row infinite collumned vector [x, x, x; y, y, y; z, z, z]
        
        println("Iteration ", c, " coordintes: ", coordinate)
        println("Iteration ", c, " coordinates reshaped: ", coordinate_reshaped)

        x_v = coordinate_reshaped[1, :] #go through all collumns in 1st row
        y_v = coordinate_reshaped[2, :]
        z_v = coordinate_reshaped[3, :]

        if c == length(iteration_list)
            println("last iteration coordinates = ", coordinate)
            println("last iteration coordinates reshaped = ", coordinate_reshaped)
        end

        step_plot = scatter(x_v, y_v, z_v, label = "atoms", color = "magenta3", background_color = "seashell", title = "Iteration $c Configuration")
        display(step_plot)
        
    end

end
