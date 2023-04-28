using JuMP
using Ipopt
using PowerModels
using Base
using Random
@time begin

# include
include("./src/custom_opfs.jl")
include("./src/compute_limits.jl")
include("./src/compute_binding_constraints.jl")

"""
Main script to generate high-quality OPF datasets.

The routine has only been tested with the IEEE 30-, 57-, and 118-bus test systems. More systems will follow soon.
Untested systems may require additional changes in the source code to adapt to different characteristics, e.g. HVDC lines, non-active generators, ...

Inputs
---------------------------------------------------------------------------------------
case_name : string
    Name and directory of the case data.

epsilon_value: int, absolute
    Relaxation value for the complementary slackness binding_constraints_p

N: int, absolute
    Number of OPF solutions in the generated datasets

u_bound: int, pu
    The upper bound of the loads. It is calculated from the nominal loads.

l_bound: int, pu
    The lower bound of the loads. It is calculated from the nominal loads.

Outputs
----------------------------------------------------------------------------------------
opf_dataset_final : Dict
    Contains the matrices of the representative generated dataset by the simple OPF. They are divided into voltage, power and load.

opf_data : Dict
    Contains the matrices of the representative generated dataset by the bilevel optimization. The solutions are not guaranteed to be optimal. They are divided into voltage, power and load.

Verification statistics: Arrays
    Statistics that compute the average, maximum, and minimum percentage difference between loads, power generation and voltages.

Number of active sets: Arrays
    The number of times the constraint is active throughout the N iterations. The arrays are numerically ordered and calculated per element.

Visualization tools: Plots
    Plots that show the voltage, generation, flow and load profiles of the generated dataset.

"""

PowerModels.silence() 

# define the case
case_name = "pglib_opf_case30_ieee"  
case_file = case_name * ".m"
test_system_final = "./pglib_opf/" * case_file
case = PowerModels.parse_file(test_system_final)

# initialize
P_gen, Q_gen, Vm, Va, Pd, Qd, Pbranch, Qbranch = initialize_bilevel_solver(case)

# define relaxation of complementary slackness and iterations
epsilon_value = 1e-4
N             = 3

# define the bounds of the loads
u_bound = 1.5
l_bound = 0.5
bounds  = [u_bound,l_bound]

feasible_count = 0
feasible_count_opf = 0
infeasible_count_opf = 0

# run
opf_data, opf_dataset_final, init_soln = update_bilevel_result(case, P_gen, Q_gen, Vm, Va, Pd, Qd, Pbranch, Qbranch, N, epsilon_value, bounds)

num_bus = collect(1:length(case["bus"]))
num_load = collect(1:length(case["load"]))

println(" ")
println("********************** VERIFICATION STATISTICS ****************************")
println(" ")

# %%
average_precision_vm=[]
for (i,data) in opf_data[:Verification_vm]
    append!(average_precision_vm, data[2])
end
average_precision_value_vm = sum(average_precision_vm)/length(average_precision_vm)

println("Average value of voltage precision on all the iterations ",average_precision_value_vm, " %")
println("Maximum difference between voltages ", maximum(average_precision_vm), " %")
println("Minimum difference between voltages ", minimum(average_precision_vm), " %")

average_precision_generators = Dict()
for (i,data) in opf_data[:Verification_pg]
    for ii in data[2]
        if ii[1] in keys(average_precision_generators)
            append!(average_precision_generators[ii[1]], ii[2])
        else 
            average_precision_generators[ii[1]]= [ii[2]]
        end
    end
end

for i in collect(keys(average_precision_generators))
    average_precision_gen = sum(average_precision_generators[i])/length(average_precision_generators[i])
    println("Average value of generation precision at bus ", i, " on all the iterations ", average_precision_gen, " %")
    println("Maximum difference between generation at bus ", i, " ", maximum(abs.(average_precision_generators[i])), " %")
    println("Minimum difference between generation at bus ", i, " ", minimum(abs.(average_precision_generators[i])), " %")
end
   
ra, rb, rc, vmax_v, vmin_v, pgmax_v, qgmax_v = compute_limits(case)

binding_constraints_p, binding_constraints_vmax, 
binding_constraints_vmin, binding_constraints_pgmax, 
binding_constraints_qgmax = compute_binding_constraints(opf_dataset_final[:Pbranch_opf][:,2:end],opf_dataset_final[:Qbranch_opf][:,2:end],ra,
                                                                                                        opf_dataset_final[:Vm_opf][:,2:end], vmax_v, vmin_v,
                                                                                                        opf_dataset_final[:Pg_opf][:,2:end], opf_dataset_final[:Qg_opf][:,2:end], pgmax_v, qgmax_v)
                                                                                            

# show results
using Plots
gr(size=(1400,1400), legend=false)
p1=plot(opf_dataset_final[:Pg_opf][:,2:end], title="Active Power Generation")
p2=plot(opf_dataset_final[:Pd_opf][:,2:end],xaxis = (range(1, stop=length(num_load), step=1), range(1, stop=length(num_load), step=1)),
    yaxis = (range(-1, stop=1, step=0.2), range(-1, stop=1, step=0.2)), legend=false, title = "Active Power Demand")
p3 = plot(opf_dataset_final[:Vm_opf][:,2:end],xaxis = (range(1, stop=length(num_bus), step=1), range(1, stop=length(num_bus), step=1)), legend=false, title="Voltage Magnitude")
p4 = plot(opf_dataset_final[:Pbranch_opf][:,2:end], title="Branch Power Flow")
plot(p1, p2, p3, p4, layout=(4,1))

end
