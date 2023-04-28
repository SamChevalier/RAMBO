using JuMP
using Ipopt
using PowerModels
using Base
using Random
using Profile

# include
include("./src/custom_opfs.jl")
include("./src/compute_limits.jl")
include("./src/compute_binding_constraints.jl")
include("./src/generate_random_points.jl")

"""
Main script to generate high-quality OPF datasets.

The routine has only been tested with the IEEE 30-, 57-, and 118-bus test systems. More systems will follow soon.
Untested systems may require additional changes in the source code to adapt to different characteristics, e.g. HVDC lines, non-active generators, ...

Inputs
---------------------------------------------------------------------------------------
case_name : string
    Name and directory of the case data.

N: int, absolute
    Number of OPF solutions in the generated datasets

u_bound: int, pu
    The upper bound of the loads. It is calculated from the nominal loads.

l_bound: int, pu
    The lower bound of the loads. It is calculated from the nominal loads.

Outputs
----------------------------------------------------------------------------------------
random_dataset_final : Dict
    Contains the matrices of the representative generated dataset by the simple OPF using a uniform distribution. They are divided into voltage, power and load.

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

load_length = length(case["load"])

# initialize
P_gen, Q_gen, Vm, Va, Pd, Qd, Pbranch, Qbranch = initialize_bilevel_solver(case)

# define relaxation of complementary slackness and iterations
N             = 3

# define the bounds of the loads
u_bound = 1.5
l_bound = 0.5
bounds  = [u_bound,l_bound]

num_bus = collect(1:length(case["bus"]))
num_load = collect(1:length(case["load"]))

random_dataset_final, count_iter = random_sampling_opf(case, Pd, Qd,P_gen, Q_gen, Vm, Va, Pbranch, Qbranch , bounds, N, load_length)


ra, rb, rc, vmax_v, vmin_v, pgmax_v, qgmax_v = compute_limits(case)

binding_constraints_p, binding_constraints_vmax, 
binding_constraints_vmin, binding_constraints_pgmax, 
binding_constraints_qgmax = compute_binding_constraints(random_dataset_final[:Pbranch_random],random_dataset_final[:Qbranch_random],ra,
                                                                                                        random_dataset_final[:Vm_random], vmax_v, vmin_v,
                                                                                                        random_dataset_final[:Pg_random], random_dataset_final[:Qg_random], pgmax_v, qgmax_v)

# show results
using Plots
gr(size=(1400,1400), legend=false)
p1=plot(random_dataset_final[:Pg_random], title="Active Power Generation")
p2=plot(random_dataset_final[:Pd_random],xaxis = (range(1, stop=length(num_load), step=1), range(1, stop=length(num_load), step=1)),
yaxis = (range(-1, stop=1, step=0.2), range(-1, stop=1, step=0.2)), legend=false, title = "Active Power Demand")
p3 = plot(random_dataset_final[:Vm_random],xaxis = (range(1, stop=length(num_bus), step=1), range(1, stop=length(num_bus), step=1)), legend=false, title="Voltage Magnitude")
p4 = plot(random_dataset_final[:Pbranch_random], title="Branch Power Flow")
plot(p1, p2, p3, p4, layout=(4,1))

