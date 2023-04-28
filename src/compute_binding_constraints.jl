function compute_binding_constraints(Pbranch::Matrix{Float64}, Qbranch::Matrix{Float64}, rate_A::Vector{Float64}, Vm_finals::Matrix{Float64}, Vmaxs::Vector{Float64}, Vmins::Vector{Float64}, Pg_finals::Matrix{Float64}, Qg_finals::Matrix{Float64}, Pgmaxs::Vector{Float64}, Qgmaxs::Vector{Float64})

    num_lines = length(rate_A)
    binding_constraints_p = zeros(num_lines)

    for i = 1:size(Pbranch, 2)
        distances = (hypot.(Pbranch[:,i], Qbranch[:,i]) - rate_A)./rate_A

        for ii = 1:num_lines
            if abs(distances[ii]*100)< 1
                binding_constraints_p[ii] +=1
            end
        end
    end

    num_buses = length(Vmaxs)

    binding_constraints_vmax = zeros(num_buses)
    binding_constraints_vmin = zeros(num_buses)

    for i = 1:size(Vm_finals, 2)
        range_v = Vmaxs-Vmins
        rel_value_max = Vmaxs - Vm_finals[:,i]
        rel_value_min = Vm_finals[:,i] - Vmins

        distances_max = rel_value_max./range_v 
        distances_min = rel_value_min./range_v 

        for ii = 1:num_buses
            if abs(distances_max[ii]*100)< 1
                binding_constraints_vmax[ii] +=1
            end
            if abs(distances_min[ii]*100)< 1
                binding_constraints_vmin[ii] +=1
            end
        end
    end


    num_gens = length(Pgmaxs)

    binding_constraints_pgmax = zeros(num_gens)
    binding_constraints_qgmax = zeros(num_gens)

    for i = 1:size(Pg_finals, 2)
        distances_max_p = (Pg_finals[:,i] - Pgmaxs)./Pgmaxs
        distances_max_q = (Qg_finals[:,i] - Qgmaxs)./Qgmaxs

        for ii = 1:num_gens
            if abs(distances_max_p[ii]*100)< 1
                binding_constraints_pgmax[ii] +=1
            end
            if abs(distances_max_q[ii]*100)< 1
                binding_constraints_qgmax[ii] +=1
            end
        end
    end

    return binding_constraints_p, binding_constraints_vmax, binding_constraints_vmin, binding_constraints_pgmax, binding_constraints_qgmax
end