include("add_opf_final_data.jl")
function random_sampling_opf(case, Pd_final, Qd_final,P_gen, Q_gen, Vm, Va, Pbranch, Qbranch , bounds, N, load_length)
    Pg_random = copy(P_gen)
    Qg_random = copy(Q_gen)
    Vm_random = copy(Vm)
    Va_random = copy(Va)
    Pd_random = copy(Pd_final)
    Qd_random = copy(Qd_final)
    Pbranch_random = copy(Pbranch)
    Qbranch_random = copy(Qbranch)

    global feasible_count_opf_random =0

    pd_max = Dict{Int, Float64}()
    pd_min = Dict{Int, Float64}()
    pd_sample = Dict{Int, Float64}()
    qd_max = Dict{Int, Float64}()
    qd_min = Dict{Int, Float64}()
    qd_sample = Dict{Int, Float64}()

    # parse bounds
    superior_bound=bounds[1]
    inferior_bound=bounds[2]

    # get the loading limits
    for ii in 1:load_length  # use load_length because Pd_final and Qd_final could have different sizes
        if Pd_final[ii, 1] > 0.0
            pd_max[ii] = Pd_final[ii]*superior_bound
            pd_min[ii] = Pd_final[ii]*inferior_bound
        elseif Pd_final[ii, 1] < 0.0
            pd_max[ii] = Pd_final[ii]*inferior_bound
            pd_min[ii] = Pd_final[ii]*superior_bound
        else
            pd_max[ii] = 0.0
            pd_min[ii] = 0.0
        end
        if Qd_final[ii,1] > 0.0
            qd_max[ii] = Qd_final[ii]*superior_bound
            qd_min[ii] = Qd_final[ii]*inferior_bound
        elseif Qd_final[ii, 1] < 0.0
            qd_max[ii] = Qd_final[ii]*inferior_bound
            qd_min[ii] = Qd_final[ii]*superior_bound
        else
            qd_max[ii] = 0.0
            qd_min[ii] = 0.0
        end
    end


    while feasible_count_opf_random < N
        pd_mean = 0
        qd_mean = 0
        pd_range = 0
        qd_range = 0
        for i in 1:load_length
            # active power
            pd_mean      = (pd_min[i] + pd_max[i])/2
            pd_range     = (pd_max[i] - pd_min[i])
            pd_sample[i] = pd_mean + pd_range*(rand(1)[1] - 0.5)


            # reactive power
            qd_mean      = (qd_min[i] + qd_max[i])/2
            qd_range     = (qd_max[i] - qd_min[i])
            qd_sample[i] = qd_mean + qd_range*(rand(1)[1] - 0.5)
        end

        for ii=1:load_length
            case["load"][string(ii)]["pd"]=pd_sample[ii]
            case["load"][string(ii)]["qd"]=qd_sample[ii]
        end

        T_max     = 150.0
        optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "max_cpu_time" => T_max)
        OPF_soln=PowerModels.solve_opf(case, ACPPowerModel, optimizer)
        if OPF_soln["termination_status"] in [LOCALLY_SOLVED, ALMOST_LOCALLY_SOLVED]
            feasible_count_opf_random += 1
            println("Number of solutions = ", feasible_count_opf_random)

            Pg_opf_v, Qg_opf_v, Vm_opf_v, Va_opf_v, Pd_opf_v, Qd_opf_v, Pbranch_opf_v, Qbranch_opf_v = create_final_opf_matrix(case, OPF_soln)

            Pg_random = hcat(Pg_random, Pg_opf_v)
            Qg_random = hcat(Qg_random, Qg_opf_v)
            Vm_random = hcat(Vm_random, Vm_opf_v)
            Va_random = hcat(Va_random, Va_opf_v)
            Pd_random = hcat(Pd_random, Pd_opf_v)
            Qd_random = hcat(Qd_random, Qd_opf_v)
            Pbranch_random = hcat(Pbranch_random, Pbranch_opf_v)
            Qbranch_random = hcat(Qbranch_random, Qbranch_opf_v)
        end
    end
    random_dataset_final = Dict(:Pg_random => Pg_random[:,2:end],
                    :Qg_random => Qg_random[:,2:end],
                    :Vm_random    => Vm_random[:,2:end],
                    :Va_random    => Va_random[:,2:end],
                    :Pd_random    => Pd_random[:,2:end],
                    :Qd_random    => Qd_random[:,2:end],
                    :Pbranch_random => Pbranch_random[:,2:end],
                    :Qbranch_random => Qbranch_random[:,2:end])
    return random_dataset_final, feasible_count_opf_random

end


