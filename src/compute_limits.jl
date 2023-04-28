function compute_limits(case)

    branch_list = []
    rate_a = []
    rate_b = []
    rate_c = []
    bus_list = []
    vmax = []
    vmin = []
    gen_list = []
    max_p_gen = []
    max_q_gen = []

    for i in collect(keys(case["branch"]))
        append!(branch_list, convert(AbstractFloat, case["branch"][i]["index"]))
        append!(rate_a, convert(AbstractFloat, case["branch"][i]["rate_a"]))
        append!(rate_b, convert(AbstractFloat, case["branch"][i]["rate_b"]))
        append!(rate_c, convert(AbstractFloat, case["branch"][i]["rate_c"]))
    end

    for i in collect(keys(case["bus"]))
        append!(bus_list, convert(AbstractFloat, case["bus"][i]["bus_i"]))
        append!(vmax, convert(AbstractFloat, case["bus"][i]["vmax"]))
        append!(vmin, convert(AbstractFloat, case["bus"][i]["vmin"]))
    end

    for i in collect(keys(case["gen"]))
        append!(gen_list, convert(AbstractFloat, case["gen"][i]["gen_bus"]))
        append!(max_p_gen, convert(AbstractFloat, case["gen"][i]["pmax"]))
        append!(max_q_gen, convert(AbstractFloat, case["gen"][i]["qmax"]))
    end

    Rate_ainter = hcat(rate_a, branch_list)
    Rate_ainter1 = Rate_ainter[sortperm(Rate_ainter[:, 2]), :]
    Rate_a_def = Rate_ainter1[:, 1]

    Rate_binter = hcat(rate_b, branch_list)
    Rate_binter1 = Rate_binter[sortperm(Rate_binter[:, 2]), :]
    Rate_b_def = Rate_binter1[:, 1]

    Rate_cinter = hcat(rate_c, branch_list)
    Rate_cinter1 = Rate_cinter[sortperm(Rate_cinter[:, 2]), :]
    Rate_c_def = Rate_cinter1[:, 1]

    Vmaxinter = hcat(vmax, bus_list)
    Vmaxinter1 = Vmaxinter[sortperm(Vmaxinter[:, 2]), :]
    Vmax_def = Vmaxinter1[:, 1]

    Vmininter = hcat(vmin, bus_list)
    Vmininter1 = Vmininter[sortperm(Vmininter[:, 2]), :]
    Vmin_def   = Vmininter1[:, 1]

    Pgmaxinter = hcat(max_p_gen, gen_list)
    Pgmaxinter1 = Pgmaxinter[sortperm(Pgmaxinter[:, 2]), :]
    Pgmax_def   = Pgmaxinter1[:, 1]

    Qgmaxinter = hcat(max_q_gen, gen_list)
    Qgmaxinter1 = Qgmaxinter[sortperm(Qgmaxinter[:, 2]), :]
    Qgmax_def   = Qgmaxinter1[:, 1]

    return Rate_a_def, Rate_b_def, Rate_c_def, Vmax_def, Vmin_def, Pgmax_def, Qgmax_def

end

