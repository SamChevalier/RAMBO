function opf_data_collection_bilevel(init_soln::Dict, data::Dict{String,Any},Pgen_final::Matrix{Float64},Qgen_final::Matrix{Float64},Voltage_mag_final::Matrix{Float64},Voltage_ang_final::Matrix{Float64},Pd_final::Matrix{Float64},Qd_final::Matrix{Float64},epsilon_value::Float64,bounds::Vector{Float64},opt_set,model=Model())

    # soften complementary slackness constraints
    eps_val = -epsilon_value

    # penalize complementarity constraints?
    penalize       = false
    penalty_weight = 25.0

    # initialize network
    @assert !haskey(data, "multinetwork")
    @assert !haskey(data, "conductors")

    # Add zeros to turn linear objective functions into quadratic ones
    # so that additional parameter checks are not required
    PowerModels.standardize_cost_terms!(data, order=2)

    # use build_ref to filter out inactive components
    ref = PowerModels.build_ref(data)[:it][pm_it_sym][:nw][nw_id_default]

    # ref contains all the relevant system parameters needed to build the OPF model
    # when we introduce constraints and variable bounds below, we use the parameters in ref.
    gen_index=collect(keys(ref[:gen]))
    bus_index=collect(keys(ref[:bus]))
    load_index=collect(keys(ref[:load]))
    branch_index=collect(keys(ref[:branch]))
    arcs_index=collect(keys(ref[:arcs]))
    ref_index=collect(keys(ref[:ref_buses]))

    # add voltage angles va for each bus
    @variable(model, va[i in keys(ref[:bus])])

    # Add voltage magnitudes vm for each bus
    @variable(model, ref[:bus][i]["vmin"] <= vm[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start=init_soln[:vm][i])

    # Add active power generation variable pg for each generator (including limits) + random start
    @variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"], start=init_soln[:pg][i])

    # Add reactive power generation variable qg for each generator (including limits)
    @variable(model, ref[:gen][i]["qmin"] <= qg[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"], start=init_soln[:qg][i])

    # Add power flow variables p to represent the active power flow for each branch
    @variable(model, -Inf <= p[(l,i,j) in ref[:arcs]] <= Inf, start=init_soln[:p][(l,i,j)])

    # Add power flow variables q to represent the reactive power flow for each branch
    @variable(model, -Inf <= q[(l,i,j) in ref[:arcs]] <= Inf, start=init_soln[:q][(l,i,j)])

    pd_max = Dict{Int, Float64}()
    pd_min = Dict{Int, Float64}()
    qd_max = Dict{Int, Float64}()
    qd_min = Dict{Int, Float64}()

    superior_bound=bounds[1]
    inferior_bound=bounds[2]

    for ii in 1:length(keys(ref[:load]))
        if Pd_final[ii, 1] > 0.0
            pd_max[ii] = Pd_final[ii, 1]*superior_bound
            pd_min[ii] = Pd_final[ii, 1]*inferior_bound
        elseif Pd_final[ii, 1] < 0.0
            pd_max[ii] = Pd_final[ii, 1]*inferior_bound
            pd_min[ii] = Pd_final[ii, 1]*superior_bound
        else
            pd_max[ii] = 0.0
            pd_min[ii] = 0.0
        end
        if Qd_final[ii,1] > 0.0
            qd_max[ii] = Qd_final[ii, 1]*superior_bound
            qd_min[ii] = Qd_final[ii, 1]*inferior_bound
        elseif Qd_final[ii, 1] < 0.0
            qd_max[ii] = Qd_final[ii, 1]*inferior_bound
            qd_min[ii] = Qd_final[ii, 1]*superior_bound
        else
            qd_max[ii] = 0.0
            qd_min[ii] = 0.0
        end
    end

    @variable(model, pd_min[i] <= pd[i in keys(ref[:load])] <= pd_max[i], start=init_soln[:pd][i])
    @variable(model, qd_min[i] <= qd[i in keys(ref[:load])] <= qd_max[i], start=init_soln[:qd][i])

    # generator limit complementarity constraints (mu1, mu2, mu3, mu4)
    @variable(model, 0 <= mu1[gen_index], base_name="mu1", start=0.0)
    @variable(model, 0 <= mu2[gen_index], base_name="mu2", start=0.0)
    @variable(model, 0 <= mu3[gen_index], base_name="mu3", start=0.0)
    @variable(model, 0 <= mu4[gen_index], base_name="mu4", start=0.0)

    if penalize == false
        for (i,gen) in ref[:gen]
            @NLconstraint(model, mu1[i]*(pg[i]-ref[:gen][i]["pmax"]) >= eps_val)
            @NLconstraint(model, mu2[i]*(ref[:gen][i]["pmin"]-pg[i]) >= eps_val)
            @NLconstraint(model, mu3[i]*(qg[i]-ref[:gen][i]["qmax"]) >= eps_val)
            @NLconstraint(model, mu4[i]*(ref[:gen][i]["qmin"]-qg[i]) >= eps_val)
        end
    end

    # branch phase angle limit complementarity constraints (mu5, mu6)
    @variable(model, 0 <= mu5[branch_index], base_name="mu5", start=0.0)
    @variable(model, 0 <= mu6[branch_index], base_name="mu6",  start=0.0)

    if penalize == false
        for (i,branch) in ref[:branch]
            va_fr = va[branch["f_bus"]]  
            va_to = va[branch["t_bus"]]  

            # Phase Angle Difference Limit
            @NLconstraint(model, mu5[i]*(branch["angmin"] - va_fr + va_to) >= eps_val)  
            @NLconstraint(model, mu6[i]*(va_fr - va_to - branch["angmax"]) >= eps_val)  
        end
    end

    # branch flow limit complementarity constraints (mu7, mu8)
    @variable(model, 0 <= mu7[branch_index], base_name="mu7", start=0.0)
    @variable(model, 0 <= mu8[branch_index], base_name="mu8", start=0.0)

    if penalize == false
        for (i,branch) in ref[:branch]
            f_idx = (i, branch["f_bus"], branch["t_bus"])
            t_idx = (i, branch["t_bus"], branch["f_bus"])
            
            if haskey(branch, "rate_a")  
                @NLconstraint(model, mu7[i]*(p[f_idx]^2 + q[f_idx]^2 - branch["rate_a"]^2) >= eps_val)
                @NLconstraint(model, mu8[i]*(p[t_idx]^2 + q[t_idx]^2 - branch["rate_a"]^2) >= eps_val)
            end
        end
    end

    # voltage magnitude limit complementarity constraints (mu9, mu10)
    @variable(model, 0 <= mu9[bus_index], base_name="mu9",   start=0.0)
    @variable(model, 0 <= mu10[bus_index], base_name="mu10", start=0.0)

    # set start values
    for i in keys(ref[:gen])
        set_start_value.(mu1[i], init_soln[:mu1][i])
        set_start_value.(mu2[i], init_soln[:mu2][i])
        set_start_value.(mu3[i], init_soln[:mu3][i])
        set_start_value.(mu4[i], init_soln[:mu4][i])

    end
    for (i,branch) in ref[:branch]
        set_start_value.(mu5[i], init_soln[:mu5][i])
        set_start_value.(mu6[i], init_soln[:mu6][i])
        set_start_value.(mu7[i], init_soln[:mu7][i])
        set_start_value.(mu8[i], init_soln[:mu8][i])

    end
    for i in keys(ref[:bus])
        set_start_value.(mu9[i],  init_soln[:mu9][i])
        set_start_value.(mu10[i], init_soln[:mu10][i])

    end

    if penalize == false
        for (i,bus) in ref[:bus]
            @NLconstraint(model,  mu9[i]*(vm[i]-ref[:bus][i]["vmax"]) >= eps_val)
            @NLconstraint(model, mu10[i]*(ref[:bus][i]["vmin"]-vm[i]) >= eps_val)
        end
    end

    # penalize
    if penalize == true
        # compute a penalized term
        @NLexpression(model, penalty_terms, sum(mu1[i]*(pg[i]-ref[:gen][i]["pmax"])+mu2[i]*(ref[:gen][i]["pmin"]-pg[i])+
                        mu3[i]*(qg[i]-ref[:gen][i]["qmax"])+mu4[i]*(ref[:gen][i]["qmin"]-qg[i]) for (i,gen) in ref[:gen] )+
                        sum(mu7[i]*(p[(i, branch["f_bus"], branch["t_bus"])]^2 + q[(i, branch["f_bus"], branch["t_bus"])]^2 - branch["rate_a"]^2) +
                        mu8[i]*(p[(i, branch["t_bus"], branch["f_bus"])]^2 + q[(i, branch["t_bus"], branch["f_bus"])]^2 - branch["rate_a"]^2) for (i,branch) in ref[:branch])+
                        sum(mu5[i]*(branch["angmin"] - va[branch["f_bus"]]   + va[branch["t_bus"]]  ) + 
                        mu6[i]*(va[branch["f_bus"]]   - va[branch["t_bus"]]   - branch["angmax"]) for (i,branch) in ref[:branch])+
                        sum(mu9[i]*(vm[i]-ref[:bus][i]["vmax"]) + mu10[i]*(ref[:bus][i]["vmin"]-vm[i]) for (i,bus) in ref[:bus]))
    end


    # equality constraints
    @variable(model, lambda1[i in keys(ref[:ref_buses])])
    @variable(model, lambda2[bus_index],    base_name="lambda2")
    @variable(model, lambda3[bus_index],    base_name="lambda3")
    @variable(model, lambda4[branch_index], base_name="lambda4")
    @variable(model, lambda5[branch_index], base_name="lambda5")    
    @variable(model, lambda6[branch_index], base_name="lambda6")    
    @variable(model, lambda7[branch_index], base_name="lambda7") 

    # reference bus
    for (i,bus) in ref[:bus]
        if isempty(ref[:bus_gens][i])
            # No action needed
        else
            terms_quadratic_cost=3 
            cost = ref[:gen][ref[:bus_gens][i][1]]["cost"]
            quadratic_cost::Array{Float64,1} = vcat(zeros(terms_quadratic_cost-length(cost)),cost)
            per_unit_powers=[2,1,0]
            quadratic_cost_opt::Array{Float64,1} = quadratic_cost ./(100 .^per_unit_powers)

            @NLconstraint(model,
                -lambda2[i] + quadratic_cost_opt[1]*pg[ref[:bus_gens][i][1]] + quadratic_cost_opt[2] + mu1[ref[:bus_gens][i][1]] 
                -mu2[ref[:bus_gens][i][1]] ==0)

        end
        
    end
    for (i,bus) in ref[:bus]
        if isempty(ref[:bus_gens][i])
            # No action needed
        else
            @NLconstraint(model,
                -lambda3[i]+ mu3[ref[:bus_gens][i][1]] -mu4[ref[:bus_gens][i][1]] ==0
            )
        end
    end
    branchesfinal=ref[:arcs][1:length(ref[:branch])]
    
    for (i,branch) in ref[:branch]

        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])
    
        @NLconstraint(model, lambda2[branch["f_bus"]]  + p[f_idx]*2*mu7[i] + lambda4[i] == 0)
        @NLconstraint(model, lambda3[branch["f_bus"]]  + q[f_idx]*2*mu7[i] + lambda6[i] == 0)
        @NLconstraint(model, lambda2[branch["t_bus"]]  + p[t_idx]*2*mu8[i] + lambda5[i] == 0)
        @NLconstraint(model, lambda3[branch["t_bus"]]  + q[t_idx]*2*mu8[i] + lambda7[i] == 0)
    end
    
    slack_bus = Dict{Int64, Vector{Int64}}()
    for (i, bus) in ref[:bus]
        slack_bus[i] = []
        if i in keys(ref[:ref_buses])
            slack_bus[i] = [i]
        end
    end

    global f_branch_index = []
    global f_bus_indexes_to = Dict()
    global f_bus_indexes = Dict()
    global f_g =    Dict()
    global f_b =    Dict()
    global f_tr =   Dict()
    global f_ti =   Dict()
    global f_g_fr = Dict()
    global f_b_fr = Dict()
    global f_g_to = Dict()
    global f_b_to = Dict()
    global f_tm =   Dict()

    global t_branch_index = []
    global t_bus_indexes_fr = Dict()
    global t_bus_indexes = Dict()
    global t_g =    Dict()
    global t_b =    Dict()
    global t_tr =   Dict()
    global t_ti =   Dict()
    global t_g_fr = Dict()
    global t_b_fr = Dict()
    global t_g_to = Dict()
    global t_b_to = Dict()
    global t_tm =   Dict()

        for (i,bus) in ref[:bus]
            bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]
            f_branch_index = []
            f_bus_indexes_to = Dict()
            f_bus_indexes = Dict()
            f_g =    Dict()
            f_b =    Dict()
            f_tr =   Dict()
            f_ti =   Dict()
            f_g_fr = Dict()
            f_b_fr = Dict()
            f_g_to = Dict()
            f_b_to = Dict()
            f_tm =   Dict()

            t_branch_index = []
            t_bus_indexes_fr = Dict()
            t_bus_indexes = Dict()
            t_g =    Dict()
            t_b =    Dict()
            t_tr =   Dict()
            t_ti =   Dict()
            t_g_fr = Dict()
            t_b_fr = Dict()
            t_g_to = Dict()
            t_b_to = Dict()
            t_tm =   Dict()

            for ii in branchesfinal
                if i ==ii[2]
                    append!(f_branch_index, ii[1])
                    f_bus_indexes[ii[1]] = ii[2]
                    f_bus_indexes_to[ii[1]] = ii[3]
                    g, b = PowerModels.calc_branch_y(ref[:branch][ii[1]])
                    tr, ti = PowerModels.calc_branch_t(ref[:branch][ii[1]])
                    f_g[ii[1]]    = g
                    f_b[ii[1]]    = b
                    f_tr[ii[1]]   = tr
                    f_ti[ii[1]]   = ti
                    f_g_fr[ii[1]] = ref[:branch][ii[1]]["g_fr"]
                    f_b_fr[ii[1]] = ref[:branch][ii[1]]["b_fr"]
                    f_g_to[ii[1]] = ref[:branch][ii[1]]["g_to"]
                    f_b_to[ii[1]] = ref[:branch][ii[1]]["b_to"]
                    f_tm[ii[1]] =   ref[:branch][ii[1]]["tap"]^2
                end
                if i ==ii[3]
                    append!(t_branch_index, ii[1])
                    t_bus_indexes[ii[1]] = ii[3]
                    t_bus_indexes_fr[ii[1]] = ii[2]
                    g, b = PowerModels.calc_branch_y(ref[:branch][ii[1]])
                    tr, ti = PowerModels.calc_branch_t(ref[:branch][ii[1]])
                    t_g[ii[1]]    = g
                    t_b[ii[1]]    = b
                    t_tr[ii[1]]   = tr
                    t_ti[ii[1]]   = ti
                    t_g_fr[ii[1]] = ref[:branch][ii[1]]["g_fr"]
                    t_b_fr[ii[1]] = ref[:branch][ii[1]]["b_fr"]
                    t_g_to[ii[1]] = ref[:branch][ii[1]]["g_to"]
                    t_b_to[ii[1]] = ref[:branch][ii[1]]["b_to"]
                    t_tm[ii[1]] =   ref[:branch][ii[1]]["tap"]^2
                end
            end

            
            @NLconstraint(model,sum(lambda1[k] for k in slack_bus[i]) +sum(mu5[a]-mu6[a] for a in t_branch_index)+ 
            sum(-mu5[s]+mu6[s] for s in f_branch_index)+
            sum(lambda4[l]*((-1)*(-t_g[l]*t_tr[l]+t_b[l]*t_ti[l])/t_tm[l]*(-vm[t_bus_indexes_fr[l]]*vm[t_bus_indexes[l]]*sin(va[t_bus_indexes_fr[l]]-va[t_bus_indexes[l]])*(-1)) -
            (-t_b[l]*t_tr[l]-t_g[l]*t_ti[l])/t_tm[l]*(vm[t_bus_indexes_fr[l]]*vm[t_bus_indexes[l]]*cos(va[t_bus_indexes_fr[l]]-va[t_bus_indexes[l]])*(-1))) for l in t_branch_index) +
            sum(lambda5[l]*((-1)*(-t_g[l]*t_tr[l]-t_b[l]*t_ti[l])/t_tm[l]*(-vm[t_bus_indexes[l]]*vm[t_bus_indexes_fr[l]]*sin(va[t_bus_indexes[l]]-va[t_bus_indexes_fr[l]])) -
            (-t_b[l]*t_tr[l]+t_g[l]*t_ti[l])/t_tm[l]*(vm[t_bus_indexes_fr[l]]*vm[t_bus_indexes[l]]*cos(va[t_bus_indexes[l]]-va[t_bus_indexes_fr[l]]))) for l in t_branch_index) +
            sum( lambda6[l]*((-t_b[l]*t_tr[l]-t_g[l]*t_ti[l])/t_tm[l]*(-vm[t_bus_indexes_fr[l]]*vm[t_bus_indexes[l]]*sin(va[t_bus_indexes_fr[l]]-va[t_bus_indexes[l]])*(-1)) -
            (-t_g[l]*t_tr[l]+t_b[l]*t_ti[l])/t_tm[l]*(vm[t_bus_indexes_fr[l]]*vm[t_bus_indexes[l]]*cos(va[t_bus_indexes_fr[l]]-va[t_bus_indexes[l]])*(-1))) for l in t_branch_index) +
            sum( lambda7[l]*((-t_b[l]*t_tr[l]+t_g[l]*t_ti[l])/t_tm[l]*(vm[t_bus_indexes[l]]*vm[t_bus_indexes_fr[l]]*sin(va[t_bus_indexes_fr[l]]-va[t_bus_indexes[l]])) -
            (-t_g[l]*t_tr[l]-t_b[l]*t_ti[l])/t_tm[l]*(vm[t_bus_indexes_fr[l]]*vm[t_bus_indexes[l]]*cos(va[t_bus_indexes[l]]-va[t_bus_indexes_fr[l]]))) for l in t_branch_index) +
            sum(lambda4[l]*((-1)*(-f_g[l]*f_tr[l]+f_b[l]*f_ti[l])/f_tm[l]*(-vm[f_bus_indexes[l]]*vm[f_bus_indexes_to[l]]*sin(va[f_bus_indexes[l]]-va[f_bus_indexes_to[l]]))-
            (-f_b[l]*f_tr[l]-f_g[l]*f_ti[l])/f_tm[l]*(vm[f_bus_indexes[l]]*vm[f_bus_indexes_to[l]]*cos(va[f_bus_indexes[l]]-va[f_bus_indexes_to[l]]))) for l in f_branch_index) +
            sum(lambda5[l]*((-1)*(-f_g[l]*f_tr[l]-f_b[l]*f_ti[l])/f_tm[l]*(-vm[f_bus_indexes_to[l]]*vm[f_bus_indexes[l]]*sin(va[f_bus_indexes_to[l]]-va[f_bus_indexes[l]])*(-1))-
            (-f_b[l]*f_tr[l]+f_g[l]*f_ti[l])/f_tm[l]*(vm[f_bus_indexes_to[l]]*vm[f_bus_indexes[l]]*cos(va[f_bus_indexes_to[l]]-va[f_bus_indexes[l]])*(-1))) for l in f_branch_index) +
            sum(lambda6[l]*((-f_b[l]*f_tr[l]-f_g[l]*f_ti[l])/f_tm[l]*(-vm[f_bus_indexes[l]]*vm[f_bus_indexes_to[l]]*sin(va[f_bus_indexes[l]]-va[f_bus_indexes_to[l]]))-
            (-f_g[l]*f_tr[l]+f_b[l]*f_ti[l])/f_tm[l]*(vm[f_bus_indexes[l]]*vm[f_bus_indexes_to[l]]*cos(va[f_bus_indexes[l]]-va[f_bus_indexes_to[l]]))) for l in f_branch_index) +
            sum(lambda7[l]*((-f_b[l]*f_tr[l]+f_g[l]*f_ti[l])/f_tm[l]*(-vm[f_bus_indexes_to[l]]*vm[f_bus_indexes[l]]*sin(va[f_bus_indexes[l]]-va[f_bus_indexes_to[l]]))-
            (-f_g[l]*f_tr[l]-f_b[l]*f_ti[l])/f_tm[l]*(vm[f_bus_indexes_to[l]]*vm[f_bus_indexes[l]]*cos(va[f_bus_indexes_to[l]]-va[f_bus_indexes[l]])*(-1))) for l in f_branch_index) == 0)
            
            @NLconstraint(model, sum(shunt["gs"] for shunt in bus_shunts)*vm[i]*2*lambda2[i] - sum(shunt["bs"] for shunt in bus_shunts)*vm[i]*2*lambda3[i]+mu9[i]-mu10[i]+
                sum(lambda4[l]*(-(-t_g[l]*t_tr[l]+t_b[l]*t_ti[l])/t_tm[l]*(vm[t_bus_indexes_fr[l]]*cos(va[t_bus_indexes_fr[l]]-va[t_bus_indexes[l]])) -
                (-t_b[l]*t_tr[l]-t_g[l]*t_ti[l])/t_tm[l]*(vm[t_bus_indexes_fr[l]]*sin(va[t_bus_indexes_fr[l]]-va[t_bus_indexes[l]]))) for l in t_branch_index) +
                sum(lambda5[l]*(-(t_g[l]+t_g_to[l])*vm[t_bus_indexes[l]]*2 -
                (-t_g[l]*t_tr[l]-t_b[l]*t_ti[l])/t_tm[l]*(vm[t_bus_indexes_fr[l]]*cos(va[t_bus_indexes[l]]-va[t_bus_indexes_fr[l]])) -
                (-t_b[l]*t_tr[l]+t_g[l]*t_ti[l])/t_tm[l]*(vm[t_bus_indexes_fr[l]]*sin(va[t_bus_indexes[l]]-va[t_bus_indexes_fr[l]])) ) for l in t_branch_index) +
                sum(lambda6[l]*((-t_b[l]*t_tr[l]-t_g[l]*t_ti[l])/t_tm[l]*(vm[t_bus_indexes_fr[l]]*cos(va[t_bus_indexes_fr[l]]-va[t_bus_indexes[l]])) -
                (-t_g[l]*t_tr[l]+t_b[l]*t_ti[l])/t_tm[l]*(vm[t_bus_indexes_fr[l]]*sin(va[t_bus_indexes_fr[l]]-va[t_bus_indexes[l]]))) for l in t_branch_index) +
                sum(lambda7[l]*((t_b[l]+t_b_to[l])*vm[t_bus_indexes[l]]*2 +
                (-t_b[l]*t_tr[l]+t_g[l]*t_ti[l])/t_tm[l]*(vm[t_bus_indexes_fr[l]]*cos(va[t_bus_indexes_fr[l]]-va[t_bus_indexes[l]])) -
                (-t_g[l]*t_tr[l]-t_b[l]*t_ti[l])/t_tm[l]*(vm[t_bus_indexes_fr[l]]*sin(va[t_bus_indexes[l]]-va[t_bus_indexes_fr[l]])) ) for l in t_branch_index) +
                sum(lambda4[j]*(-(f_g[j]+f_g_fr[j])/f_tm[j]*vm[f_bus_indexes[j]]*2 -
                (-f_g[j]*f_tr[j]+f_b[j]*f_ti[j])/f_tm[j]*(vm[f_bus_indexes_to[j]]*cos(va[f_bus_indexes[j]]-va[f_bus_indexes_to[j]])) -
                (-f_b[j]*f_tr[j]-f_g[j]*f_ti[j])/f_tm[j]*(vm[f_bus_indexes_to[j]]*sin(va[f_bus_indexes[j]]-va[f_bus_indexes_to[j]]))) for j in f_branch_index) +
                sum(lambda5[j]*(-(-f_g[j]*f_tr[j]-f_b[j]*f_ti[j])/f_tm[j]*(vm[f_bus_indexes_to[j]]*cos(va[f_bus_indexes_to[j]]-va[f_bus_indexes[j]])) -
                (-f_b[j]*f_tr[j]+f_g[j]*f_ti[j])/f_tm[j]*(vm[f_bus_indexes_to[j]]*sin(va[f_bus_indexes_to[j]]-va[f_bus_indexes[j]])) ) for j in f_branch_index) +
                sum(lambda6[j]*((f_b[j]+f_b_fr[j])*vm[f_bus_indexes[j]]*2 +
                (-f_b[j]*f_tr[j]-f_g[j]*f_ti[j])/f_tm[j]*(vm[f_bus_indexes_to[j]]*cos(va[f_bus_indexes[j]]-va[f_bus_indexes_to[j]])) -
                (-f_g[j]*f_tr[j]+f_b[j]*f_ti[j])/f_tm[j]*(vm[f_bus_indexes_to[j]]*sin(va[f_bus_indexes[j]]-va[f_bus_indexes_to[j]]))) for j in f_branch_index) +
                sum(lambda7[j]*((-f_b[j]*f_tr[j]+f_g[j]*f_ti[j])/f_tm[j]*(vm[f_bus_indexes_to[j]]*cos(va[f_bus_indexes[j]]-va[f_bus_indexes_to[j]])) -
                (-f_g[j]*f_tr[j]-f_b[j]*f_ti[j])/f_tm[j]*(vm[f_bus_indexes_to[j]]*sin(va[f_bus_indexes_to[j]]-va[f_bus_indexes[j]])) ) for j in f_branch_index)
                == 0)
        end
        

    for arc in ref[:arcs]
        (l,i,j) = arc
        branch = ref[:branch][l]
        if haskey(branch, "rate_a")
            JuMP.set_lower_bound(p[arc], -branch["rate_a"])
            JuMP.set_upper_bound(p[arc],  branch["rate_a"])
            JuMP.set_lower_bound(q[arc], -branch["rate_a"])
            JuMP.set_upper_bound(q[arc],  branch["rate_a"])
        end
    end

    # include penalties?
    if penalize == false
        @NLobjective(model, Max,
        sum(log(sum(abs(pg[i]-Pgen_final[i,j])        for i in opt_set[:pg_varset][j])) for j in opt_set[:pg_dataset])+  #1:size(Pgen_final,2))+
        #sum(log(sum(abs(qg[i]-Qgen_final[i,j]) for (i,gen) in ref[:gen])) for j in 1:size(Qgen_final,2))+
        sum(log(sum(abs(vm[i]-Voltage_mag_final[i,j]) for i in opt_set[:vm_varset][j])) for j in opt_set[:vm_dataset])+  #1:size(Voltage_mag_final,2))+
        sum(log(sum(abs(pd[i]-Pd_final[i,j])          for i in opt_set[:pd_varset][j])) for j in opt_set[:pd_dataset]))  #1:size(Pd_final,2)))
    else
        @NLobjective(model, Max,
        sum(log(sum(abs(pg[i]-Pgen_final[i,j])        for i in opt_set[:pg_varset][j])) for j in opt_set[:pg_dataset])+  #1:size(Pgen_final,2))+
        #sum(log(sum(abs(qg[i]-Qgen_final[i,j]) for (i,gen) in ref[:gen])) for j in 1:size(Qgen_final,2))+
        sum(log(sum(abs(vm[i]-Voltage_mag_final[i,j]) for i in opt_set[:vm_varset][j])) for j in opt_set[:vm_dataset])+  #1:size(Voltage_mag_final,2))+
        sum(log(sum(abs(pd[i]-Pd_final[i,j])          for i in opt_set[:pd_varset][j])) for j in opt_set[:pd_dataset])+  #1:size(Pd_final,2)))
        penalty_weight*penalty_terms)
    end

    # Fix the voltage angle to zero at the reference bus
    for (i,bus) in ref[:ref_buses]
        # Refrence Bus
        @constraint(model, va[i] == 0)
    end

    # Nodal power balance constraints
    for (i,bus) in ref[:bus]
        # Build a list of the loads and shunt elements connected to the bus i
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        # Active power balance at node i
        @constraint(model,
            sum(p[a] for a in ref[:bus_arcs][i]) ==
            sum(pg[g] for g in ref[:bus_gens][i]) -
            sum(pd[l] for l in ref[:bus_loads][i]) -
            sum(shunt["gs"] for shunt in bus_shunts)*vm[i]^2
        )

        # Reactive power balance at node i
        @constraint(model,
            sum(q[a] for a in ref[:bus_arcs][i])  ==
            sum(qg[g] for g in ref[:bus_gens][i]) -
            sum(qd[l] for l in ref[:bus_loads][i]) +
            sum(shunt["bs"] for shunt in bus_shunts)*vm[i]^2
        )
    end

    # Branch power flow physics and limit constraints
    for (i,branch) in ref[:branch]

        # Build the from and to variable id of the i-th branch, which is a tuple given by (branch id, from bus, to bus)
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])
        # the from and to sides of a branch are useful to calculate power losses later


        p_fr = p[f_idx]  # p_fr is a reference to the optimization variable p[f_idx]
        q_fr = q[f_idx]  # q_fr is a reference to the optimization variable q[f_idx]
        p_to = p[t_idx]  # p_to is a reference to the optimization variable p[t_idx]
        q_to = q[t_idx]  # q_to is a reference to the optimization variable q[t_idx]
        # adding constraints to p_fr is equivalent to adding constraints to p[f_idx], and so on

        vm_fr = vm[branch["f_bus"]]  # vm_fr is a reference to the optimization variable vm on the from side of the branch
        vm_to = vm[branch["t_bus"]]  # vm_to is a reference to the optimization variable vm on the to side of the branch
        va_fr = va[branch["f_bus"]]  # va_fr is a reference to the optimization variable va on the from side of the branch
        va_to = va[branch["t_bus"]]  # va_fr is a reference to the optimization variable va on the to side of the branch

        # Compute the branch parameters and transformer ratios from the data
        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        tm = branch["tap"]^2
        # tap is assumed to be 1.0 on non-transformer branches

        # AC Line Flow Constraints

        # From side of the branch flow
        @NLconstraint(model, p_fr ==  (g+g_fr)/tm*vm_fr^2 + (-g*tr+b*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )
        @NLconstraint(model, q_fr == -(b+b_fr)/tm*vm_fr^2 - (-b*tr-g*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )

        # To side of the branch flow
        @NLconstraint(model, p_to ==  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/tm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )
        @NLconstraint(model, q_to == -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/tm*(vm_to*vm_fr*cos(va_fr-va_to)) + (-g*tr-b*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )

        # Phase Angle Difference Limit
        @constraint(model, va_fr - va_to <= branch["angmax"])
        @constraint(model, va_fr - va_to >= branch["angmin"])

        # Apparent Power Limit, From and To
        if haskey(branch, "rate_a") 
            @constraint(model, p[f_idx]^2 + q[f_idx]^2 <= branch["rate_a"]^2)
            @constraint(model, p[t_idx]^2 + q[t_idx]^2 <= branch["rate_a"]^2)
        end
    end

    # output
    return model, gen_index, bus_index, load_index, branch_index, arcs_index, ref_index
end

function opf_data_collection_powerflow_initialization(data::Dict{String,Any},Pgen_final::Matrix{Float64},Qgen_final::Matrix{Float64},Voltage_mag_final::Matrix{Float64},Voltage_ang_final::Matrix{Float64},Pd_final::Matrix{Float64},Qd_final::Matrix{Float64},epsilon_value::Float64,bounds::Vector{Float64},opt_set,model=Model())

    # soften complementary slackness constraints
    eps_val = -epsilon_value

    # store constraints
    lambda = Dict(:lam1 => Dict(), :lam2 => Dict(), :lam3 => Dict(),
                  :lam4 => Dict(), :lam5 => Dict(), :lam6 => Dict(),
                  :lam7 => Dict())
    mu = Dict(:mu5 => Dict(), :mu6 => Dict(), :mu7 => Dict(), :mu8 => Dict())

    # initialize network
    @assert !haskey(data, "multinetwork")
    @assert !haskey(data, "conductors")

    # Add zeros to turn linear objective functions into quadratic ones
    # so that additional parameter checks are not required
    PowerModels.standardize_cost_terms!(data, order=2)

    # use build_ref to filter out inactive components
    ref = PowerModels.build_ref(data)[:it][pm_it_sym][:nw][nw_id_default]

    # ref contains all the relevant system parameters needed to build the OPF model
    # when we introduce constraints and variable bounds below, we use the parameters in ref.
    gen_index=collect(keys(ref[:gen]))
    bus_index=collect(keys(ref[:bus]))
    load_index=collect(keys(ref[:load]))
    branch_index=collect(keys(ref[:branch]))
    arcs_index=collect(keys(ref[:arcs]))
    ref_index=collect(keys(ref[:ref_buses]))

    # add voltage angles va for each bus
    @variable(model, va[i in keys(ref[:bus])])

    # Add voltage magnitudes vm for each bus
    @variable(model, ref[:bus][i]["vmin"] <= vm[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"])

    # Add active power generation variable pg for each generator (including limits) + random start
    @variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])

    # Add reactive power generation variable qg for each generator (including limits)
    @variable(model, ref[:gen][i]["qmin"] <= qg[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"])

    # Add power flow variables p to represent the active power flow for each branch
    @variable(model, -Inf <= p[(l,i,j) in ref[:arcs]] <= Inf)

    # Add power flow variables q to represent the reactive power flow for each branch
    @variable(model, -Inf <= q[(l,i,j) in ref[:arcs]] <= Inf)

    pd_max = Dict{Int, Float64}()
    pd_min = Dict{Int, Float64}()
    qd_max = Dict{Int, Float64}()
    qd_min = Dict{Int, Float64}()

    superior_bound=bounds[1]
    inferior_bound=bounds[2]

    for ii in 1:length(keys(ref[:load]))
        if Pd_final[ii, 1] > 0.0
            pd_max[ii] = Pd_final[ii, 1]*superior_bound
            pd_min[ii] = Pd_final[ii, 1]*inferior_bound
        elseif Pd_final[ii, 1] < 0.0
            pd_max[ii] = Pd_final[ii, 1]*inferior_bound
            pd_min[ii] = Pd_final[ii, 1]*superior_bound
        else
            pd_max[ii] = 0.0
            pd_min[ii] = 0.0
        end
        if Qd_final[ii,1] > 0.0
            qd_max[ii] = Qd_final[ii, 1]*superior_bound
            qd_min[ii] = Qd_final[ii, 1]*inferior_bound
        elseif Qd_final[ii, 1] < 0.0
            qd_max[ii] = Qd_final[ii, 1]*inferior_bound
            qd_min[ii] = Qd_final[ii, 1]*superior_bound
        else
            qd_max[ii] = 0.0
            qd_min[ii] = 0.0
        end
    end

    # loads
    @variable(model, pd_min[i] <= pd[i in keys(ref[:load])] <= pd_max[i])
    @variable(model, qd_min[i] <= qd[i in keys(ref[:load])] <= qd_max[i])

    # set start values for vm -- uniform distribution between v_min and v_max
    for i in keys(ref[:bus])
        vm_mean  = (ref[:bus][i]["vmin"] + ref[:bus][i]["vmax"])/2
        vm_range = ref[:bus][i]["vmax"] - ref[:bus][i]["vmin"]
        vm_start = vm_mean + vm_range*(rand(1)[1] - 0.5)
        set_start_value.(vm[i],vm_start)
    end

    # set start values for va -- uniform distribution between -0.2 and 0.2
    for i in keys(ref[:bus])
        va_mean  = 0
        va_range = 0.4
        va_start = va_mean + va_range*(rand(1)[1] - 0.5)
        set_start_value.(va[i],va_start)
    end

    # set start values for pg and qg -- uniform distribution between min and max
    for i in keys(ref[:gen])
        # active power
        pg_mean  = (ref[:gen][i]["pmin"] + ref[:gen][i]["pmax"])/2
        pg_range = abs(ref[:gen][i]["pmax"] - ref[:gen][i]["pmin"])
        pg_start = pg_mean + pg_range*(rand(1)[1] - 0.5)
        set_start_value.(pg[i],pg_start)

        # reactive power
        qg_mean  = (ref[:gen][i]["qmin"] + ref[:gen][i]["qmax"])/2
        qg_range = abs(ref[:gen][i]["qmax"] - ref[:gen][i]["qmin"])
        qg_start = qg_mean + qg_range*(rand(1)[1] - 0.5)
        set_start_value.(qg[i],qg_start)
    end

    # set start values for loads
    for i in keys(ref[:load])
        # active power
        pd_mean  = (pd_min[i] + pd_max[i])/2
        pd_range = (pd_max[i] - pd_min[i])
        pd_start = pd_mean + pd_range*(rand(1)[1] - 0.5)
        set_start_value.(pd[i],pd_start)

        # reactive power
        qd_mean  = (qd_min[i] + qd_max[i])/2
        qd_range = (qd_max[i] - qd_min[i])
        qd_start = qd_mean + qd_range*(rand(1)[1] - 0.5)
        set_start_value.(qd[i],qd_start)
    end

    # find the slack bus
    slack_bus = Dict{Int64, Vector{Int64}}()
    for (i, bus) in ref[:bus]
        slack_bus[i] = []
        if i in keys(ref[:ref_buses])
            slack_bus[i] = [i]
        end
    end

    global f_branch_index = []
    global f_bus_indexes_to = Dict()
    global f_bus_indexes = Dict()
    global f_g =    Dict()
    global f_b =    Dict()
    global f_tr =   Dict()
    global f_ti =   Dict()
    global f_g_fr = Dict()
    global f_b_fr = Dict()
    global f_g_to = Dict()
    global f_b_to = Dict()
    global f_tm =   Dict()

    global t_branch_index = []
    global t_bus_indexes_fr = Dict()
    global t_bus_indexes = Dict()
    global t_g =    Dict()
    global t_b =    Dict()
    global t_tr =   Dict()
    global t_ti =   Dict()
    global t_g_fr = Dict()
    global t_b_fr = Dict()
    global t_g_to = Dict()
    global t_b_to = Dict()
    global t_tm =   Dict()

    branchesfinal=ref[:arcs][1:length(ref[:branch])]

        for (i,bus) in ref[:bus]
            bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]
            f_branch_index = []
            f_bus_indexes_to = Dict()
            f_bus_indexes = Dict()
            f_g =    Dict()
            f_b =    Dict()
            f_tr =   Dict()
            f_ti =   Dict()
            f_g_fr = Dict()
            f_b_fr = Dict()
            f_g_to = Dict()
            f_b_to = Dict()
            f_tm =   Dict()

            t_branch_index = []
            t_bus_indexes_fr = Dict()
            t_bus_indexes = Dict()
            t_g =    Dict()
            t_b =    Dict()
            t_tr =   Dict()
            t_ti =   Dict()
            t_g_fr = Dict()
            t_b_fr = Dict()
            t_g_to = Dict()
            t_b_to = Dict()
            t_tm =   Dict()

            for ii in branchesfinal
                if i ==ii[2]
                    append!(f_branch_index, ii[1])
                    f_bus_indexes[ii[1]] = ii[2]
                    f_bus_indexes_to[ii[1]] = ii[3]
                    g, b = PowerModels.calc_branch_y(ref[:branch][ii[1]])
                    tr, ti = PowerModels.calc_branch_t(ref[:branch][ii[1]])
                    f_g[ii[1]]    = g
                    f_b[ii[1]]    = b
                    f_tr[ii[1]]   = tr
                    f_ti[ii[1]]   = ti
                    f_g_fr[ii[1]] = ref[:branch][ii[1]]["g_fr"]
                    f_b_fr[ii[1]] = ref[:branch][ii[1]]["b_fr"]
                    f_g_to[ii[1]] = ref[:branch][ii[1]]["g_to"]
                    f_b_to[ii[1]] = ref[:branch][ii[1]]["b_to"]
                    f_tm[ii[1]] =   ref[:branch][ii[1]]["tap"]^2
                end
                if i ==ii[3]
                    append!(t_branch_index, ii[1])
                    t_bus_indexes[ii[1]] = ii[3]
                    t_bus_indexes_fr[ii[1]] = ii[2]
                    g, b = PowerModels.calc_branch_y(ref[:branch][ii[1]])
                    tr, ti = PowerModels.calc_branch_t(ref[:branch][ii[1]])
                    t_g[ii[1]]    = g
                    t_b[ii[1]]    = b
                    t_tr[ii[1]]   = tr
                    t_ti[ii[1]]   = ti
                    t_g_fr[ii[1]] = ref[:branch][ii[1]]["g_fr"]
                    t_b_fr[ii[1]] = ref[:branch][ii[1]]["b_fr"]
                    t_g_to[ii[1]] = ref[:branch][ii[1]]["g_to"]
                    t_b_to[ii[1]] = ref[:branch][ii[1]]["b_to"]
                    t_tm[ii[1]] =   ref[:branch][ii[1]]["tap"]^2
                end
            end
        end
            
    for arc in ref[:arcs]
        (l,i,j) = arc
        branch = ref[:branch][l]
        if haskey(branch, "rate_a")
            JuMP.set_lower_bound(p[arc], -branch["rate_a"])
            JuMP.set_upper_bound(p[arc],  branch["rate_a"])
            JuMP.set_lower_bound(q[arc], -branch["rate_a"])
            JuMP.set_upper_bound(q[arc],  branch["rate_a"])
        end
    end

    @NLobjective(model, Max,
        sum(log(sum(abs(pg[i]-Pgen_final[i,j])        for i in opt_set[:pg_varset][j])) for j in opt_set[:pg_dataset])+  #1:size(Pgen_final,2))+
        #sum(log(sum(abs(qg[i]-Qgen_final[i,j]) for (i,gen) in ref[:gen])) for j in 1:size(Qgen_final,2))+
        sum(log(sum(abs(vm[i]-Voltage_mag_final[i,j]) for i in opt_set[:vm_varset][j])) for j in opt_set[:vm_dataset])+  #1:size(Voltage_mag_final,2))+
        sum(log(sum(abs(pd[i]-Pd_final[i,j])          for i in opt_set[:pd_varset][j])) for j in opt_set[:pd_dataset]))  #1:size(Pd_final,2)))

    # Fix the voltage angle to zero at the reference bus
    for (i,bus) in ref[:ref_buses]
        # Refrence Bus
        lambda[:lam1][i] = @constraint(model, va[i] == 0)
    end

    # Nodal power balance constraints
    for (i,bus) in ref[:bus]
        # Build a list of the loads and shunt elements connected to the bus i
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        # Active power balance at node i
        lambda[:lam2][i] = @constraint(model,
            sum(p[a] for a in ref[:bus_arcs][i]) ==
            sum(pg[g] for g in ref[:bus_gens][i]) -
            sum(pd[l] for l in ref[:bus_loads][i]) -
            sum(shunt["gs"] for shunt in bus_shunts)*vm[i]^2
        )

        # Reactive power balance at node i
        lambda[:lam3][i] = @constraint(model,
            sum(q[a] for a in ref[:bus_arcs][i])  ==
            sum(qg[g] for g in ref[:bus_gens][i]) -
            sum(qd[l] for l in ref[:bus_loads][i]) +
            sum(shunt["bs"] for shunt in bus_shunts)*vm[i]^2
        )
    end

    # Branch power flow physics and limit constraints
    for (i,branch) in ref[:branch]

        # Build the from and to variable id of the i-th branch, which is a tuple given by (branch id, from bus, to bus)
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])
        # the from and to sides of a branch are useful to calculate power losses later

        p_fr = p[f_idx]  # p_fr is a reference to the optimization variable p[f_idx]
        q_fr = q[f_idx]  # q_fr is a reference to the optimization variable q[f_idx]
        p_to = p[t_idx]  # p_to is a reference to the optimization variable p[t_idx]
        q_to = q[t_idx]  # q_to is a reference to the optimization variable q[t_idx]
        # adding constraints to p_fr is equivalent to adding constraints to p[f_idx], and so on

        vm_fr = vm[branch["f_bus"]]  # vm_fr is a reference to the optimization variable vm on the from side of the branch
        vm_to = vm[branch["t_bus"]]  # vm_to is a reference to the optimization variable vm on the to side of the branch
        va_fr = va[branch["f_bus"]]  # va_fr is a reference to the optimization variable va on the from side of the branch
        va_to = va[branch["t_bus"]]  # va_fr is a reference to the optimization variable va on the to side of the branch

        # Compute the branch parameters and transformer ratios from the data
        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        tm = branch["tap"]^2
        # tap is assumed to be 1.0 on non-transformer branches

        # From side of the branch flow
        lambda[:lam4][i] = @NLconstraint(model, p_fr ==  (g+g_fr)/tm*vm_fr^2 + (-g*tr+b*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )
        lambda[:lam5][i] = @NLconstraint(model, q_fr == -(b+b_fr)/tm*vm_fr^2 - (-b*tr-g*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )

        # To side of the branch flow
        lambda[:lam6][i] = @NLconstraint(model, p_to ==  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/tm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )
        lambda[:lam7][i] = @NLconstraint(model, q_to == -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/tm*(vm_to*vm_fr*cos(va_fr-va_to)) + (-g*tr-b*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )

        # Phase Angle Difference Limit
        mu[:mu5][i] = @constraint(model, va_fr - va_to <= branch["angmax"])
        mu[:mu6][i] = @constraint(model, va_fr - va_to >= branch["angmin"])

        # Apparent Power Limit, From and To
        if haskey(branch, "rate_a")
            mu[:mu7][i] = @constraint(model, p[f_idx]^2 + q[f_idx]^2 <= branch["rate_a"]^2)
            mu[:mu8][i] = @constraint(model, p[t_idx]^2 + q[t_idx]^2 <= branch["rate_a"]^2)
        end
    end

    # output
    return model, gen_index, bus_index, load_index, branch_index, arcs_index, ref_index, ref, lambda, mu
end
