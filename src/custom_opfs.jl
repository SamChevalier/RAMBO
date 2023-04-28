include("./opti_relaxed_constraints.jl")
include("./add_opf_final_data.jl")


function update_bilevel_result(case, Pgen_final, Qgen_final, Voltage_mag_final, Voltage_ang_final, Pd_final, Qd_final, Pbranch, Qbranch, N, epsilon_value, bounds)
    check_simple_OPF_vm=Dict()
    check_simple_OPF_pg = Dict()
    Pg_opf = copy(Pgen_final)
    Qg_opf = copy(Qgen_final)
    Vm_opf = copy(Voltage_mag_final)
    Va_opf = copy(Voltage_ang_final)
    Pd_opf = copy(Pd_final)
    Qd_opf = copy(Qd_final)
    Pbranch_opf = copy(Pbranch)
    Qbranch_opf = copy(Qbranch)
    nb = length(Voltage_mag_final)
    nl = length(Pbranch_opf)
    ng = length(Pg_opf)
    nd = length(Pd_opf)

    # call ref to get the relevant keys
    outer_ref  = PowerModels.build_ref(case)[:it][pm_it_sym][:nw][nw_id_default]
    gen_index  = collect(keys(outer_ref[:gen]))
    nng        = length(gen_index)
    bus_index  = collect(keys(outer_ref[:bus]))
    nnb        = length(bus_index)
    load_index = collect(keys(outer_ref[:load]))
    nnl        = length(load_index)   

    # maximum number of previous solutions to optimize over
    max_solns = 30
    max_vars  = 30
    opt_set   = Dict(:pg_dataset => Int64[],
                     :pd_dataset => Int64[],
                     :vm_dataset => Int64[],
                     :pg_varset =>  Vector{Vector{Int64}}[],
                     :pd_varset =>  Vector{Vector{Int64}}[],
                     :vm_varset =>  Vector{Vector{Int64}}[])

    # initialization dictionary
    init_soln = Dict(:vm  => Dict(), :va  => Dict(), :pg  => Dict(),
        :qg  => Dict(), :pd  => Dict(), :qd  => Dict(),
        :p   => Dict(), :q   => Dict(), 
        :lam1 => Dict(), :lam2 => Dict(), :lam3 => Dict(),
        :lam4 => Dict(), :lam5 => Dict(), :lam6 => Dict(), :lam7 => Dict(),
        :mu1 => Dict(), :mu2 => Dict(), :mu3 => Dict(), :mu4 => Dict(), :mu5 => Dict(),
        :mu6 => Dict(), :mu7 => Dict(), :mu8 => Dict(), :mu9 => Dict(), :mu10 => Dict())

    while feasible_count_opf< N
        # first, get a powerflow solution to initialize the bilevel solver -- we
        # run this over and over until we find a solution
        #
        # loop until we find a solution
        soln_found = 0
        while soln_found == 0
            # only optimize over a (random subset) of previous solutions
            nd                   = size(Pgen_final,2)
            opt_set[:pg_dataset] = copy(shuffle(1:nd)[1:min(max_solns,nd)])
            opt_set[:pd_dataset] = copy(shuffle(1:nd)[1:min(max_solns,nd)])
            opt_set[:vm_dataset] = copy(shuffle(1:nd)[1:min(max_solns,nd)])

            # for each...
            nv                  = min(max_vars,nng)
            opt_set[:pg_varset] = [shuffle(gen_index)[1:min(max_vars,nng)]  for ii in 1:nd]
            opt_set[:pd_varset] = [shuffle(load_index)[1:min(max_vars,nnl)] for ii in 1:nd]
            opt_set[:vm_varset] = [shuffle(bus_index)[1:min(max_vars,nnb)]  for ii in 1:nd]

            # optimizer
            T_max     = 150.0
            optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "max_cpu_time" => T_max)
            model     = Model(optimizer)

            # set up
            opt_output = opf_data_collection_powerflow_initialization(case, Pgen_final[:, :], Qgen_final[:, :], Voltage_mag_final[:, :], Voltage_ang_final[:, :], Pd_final[:, :], Qd_final[:, :], epsilon_value, bounds, opt_set, model)
            println("=== intialization opt start ===")
            result = optimize!(opt_output[1])
            println("=== intialization opt end ===")

            if string(termination_status(opt_output[1])) in ["LOCALLY_SOLVED", "ALMOST_LOCALLY_SOLVED"]
                soln_found = 1
                # parse the solution -- updates init_soln dictionary
                parse_initialization!(init_soln, opt_output[1], opt_output[8], opt_output[9], opt_output[10])
            end
        end

        # now, solve bilevel
        T_max      = 150.0
        optimizer  = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "max_cpu_time" => T_max)
        model      = Model(optimizer)
        opt_output = opf_data_collection_bilevel(init_soln, case, Pgen_final[:, :], Qgen_final[:, :], Voltage_mag_final[:, :], Voltage_ang_final[:, :], Pd_final[:, :], Qd_final[:, :], epsilon_value, bounds, opt_set, model)
        println("=== bilevel opt start ===")
        result = optimize!(opt_output[1])
        t = solve_time(opt_output[1])
        println("=== bilevel solve time: $t ===")
        println("=== bilevel opt end ===")

        println("Objective: =============================================")
        println(objective_value(opt_output[1]))
        println("=============================================")

        # inds
        gen_index = opt_output[2]
        bus_index = opt_output[3]
        load_index = opt_output[4]


        if string(termination_status(opt_output[1])) in ["LOCALLY_SOLVED", "ALMOST_LOCALLY_SOLVED"] # , "ALMOST_LOCALLY_SOLVED"
            global feasible_count += 1
        else
            println("***************** NOT FEASIBLE ----> ", feasible_count_opf)
        end

        Objective_value = objective_value(opt_output[1])
        
        Pgen = value.(opt_output[1][:pg])
        Pgen_final_before_sortt = Pgen.data
        C = hcat(Pgen_final_before_sortt, gen_index)
        D = C[sortperm(C[:, 2]), :]
        Pgenn = D[:, 1]

        Qgen = value.(opt_output[1][:qg])
        Qgen_final_before_sortt = Qgen.data
        C1 = hcat(Qgen_final_before_sortt, gen_index)
        D1 = C1[sortperm(C1[:, 2]), :]
        Qgenn = D1[:, 1]

        Pgen_final = hcat(Pgen_final, Pgenn)
        Qgen_final = hcat(Qgen_final, Qgenn)

        Voltage_mag_before_sort = value.(opt_output[1][:vm])
        E = hcat(Voltage_mag_before_sort, bus_index)
        F = E[sortperm(E[:, 2]), :]
        Voltage_mag = F[:, 1]

        Voltage_ang_before_sort = value.(opt_output[1][:va])
        E1 = hcat(Voltage_ang_before_sort, bus_index)
        F1 = E1[sortperm(E1[:, 2]), :]
        Voltage_ang = F1[:, 1]

        Voltage_mag_final = hcat(Voltage_mag_final, Voltage_mag)
        Voltage_ang_final = hcat(Voltage_ang_final, Voltage_ang)

        println("Iteration=$feasible_count_opf")
        # Check that the solver terminated without an error
        println("The solver termination status is $(termination_status(opt_output[1]))")

        Pgen_final = Matrix{Float64}(Pgen_final)
        Qgen_final = Matrix{Float64}(Qgen_final)
        Voltage_mag_final = Matrix{Float64}(Voltage_mag_final)
        Voltage_ang_final = Matrix{Float64}(Voltage_ang_final)

        Pd_before_sort = value.(opt_output[1][:pd])
        EE1 = hcat(Pd_before_sort, load_index)
        FF1 = EE1[sortperm(EE1[:, 2]), :]
        Pdseguim = FF1[:, 1]

        Qd_before_sort = value.(opt_output[1][:qd])
        EE2 = hcat(Qd_before_sort, load_index)
        FF2 = EE2[sortperm(EE2[:, 2]), :]
        Qdseguim = FF2[:, 1]

        Pd_final = hcat(Pd_final, Pdseguim)
        Qd_final = hcat(Qd_final, Qdseguim)

        Pd_final = Matrix{Float64}(Pd_final)
        Qd_final = Matrix{Float64}(Qd_final)
        
        for ii=1:length(case["load"])
            case["load"][string(ii)]["pd"]=Pdseguim[ii]
            case["load"][string(ii)]["qd"]=Qdseguim[ii]
        end
        
        if string(termination_status(opt_output[1])) in ["LOCALLY_SOLVED", "ALMOST_LOCALLY_SOLVED"] # 
            T_max     = 150.0
            optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "max_cpu_time" => T_max)
            OPF_soln=PowerModels.solve_opf(case, ACPPowerModel, optimizer)
            if OPF_soln["termination_status"] in [LOCALLY_SOLVED, ALMOST_LOCALLY_SOLVED]
                global feasible_count_opf += 1
            end

            Pg_opf_v, Qg_opf_v, Vm_opf_v, Va_opf_v, Pd_opf_v, Qd_opf_v, Pbranch_opf_v, Qbranch_opf_v = create_final_opf_matrix(case, OPF_soln)

            Pg_opf = hcat(Pg_opf, Pg_opf_v)
            Qg_opf = hcat(Qg_opf, Qg_opf_v)
            Vm_opf = hcat(Vm_opf, Vm_opf_v)
            Va_opf = hcat(Va_opf, Va_opf_v)
            Pd_opf = hcat(Pd_opf, Pd_opf_v)
            Qd_opf = hcat(Qd_opf, Qd_opf_v)
            Pbranch_opf = hcat(Pbranch_opf, Pbranch_opf_v)
            Qbranch_opf = hcat(Qbranch_opf, Qbranch_opf_v)



            load_system=0
            for i =1:length(case["load"])
                load_system =load_system+case["load"][string(i)]["pd"]
            end

            generation_system =0
            for i=1:length(OPF_soln["solution"]["gen"])
                generation_system= generation_system+OPF_soln["solution"]["gen"][string(i)]["pg"]
            end
            
            
            println("load system ",load_system)
            println("load system bi ", sum(Pdseguim), " ", length(Pdseguim))
            println("generation system ", sum(generation_system))
            println("generation system bi ", sum(Pgenn))
        
            coincidences_vm = []
            for iii=1:length(case["bus"])
                value_coincidence_vm = 100*abs(OPF_soln["solution"]["bus"][string(iii)]["vm"]-Voltage_mag[iii])/OPF_soln["solution"]["bus"][string(iii)]["vm"]
                append!(coincidences_vm, value_coincidence_vm)
            end
            average_coin_percent= sum(coincidences_vm)/length(coincidences_vm)
            check_simple_OPF_vm[feasible_count_opf]= [OPF_soln["termination_status"], average_coin_percent]
            println(average_coin_percent)

            coincidences_pg =Dict()
            for iii=1:length(case["gen"])
                if case["gen"][string(iii)]["pmax"] != 0.0
                    value_coincidence_pg = 100*abs(OPF_soln["solution"]["gen"][string(iii)]["pg"]-Pgenn[iii])/OPF_soln["solution"]["gen"][string(iii)]["pg"]
                    if abs(OPF_soln["solution"]["gen"][string(iii)]["pg"]-Pgenn[iii]) < 5e-3
                        #comparing two very small values can lead to huge differences e.g. 2e-7/4e-9
                        coincidences_pg[iii] =0.0
                    else
                        coincidences_pg[iii]= value_coincidence_pg
                        println("coinc ",iii)
                        println("opf ", OPF_soln["solution"]["gen"][string(iii)]["pg"])
                        println("bi ",Pgenn[iii])
                        println(value_coincidence_pg)
                    end
                end
            end
            check_simple_OPF_pg[feasible_count_opf]= [OPF_soln["termination_status"], coincidences_pg]
        else
            global infeasible_count_opf += 1
        end
        println("Iteration ", feasible_count_opf, " of ", N)
        
    end
    # output dict
    opf_data = Dict(:P_gen => Pgen_final,
                    :Q_gen => Qgen_final,
                    :Vm    => Voltage_mag_final,
                    :Va    => Voltage_ang_final,
                    :Pd    => Pd_final,
                    :Qd    => Qd_final,
                    :Verification_vm => check_simple_OPF_vm,
                    :Verification_pg => check_simple_OPF_pg)

    opf_dataset_final = Dict(:Pg_opf => Pg_opf,
                    :Qg_opf => Qg_opf,
                    :Vm_opf    => Vm_opf,
                    :Va_opf    => Va_opf,
                    :Pd_opf    => Pd_opf,
                    :Qd_opf    => Qd_opf,
                    :Pbranch_opf => Pbranch_opf,
                    :Qbranch_opf => Qbranch_opf)

    if N==feasible_count
        println("**** Great ! All the bilvevel solutions where feasible")
        println("**** There were ", infeasible_count_opf, " OPF problems that were infeasible")
    else
        println("**** Check bilvevel solutions! Not all of them were feasible")
        println(feasible_count, " / ", N, " were feasible")
    end
    # output
    return opf_data, opf_dataset_final, init_soln

end

function initialize_bilevel_solver(case)

    # Solving the optimal power flow using the standard function of Power Models
    T_max     = 150.0
    optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "max_cpu_time" => T_max)
    OPF_soln  = PowerModels.solve_opf(case, ACPPowerModel, optimizer)
    
    firststatus = println(OPF_soln["termination_status"]) #answers the solution of the model


    llista_gen = []
    llista_bus = []
    llista_demand = []
    llista_branch = []

    for i in collect(keys(case["bus"]))
        append!(llista_bus, convert(AbstractFloat, case["bus"][i]["bus_i"]))
    end

    for i in collect(keys(case["gen"]))
        if case["gen"][i]["gen_status"] == 1
            append!(llista_gen, convert(AbstractFloat, case["gen"][i]["gen_bus"]))
        end
        if case["gen"][i]["gen_status"] != 1
            throw(DomainError("case data error", "There are generators with a state variabledifferent than 0"))
        end
    end

    for i in collect(keys(case["load"]))
        append!(llista_demand, convert(AbstractFloat, case["load"][i]["load_bus"]))
    end

    for i in collect(keys(case["branch"]))
        append!(llista_branch, convert(AbstractFloat, case["branch"][i]["index"]))
    end 

    Pginicial = []
    for ii in collect(keys(case["gen"]))
        if case["gen"][ii]["gen_status"] == 1
            append!(Pginicial, convert(AbstractFloat, OPF_soln["solution"]["gen"][ii]["pg"]))
        end
    end

    Qginicial = []
    for ii in collect(keys(case["gen"]))
        if case["gen"][ii]["gen_status"] == 1
            append!(Qginicial, convert(AbstractFloat, OPF_soln["solution"]["gen"][ii]["qg"]))
        end
    end

    Vminicial = []
    for ii in collect(keys(case["bus"]))
        append!(Vminicial, convert(AbstractFloat, OPF_soln["solution"]["bus"][ii]["vm"]))
    end

    Vainicial = []
    for ii in collect(keys(case["bus"]))
        append!(Vainicial, convert(AbstractFloat, OPF_soln["solution"]["bus"][ii]["va"]))
    end

    Pdinicial = []
    for ii in collect(keys(case["load"]))
        if case["load"][ii]["status"] == 1
            append!(Pdinicial, convert(AbstractFloat, case["load"][ii]["pd"]))
        end
    end

    Qdinicial = []
    for ii in collect(keys(case["load"]))
        if case["load"][ii]["status"] == 1
            append!(Qdinicial, convert(AbstractFloat, case["load"][ii]["qd"]))
        end
    end

    Pbranch_i = []
    Qbranch_i = []
    for ii in collect(keys(case["branch"]))
        qf = convert(AbstractFloat, OPF_soln["solution"]["branch"][ii]["qf"])
        qt = convert(AbstractFloat, OPF_soln["solution"]["branch"][ii]["qt"])
        pf = convert(AbstractFloat, OPF_soln["solution"]["branch"][ii]["pf"])
        pt = convert(AbstractFloat, OPF_soln["solution"]["branch"][ii]["pt"])
        s_from = pf^2 + qf^2
        s_to = pt^2 + qt^2
        if s_from > s_to
            append!(Pbranch_i, pf)
            append!(Qbranch_i, qf)
        else
            append!(Pbranch_i, pt)
            append!(Qbranch_i, qt)
        end
    end

    Vminter1 = hcat(Vminicial, llista_bus)
    Vminter2 = Vminter1[sortperm(Vminter1[:, 2]), :]
    Vm = Vminter2[:, 1]

    Vainter1 = hcat(Vainicial, llista_bus)
    Vainter2 = Vainter1[sortperm(Vainter1[:, 2]), :]
    Va = Vainter2[:, 1]

    Pginter1 = hcat(Pginicial, llista_gen)
    Pginter2 = Pginter1[sortperm(Pginter1[:, 2]), :]
    P_gen_def = Pginter2[:, 1]

    Qginter1 = hcat(Qginicial, llista_gen)
    Qginter2 = Qginter1[sortperm(Qginter1[:, 2]), :]
    Q_gen_def = Qginter2[:, 1]

    Pdinter = hcat(Pdinicial, llista_demand)
    Pdinter1 = Pdinter[sortperm(Pdinter[:, 2]), :]
    Pd_def = Pdinter1[:, 1]

    Qdinter = hcat(Qdinicial, llista_demand)
    Qdinter1 = Qdinter[sortperm(Qdinter[:, 2]), :]
    Qd_def = Qdinter1[:, 1]

    Pbranch_inter = hcat(Pbranch_i, llista_branch)
    Pbranch_inter1 = Pbranch_inter[sortperm(Pbranch_inter[:, 2]), :]
    Pbranch__def = Pbranch_inter1[:, 1]

    Qbranch_inter = hcat(Qbranch_i, llista_branch)
    Qbranch_inter1 = Qbranch_inter[sortperm(Qbranch_inter[:, 2]), :]
    Qbranch__def = Qbranch_inter1[:, 1]

    return P_gen_def, Q_gen_def, Vm, Va, Pd_def, Qd_def, Pbranch__def, Qbranch__def

end

# evaluate constraints
function eval_constraints(ref, pg, qg, vm, va, p, q)
    c1_set  = Int64[]
    c2_set  = Int64[]
    c3_set  = Int64[]
    c4_set  = Int64[]
    c5_set  = Int64[]
    c6_set  = Int64[]
    c7_set  = Int64[]
    c8_set  = Int64[]
    c9_set  = Int64[]
    c10_set = Int64[]
    total   = 0
    total_added = 0


    for (i,gen) in ref[:gen]
        if abs(value((pg[i]-ref[:gen][i]["pmax"]))) < 0.25
            push!(c1_set,i)
            total_added += 1
        end
        if abs(value((ref[:gen][i]["pmin"]-pg[i]))) < 0.25
            push!(c2_set,i)
            total_added += 1
        end
        if abs(value((qg[i]-ref[:gen][i]["qmax"]))) < 0.25
            push!(c3_set,i)
            total_added += 1
        end
        if abs(value((ref[:gen][i]["qmin"]-qg[i]))) < 0.25
            push!(c4_set,i)
            total_added += 1
        end
        total += 4
    end

    for (i,branch) in ref[:branch]
        va_fr = va[branch["f_bus"]]  
        va_to = va[branch["t_bus"]]  

        # Phase Angle Difference Limit
        if abs(value((branch["angmin"] - va_fr + va_to))  ) < 0.25
            push!(c5_set,i)
            total_added += 1
        end
        if abs(value((va_fr - va_to - branch["angmax"]))  ) < 0.25
            push!(c6_set,i)
            total_added += 1
        end
        total += 2
    end

    for (i,branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])  # 20 variables correctly sorted
        t_idx = (i, branch["t_bus"], branch["f_bus"])  # 20 variables correctly sorted
        

        if haskey(branch, "rate_a")  
            if abs(value((p[f_idx]^2 + q[f_idx]^2 - branch["rate_a"]^2))) < 0.25
                push!(c7_set,i)
                total_added += 1
            end
            if abs(value((p[t_idx]^2 + q[t_idx]^2 - branch["rate_a"]^2))) < 0.25
                push!(c8_set,i)
                total_added += 1
            end
        end
        total += 2
    end

    for (i,bus) in ref[:bus]
        if abs(value((vm[i]-ref[:bus][i]["vmax"]))) < 0.25
            push!(c9_set,i)
            total_added += 1
        end
        if abs(value((ref[:bus][i]["vmin"]-vm[i]))) < 0.25
            push!(c10_set,i)
            total_added += 1
        end
        total += 2
    end

    c_set = Dict(:c1 => c1_set,
    :c2 => c2_set,
    :c3 => c3_set,
    :c4 => c4_set,
    :c5 => c5_set,
    :c6 => c6_set,
    :c7 => c7_set,
    :c8 => c8_set,
    :c9 => c9_set,
    :c10 => c10_set)

    println("==============")
    println("total added: $total_added")
    println("total: $total")
    println("==============")

    return c_set
end

# parse
function parse_initialization!(init_soln, opf_result, ref, lambda, mu)
    for i in keys(ref[:bus])
        init_soln[:vm][i] = value.(opf_result[:vm][i])
        init_soln[:va][i] = value.(opf_result[:va][i])

        init_soln[:mu9][i]  = abs(dual(UpperBoundRef(opf_result[:vm][i])))
        init_soln[:mu10][i] = abs(dual(LowerBoundRef(opf_result[:vm][i])))
    end

    for i in keys(ref[:gen])
        init_soln[:pg][i] = value.(opf_result[:pg][i])
        init_soln[:qg][i] = value.(opf_result[:qg][i])

        init_soln[:mu1][i] = abs(dual(UpperBoundRef(opf_result[:pg][i])))
        init_soln[:mu2][i] = abs(dual(LowerBoundRef(opf_result[:pg][i])))
        init_soln[:mu3][i] = abs(dual(UpperBoundRef(opf_result[:qg][i])))
        init_soln[:mu4][i] = abs(dual(LowerBoundRef(opf_result[:qg][i])))
    end
    for i in keys(ref[:load])
        init_soln[:pd][i] = value.(opf_result[:pd][i])
        init_soln[:qd][i] = value.(opf_result[:qd][i])
    end
    for (l,i,j) in ref[:arcs]
        init_soln[:p][(l,i,j)]  = value.(opf_result[:p][(l,i,j)])
        init_soln[:q][(l,i,j)]  = value.(opf_result[:p][(l,i,j)])
    end

    # let's also get the duals associated with the equality constraints:
    # we don't know the correct signs, so we need to solve for them!

    for (i,bus) in ref[:ref_buses]
        init_soln[:lam1][i] = dual(lambda[:lam1][i])
    end

    for (i,bus) in ref[:bus]
        init_soln[:lam2][i] = dual(lambda[:lam2][i])
        init_soln[:lam3][i] = dual(lambda[:lam3][i])
    end

    for (i,branch) in ref[:branch]
        init_soln[:lam4][i] = dual(lambda[:lam4][i])
        init_soln[:lam5][i] = dual(lambda[:lam5][i])
        init_soln[:lam6][i] = dual(lambda[:lam6][i])
        init_soln[:lam7][i] = dual(lambda[:lam7][i])

        init_soln[:mu5][i] = abs(dual(mu[:mu5][i]))
        init_soln[:mu6][i] = abs(dual(mu[:mu6][i]))
        init_soln[:mu7][i] = abs(dual(mu[:mu7][i]))
        init_soln[:mu8][i] = abs(dual(mu[:mu8][i]))
    end
end