function create_final_opf_matrix(case, OPF_soln)

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
        sf = sqrt((pf^2) + (qf^2))
        st = sqrt((pt^2) + (qt^2))
        if sf > st
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