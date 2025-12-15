function variable_mc_load_power_VoltVar(pm::AbstractExplicitNeutralIVRModelVoltVar; nw=nw_id_default, report::Bool=true)
    #Get number of phases connected to load and load_ids with VoltVAr control
    int_dim = Dict(i => PMD._infer_int_dim_unit(load, false) for (i,load) in PMD.ref(pm, nw, :load))
    load_ids_VoltVAr = [id for (id,load) in PMD.ref(pm, nw, :load) if load["PV_setpoint"]=="VoltVAr"]
    #Define variables with starting values
    pd = PMD.var(pm, nw)[:pd] = Dict{Int,Any}(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_pd_$(i)",
            start = PMD.comp_start_value(PMD.ref(pm, nw, :load, i), "pd_start", c, 0.0)
        ) for i in load_ids_VoltVAr
    )
    qd = var(pm, nw)[:qd] = Dict{Int,Any}(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_qd_$(i)",
            start = PMD.comp_start_value(ref(pm, nw, :load, i), "qd_start", c, 0.0)
        ) for i in load_ids_VoltVAr
    )

    lambda = PMD.var(pm, nw)[:lambda] = Dict{Int,Any}(i => Dict{Int,Any}(
        c => JuMP.@variable(pm.model,
            [k in 1:4],
            base_name="$(nw)_lambda_$(i)_c$(c)",
            start = PMD.comp_start_value(PMD.ref(pm, nw, :load, i), "lambda_start", k, 0.0))
        for c in 1:int_dim[i]
    ) for i in load_ids_VoltVAr
    )
    #Ensure variables are reported in the solution
    report && IM.sol_component_value(pm, pmd_it_sym, nw, :load, :pd, load_ids_VoltVAr, pd)
    report && IM.sol_component_value(pm, pmd_it_sym, nw, :load, :qd, load_ids_VoltVAr, qd)
    report && IM.sol_component_value(pm, pmd_it_sym, nw, :load, :lambda, load_ids_VoltVAr, lambda)
end

function constraint_mc_load_current(pm::AbstractExplicitNeutralIVRModelVoltVar, id::Int; nw::Int=nw_id_default, report::Bool=true)
    load = PMD.ref(pm, nw, :load, id)
    bus = PMD.ref(pm, nw, :bus, load["load_bus"])

    PV_setpoint = load["PV_setpoint"]

    if PV_setpoint != "VoltVAr"
        a, alpha,b, beta = PMD._load_expmodel_params(load, bus)
        PMD.constraint_mc_load_current_wye(pm, nw, id, load["load_bus"], load["connections"], a, alpha, b, beta; report=report)
    else
        constraint_mc_load_current_wye_VoltVar(pm, id, load["load_bus"], load["connections"]; nw=nw)
    end
end

function constraint_mc_load_current_wye_VoltVar(pm::AbstractExplicitNeutralIVRModelVoltVar, id::Int, bus::Int, connections::Vector{Int}; nw::Int=nw_id_default, report::Bool=true) # id =  load_id
    vr = PMD.var(pm, nw, :vr, bus)
    vi = PMD.var(pm, nw, :vi, bus)
    pd = PMD.var(pm, nw, :pd, id)
    qd = PMD.var(pm, nw, :qd, id)

    crd = JuMP.NonlinearExpr[]
    cid = JuMP.NonlinearExpr[]

    phases = connections[1:end-1]
    n = connections[end] #neutral conductor

    for (idx, p) in enumerate(phases)
        push!(crd, JuMP.@expression(pm.model, ((vr[p]-vr[n])*pd[idx] + (vi[p]-vi[n])*qd[idx])/((vr[p]-vr[n])^2 + (vi[p]-vi[n])^2)))
        push!(cid, JuMP.@expression(pm.model, ((vi[p]-vi[n])*pd[idx] - (vr[p]-vr[n])*qd[idx])/((vr[p]-vr[n])^2 + (vi[p]-vi[n])^2)))
    end

    PMD.var(pm, nw, :crd)[id] = crd
    PMD.var(pm, nw, :cid)[id] = cid

    crd_bus_n = JuMP.@expression(pm.model, -sum(crd[i] for i in 1:length(phases)))
    cid_bus_n = JuMP.@expression(pm.model, -sum(cid[i] for i in 1:length(phases)))

    PMD.var(pm, nw, :crd_bus)[id] = crd_bus = PMD._merge_bus_flows(pm, [crd..., crd_bus_n], connections)
    PMD.var(pm, nw, :cid_bus)[id] = cid_bus = PMD._merge_bus_flows(pm, [cid..., cid_bus_n], connections)

    if report
        pd_bus = JuMP.NonlinearExpr[]
        qd_bus = JuMP.NonlinearExpr[]
        for (idx, c) in enumerate(connections)
            push!(pd_bus, JuMP.@expression(pm.model, crd_bus[idx]*vr[c] + cid_bus[idx]*vi[c]))
            push!(qd_bus, JuMP.@expression(pm.model, -crd_bus[idx]*vi[c] + cid_bus[idx]*vr[c]))
        end
        PMD.sol(pm, nw, :load, id)[:pd_bus] = JuMP.Containers.DenseAxisArray(pd_bus, connections)
        PMD.sol(pm, nw, :load, id)[:qd_bus] = JuMP.Containers.DenseAxisArray(qd_bus, connections)

        PMD.sol(pm, nw, :load, id)[:crd] = JuMP.Containers.DenseAxisArray(crd, connections) #Why connections and not connections[1:end-1]?
        PMD.sol(pm, nw, :load, id)[:cid] = JuMP.Containers.DenseAxisArray(cid, connections)

        PMD.sol(pm, nw, :load, id)[:crd_bus] = crd_bus
        PMD.sol(pm, nw, :load, id)[:cid_bus] = cid_bus
    end
end

function variable_mc_load_current_VoltVar(pm::AbstractExplicitNeutralIVRModelVoltVar; nw=nw_id_default,bounded::Bool=true, report::Bool=true)
    #Get number of phases connected to load and load_ids with VoltVAr control
    connections = Dict(i => load["connections"] for (i,load) in PMD.ref(pm, nw, :load))
    load_ids_VoltVAr = [id for (id,load) in PMD.ref(pm, nw, :load) if load["PV_setpoint"]=="VoltVAr"]
    #Define variables with starting values
    crd = PMD.var(pm, nw)[:crd] = Dict{Int,Any}(i => JuMP.@variable(pm.model,
            [c in 1:connections[i]], base_name="$(nw)_crd_$(i)",
            start = PMD.comp_start_value(PMD.ref(pm, nw, :load, i), "crd_start", c, 0.0)
        ) for i in load_ids_VoltVAr
    )
    cid = PMD.var(pm, nw)[:cid] = Dict{Int,Any}(i => JuMP.@variable(pm.model,
            [c in 1:connections[i]], base_name="$(nw)_cid_$(i)",
            start = PMD.comp_start_value(PMD.ref(pm, nw, :load, i), "cid_start", c, 0.0)
        ) for i in load_ids_VoltVAr
    )

    if bounded
        for (i, l) in PMD.ref(pm, nw, :load)
            cmax = _calc_load_current_max(l, PMD.ref(pm, nw, :bus, l["load_bus"]))
            for (idx,c) in enumerate(connections[i])
                PMD.set_lower_bound(crd[i][c], -cmax[idx])
                PMD.set_upper_bound(crd[i][c],  cmax[idx])
                PMD.set_lower_bound(cid[i][c], -cmax[idx])
                PMD.set_upper_bound(cid[i][c],  cmax[idx])
            end
        end
    end
    #Ensure variables are reported in the solution
    report && IM.sol_component_value(pm, pmd_it_sym, nw, :load, :crd, load_ids_VoltVAr, crd)
    report && IM.sol_component_value(pm, pmd_it_sym, nw, :load, :cid, load_ids_VoltVAr, cid)
end

function constraint_mc_load_current_vars(pm::AbstractExplicitNeutralIVRModelVoltVar, id::Int; nw::Int=nw_id_default, report::Bool=true)
    load = PMD.ref(pm, nw, :load, id)
    bus = PMD.ref(pm, nw, :bus, load["load_bus"])

    PV_setpoint = load["PV_setpoint"]

    if PV_setpoint != "VoltVAr"
        a, alpha,b, beta = PMD._load_expmodel_params(load, bus)
        PMD.constraint_mc_load_current_wye(pm, nw, id, load["load_bus"], load["connections"], a, alpha, b, beta; report=report)
    else
        constraint_mc_load_current_wye_VoltVAr_vars(pm, id, load["load_bus"], load["connections"]; nw=nw, report=report)
    end

end

function constraint_mc_load_current_wye_VoltVAr_vars(pm::AbstractExplicitNeutralIVRModelVoltVar, id::Int, bus::Int, connections::Vector{Int}; nw::Int=nw_id_default, report::Bool=true) # id =  load_id
    vr = PMD.var(pm, nw, :vr, bus)
    vi = PMD.var(pm, nw, :vi, bus)
    crd = PMD.var(pm, nw, :crd, id)
    cid = PMD.var(pm, nw, :cid, id)

    phases = connections[1:end-1]
    n = connections[end] #neutral conductor

    pd = PMD.var(pm, nw, :pd, id)
    qd = PMD.var(pm, nw, :qd, id)

    for (idx, p) in enumerate(phases)
        JuMP.@constraint(pm.model, pd[idx] == (vr[p] - vr[n])*crd[idx] + (vi[p] - vi[n])*cid[idx])
        JuMP.@constraint(pm.model, qd[idx] == -(vr[p] - vr[n])*cid[idx] + (vi[p] - vi[n])*crd[idx])
    end

    crd_bus_n = JuMP.@expression(pm.model, -sum(crd[i] for i in 1:length(phases)))
    cid_bus_n = JuMP.@expression(pm.model, -sum(cid[i] for i in 1:length(phases)))
    PMD.var(pm, nw, :crd_bus)[id] = crd_bus = PMD._merge_bus_flows(pm, [crd..., crd_bus_n], connections)
    PMD.var(pm, nw, :cid_bus)[id] = cid_bus = PMD._merge_bus_flows(pm, [cid..., cid_bus_n], connections)

    if report
        pd_bus = JuMP.NonlinearExpr[]
        qd_bus = JuMP.NonlinearExpr[]
        for (idx, c) in enumerate(connections)
            push!(pd_bus, JuMP.@expression(pm.model, crd_bus[idx]*vr[c] + cid_bus[idx]*vi[c]))
            push!(qd_bus, JuMP.@expression(pm.model, -crd_bus[idx]*vi[c] + cid_bus[idx]*vr[c]))
        end
        PMD.sol(pm, nw, :load, id)[:pd_bus] = JuMP.Containers.DenseAxisArray(pd_bus, connections)
        PMD.sol(pm, nw, :load, id)[:qd_bus] = JuMP.Containers.DenseAxisArray(qd_bus, connections)

        PMD.sol(pm, nw, :load, id)[:crd_bus] = crd_bus
        PMD.sol(pm, nw, :load, id)[:cid_bus] = cid_bus
    end
end

function constraint_mc_load_power(pm::AbstractExplicitNeutralIVRModelVoltVar, id::Int; nw::Int=nw_id_default, report::Bool=true) #TODO needs to be fixed still
    load = PMD.ref(pm, nw, :load, id)
    bus = PMD.ref(pm, nw,:bus, load["load_bus"])
    lambda = PMD.var(pm, nw, :lambda, id)
    breakpoints = load["VV_breakpoints"] #TODO need to add these still
    Q_values = load["VV_Q_values"] #TODO need to add these still
    S = load["S_rating"] #TODO need to add these still
    Q_min = Q_values[1]
    Q_max = Q_values[end]

    
    vr = PMD.var(pm, nw, :vr, bus)
    vi = PMD.var(pm, nw, :vi, bus)
    qd = PMD.var(pm, nw, :qd, id)
    pd = PMD.var(pm, nw, :pd, id)

    phases = connections[1:end-1]
    n = connections[end] #neutral conductor

    for (c, p) in phases  #maybe with epigraph to make convex?
        lambda_c = lambda[c]
        for (k, lambda_c_k) in lambda_c
            JuMP.@constraint(pm.model, lambda_c_k >= 0)
        end
        JuMP.@constraint(pm.model, sum(lambda_c[k] for k in 1:4) == 1)
        #JuMP.@constraint(pm.model, JuMP.SOS2(lambda_c)) #not supported
        JuMP.@constraint(pm.model, (vr[p]-vr[n])^2 + (vi[p]-vi[n])^2 == sum(lambda_c[k]*breakpoints[k] for k in 1:4))
        JuMP.@constraint(pm.model, qd[c] == sum(lambda_c[k]*Q_values[k] for k in 1:4))
        JuMP.@constraint(pm.model, Q_min <= qd[c] <= Q_max)
        JuMP.@constraint(pm.model, pd[c]^2 + qd[c]^2 <= S^2)
    end
end

