using PowerModels
using Ipopt
using JuMP
PowerModels.silence()


function BPFRV(file_name)
    data = PowerModels.parse_file(file_name)

    PowerModels.calc_thermal_limits!(data)
    PowerModels.standardize_cost_terms!(data, order=2)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

    model = Model(Ipopt.Optimizer)

    set_optimizer_attribute(model, "print_level", 0)

    @variable(model, vr[i in keys(ref[:bus])], start=1.0)
    @variable(model, vi[i in keys(ref[:bus])], start=0.0)
    @variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])
    @variable(model, ref[:gen][i]["qmin"] <= qg[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"])

    @variable(model, -ref[:branch][l]["rate_a"] <= p[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])
    @variable(model, -ref[:branch][l]["rate_a"] <= q[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])

    @objective(model, Min,
        sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen])
    )

    for (i,bus) in ref[:ref_buses]
        @constraint(model, vi[i] == 0)
    end

    for (i, bus) in ref[:bus]
        @constraint(model, bus["vmin"]^2 <= vr[i]^2 + vi[i]^2)
        @constraint(model, vr[i]^2 + vi[i]^2 <= bus["vmax"]^2)
    end

    for (i,bus) in ref[:bus]
        # Build a list of the loads and shunt elements connected to the bus i
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        # Active power balance at node i
        @constraint(model, 
            sum(pg[g] for g in ref[:bus_gens][i]) - sum(load["pd"] for load in bus_loads) ==
            sum(p[a] for a in ref[:bus_arcs][i])  + sum(shunt["gs"] for shunt in bus_shunts) * (vr[i]^2 + vi[i]^2)
        )

        @constraint(model,
            sum(qg[g] for g in ref[:bus_gens][i]) - sum(load["pd"] for load in bus_loads) == 
            sum(q[a] for a in ref[:bus_arcs][i]) - sum(shunt["bs"] for shunt in bus_shunts) * (vr[i]^2 + vi[i]^2)
        )
    end

    for (i, branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])

        p_fr = p[f_idx]
        q_fr = q[f_idx]
        p_to = p[t_idx]
        q_to = q[t_idx]

        vr_fr = vr[branch["f_bus"]]
        vr_to = vr[branch["t_bus"]]
        vi_fr = vi[branch["f_bus"]]
        vi_to = vi[branch["t_bus"]]

        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        tm = branch["tap"]^2 

        g_net_fr = (g+g_fr)/tm
        b_net_fr = (b+b_fr)/tm
        g_net_to = (g+g_to)
        b_net_to = (b+b_to)

        G_fr = (-g*tr+b*ti)/tm
        B_fr = (-b*tr-g*ti)/tm

        G_to = (-g*tr-b*ti)/tm
        B_to = (-b*tr+g*ti)/tm

        @constraint(model, p_fr == g_net_fr * (vr_fr^2 + vi_fr^2) + G_fr * (vr_fr * vr_to + vi_fr * vi_to) + B_fr * (vi_fr * vr_to - vr_fr * vi_to))
        @constraint(model, q_fr == -b_net_fr * (vr_fr^2 + vi_fr^2) + G_fr * (vi_fr * vr_to - vr_fr * vi_to) - B_fr * (vr_fr * vr_to + vi_fr * vi_to))

        @constraint(model, p_to == g_net_to * (vr_to^2 + vi_to^2) + G_to * (vr_to * vr_fr + vi_to * vi_fr) + B_to * (vi_to * vr_fr - vr_to * vi_fr))
        @constraint(model, q_to == -b_net_to * (vr_to^2 + vi_to^2) + G_to * (vi_to * vr_fr - vr_to * vi_fr) - B_to * (vr_to * vr_fr + vi_to * vi_fr))

        # Apparent power limit, from side and to side
        @constraint(model, p_fr^2 + q_fr^2 <= branch["rate_a"]^2)
        @constraint(model, p_to^2 + q_to^2 <= branch["rate_a"]^2) 
    end

    optimize!(model)

    return model
end


# m = BPFRV("/home/jinxin/xjx/SRIBD/ACOPF/code/PowerModelsAnnex.jl/src/model/case9.m")
# println(solve_time(m))
# println(objective_value(m))
