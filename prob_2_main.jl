using Gurobi, JuMP

#=
- **User cuts**: Ajuda a fortalecer um MIP removendo soluções fracionárias. Os cortes **não** são necessários para o modelo, mas podem ajudar um MIP a resolver mais rápido.

- MOI.set(model, MOI.RawParameter("Precrush"), 1)

- **Lazy constraints**: São necessárias para o modelo: o modelo estaria incorreto sem essas restrições. Normalmente são utilizadas em modelos que contém um número relativamente grande de restrições que não são facilmente satisfeitas. 

- MOI.set(model, MOI.RawParameter("LazyConstraints"), 1)

Parâmetros: https://www.gurobi.com/documentation/9.1/refman/parameters.html
=#

# 

include("graph_helper.jl")
include("sp_cut_callback.jl")

function solve_set_packing(with_cut=true)
    m, n, edges = G
    I = 1:m
    J = 1:n
    J_i(i) = edges[i]

    global explored = []

    model = Model(Gurobi.Optimizer)

    MOI.set(model, MOI.RawParameter("PreCrush"), 1) # Habilitar cortes do tipo UserCuts
    MOI.set(model, MOI.RawParameter("Cuts"), 0) # Desabilitar cortes
    MOI.set(model, MOI.RawParameter("Presolve"), 0) # Desabilitar presolve
    MOI.set(model, MOI.RawParameter("Heuristics"), 0) # Desabilitar heurísticas
    MOI.set(model, MOI.RawParameter("OutputFlag"), 0) # Desabilitar log  

    @variable(model, x[j in J], Bin)
    @constraint(model, mutual_exclusion[i in I], sum(x[j] for j in J_i(i)) <= 1)

    @objective(model, Max, sum(x[j] for j in J))

    if with_cut
        MOI.set(model, MOI.UserCutCallback(), sp_cut_callback)
    end
    @time optimize!(model)

    stable_set = [j for j in J if value(x[j]) > 10^-10]
    println("X: $([value(x[j]) for j in stable_set])")
    println("found stable set of size $(length(stable_set)): $stable_set")
    println("Explored: $(node_count(model))")
    # println("Conflict Status: $(MOI.get(model, MOI.ConflictStatus()))")
    # plot_graph((m, n, edges), stable_set)
end

for data in [5, 8, 15, 64, 256, 512]
    println("")
    println("-----------------------------------------------------------")
    println("Solving for data: $(data)")
    global G = read_graph_from_txt("data/1dc_$(data).txt")
    println("")
    println("Solving with cuts:")
    solve_set_packing(true)
    println("")
    println("Solving raw:")
    solve_set_packing(false)
end
