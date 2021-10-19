using Gurobi, JuMP

function read_graph_from_txt(filename)
    f = open(filename);
    lines = readlines(f)
    header, edges = lines[1], lines[2:end]

    header = split(header, " ")
    n, m = parse.(Int, header[3:4])

    edges = [split(edge, " ")[2:end] for edge in edges]
    edges = [parse.(Int, edge) for edge in edges]

    n, m, edges
end

n, m, edges = read_graph_from_txt("data/1dc_15.txt")
# n, m, edges = read_graph_from_txt("data/1dc_64.txt")
I = 1:m
J = 1:n
J_i(i) = edges[i]

model = Model(Gurobi.Optimizer)
@variable(model, x[j in J], Bin)
@constraint(model, mutual_exclusion[i in I], sum(x[j] for j in J_i(i)) <= 1)
@objective(model, Max, sum(x[j] for j in J))

optimize!(model)

stable_set = [j for j in J if value(x[j]) > 0]
println("found stable set of size $(length(stable_set)): $stable_set")