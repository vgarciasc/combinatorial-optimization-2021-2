using Gurobi, JuMP

m, n = 6, 5
M, N = 1:m, 1:n
f = [4 8 11 7 5]
c = [
    6 2 1 3 5
    4 10 2 6 1
    3 2 4 1 3
    2 0 4 1 4
    1 8 6 2 5
    3 2 4 8 1
]

model = Model(Gurobi.Optimizer)
@variable(model, x[i in M, j in N] >= 0)
@variable(model, y[j in N], Bin)
@constraint(model, all_points_served[i in M], sum(x[i, j] for j in N) == 1)
@constraint(model, cluster_serves_point[i in M, j in N], x[i, j] <= y[j])
@objective(model, Min, sum(c[i, j] * x[i, j] for i in M, j in N) + sum(y[j] * f[j] for j in N))

optimize!(model)
print(model)
@show(value.(x))
@show(value.(y))
print(objective_value(model))