using Gurobi, JuMP
using Plots

include("lagrangian_subproblem.jl")
include("plotting_helper.jl")

dist(a, b) = sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2)

K = 2
# pts = [[1 1], [1 4], [1 5], [2 3], [2 5], [3 4], [4 2], [5 1], [6 2], [6 6]]
pts = [[rand() * 10, rand() * 10] for _ in 1:50]
n = length(pts)
N = 1:n
d = [dist(pts[i], pts[j]) for i in N, j in N]

α = 0.5
C = floor(α * n)

println(
    "Lagrangian: ", lagrangian_subproblem(1:n, d, K)
)

model = Model(Gurobi.Optimizer)
@variable(model, x[i in N, j in N] >= 0)
@variable(model, y[j in N], Bin)
@variable(model, w[i in N], Bin)
@constraint(model, all_points_served[i in N], sum(x[i, j] for j in N) == 1 - w[i])
@constraint(model, K_clusters, sum(y[j] for j in N) == K)
@constraint(model, cluster_serves_point[i in N, j in N], x[i, j] <= y[j])
@constraint(model, maximum_outliers, sum(w[i] for i in N) <= C)
@objective(model, Min, sum(d[i, j] * x[i, j] for i in N, j in N))

@time optimize!(model)
# print(model)
# @show(value.(x))
# @show(value.(y))
# @show(value.(w))
# print(objective_value(model))

# pyplot()
plot_2d_solution(pts, x, y, w, K)

# K_centers_idx = [j for j in N if value(y[j]) != 0]
# p = plot()
# for k in 1:K
#     cluster_center = K_centers_idx[k]
#     scatter!([pts[i][1] for i in N if value(x[i, cluster_center]) != 0],
#              [pts[i][2] for i in N if value(x[i, cluster_center]) != 0],
#              label="cluster $k")
# end
# current()