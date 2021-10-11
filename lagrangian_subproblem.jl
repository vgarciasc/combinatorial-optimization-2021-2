function lagrangian_subproblem(u, c, K)
    M = length(u)
    N = length(c[1, :])

    y = zeros(M)

    lagrangian_c = [c[i, j] - u[i] for i in 1:M, j in 1:N]

    F = [
        sum(
            min(0, lagrangian_c[i, j]) for i in 1:M
        ) for j in 1:N
    ]

    ind = sortperm(F)[1:K] # Get the indices of the K least F values
    y[ind] .= 1 
    x = [lagrangian_c[i, j] < 0 && y[j] == 1 for i in 1:M, j in 1:N]

    lower_bound = sum(lagrangian_c .* x) + sum(u)

    return x, y, lower_bound
end

function lagrangian_subproblem_JuMP(u, c, K)
    model = Model(Gurobi.Optimizer)

    @variable(model, x[i in N, j in N] >= 0)
    @variable(model, y[j in N], Bin)
    # Relaxing demand constraints
    # @constraint(model, all_points_served[i in N], sum(x[i, j] for j in N) == 1)
    @constraint(model, K_clusters, sum(y[j] for j in N) == K)
    @constraint(model, cluster_serves_point[i in N, j in N], x[i, j] <= y[j])
    @objective(model, Min, sum(d[i, j] * x[i, j] for i in N, j in N) + sum(u[i] * (1 - sum(x[i, j] for j in N)) for i in N))

    optimize!(model)
    x, y = isone.(value.(x)), isone.(value.(y))
    z = objective_value(model)

    x, y, z
end

if abspath(PROGRAM_FILE) == @__FILE__
    using Gurobi, JuMP

    dist(a, b) = sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2)
    K = 3
    pts = [[1 1], [1 4], [1 5], [2 3], [2 5], [3 4], [4 2], [5 1], [6 2], [6 6]]
    n = length(pts)
    N = 1:n
    d = [dist(pts[i], pts[j]) for i in N, j in N]
    u = rand(-10:10, 10)

    x, y, z = lagrangian_subproblem_JuMP(u, d, K)
    x_u, y_u, z_u = lagrangian_subproblem(u, d, K)
    print(x_u == x)
    print(y_u == y)
    print(abs(z_u - z) < 0.01)
end