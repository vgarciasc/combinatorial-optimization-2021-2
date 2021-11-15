using Gurobi, JuMP

include("plotting_helper.jl")
    
function lagrangian_subproblem(u, c, K, α=0)
    m = length(u)
    n = length(c[1, :])
    C = floor(Int, α * n)

    y = zeros(m)
    w = zeros(m)

    lagrangian_c = [c[i, j] - u[i] for i in 1:m, j in 1:n]

    F = [
        sum(
            min(0, lagrangian_c[i, j]) for i in 1:m
        ) for j in 1:n
    ]

    ind = sortperm(F)[1:K] # Get the indices of the K least F values
    y[ind] .= 1 
    x = [lagrangian_c[i, j] < 0 && y[j] == 1 for i in 1:m, j in 1:n]

    if α > 0
        ind = sortperm(u, rev=true)[1:C] # Get the indices of the C largest u values
        for i in ind
            w[i] = (u[i] > 0) ? 1 : 0
        end
    end

    lower_bound = sum(lagrangian_c .* x) + sum(u) - sum(u[i] * w[i] for i in 1:m)

    return x, y, w, lower_bound
end

function lagrangian_subproblem_JuMP(u, c, K, α=0)
    n = length(c[1, :])
    N = 1:n
    C = floor(α * n)

    model = Model(Gurobi.Optimizer)

    @variable(model, x[i in N, j in N] >= 0)
    @variable(model, y[j in N], Bin)
    @variable(model, w[i in N], Bin)

    # Relaxing demand constraints
    # @constraint(model, all_points_served[i in N], sum(x[i, j] for j in N) == 1)
    @constraint(model, cluster_serves_point[i in N, j in N], x[i, j] <= y[j])
    @constraint(model, K_clusters, sum(y[j] for j in N) == K)
    @constraint(model, max_C_outliers, sum(w[i] for i in N) <= C)

    @objective(model, Min, sum(c[i, j] * x[i, j] for i in N, j in N) + sum(u[i] * (1 - w[i] - sum(x[i, j] for j in N)) for i in N))

    optimize!(model)
    
    x, y = isone.(value.(x)), isone.(value.(y))
    z = objective_value(model)
    w = (value.(w))

    x, y, w, z
end

if abspath(PROGRAM_FILE) == @__FILE__
    dist(a, b) = sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2)
    K = 2
    pts = [[rand() * 10, rand() * 10] for _ in 1:10]
    n = length(pts)
    N = 1:n
    d = [dist(pts[i], pts[j]) for i in N, j in N]
    u = rand(n) * 100
    α = 0.5
    
    # scatter([pt[1] for pt in pts], [pt[2] for pt in pts], label="", msw=0, color=:lightgrey, markersize=5, xlabel="x1", ylabel="x2")
    
    x, y, w, z = lagrangian_subproblem_JuMP(u, d, K, α)
    x_u, y_u, w_u, z_u = lagrangian_subproblem(u, d, K, α)
    
    # plot_2d_solution(pts, x, y, w, K)

    print(x_u == x)
    print(y_u == y)
    print(w_u == w)
    print(abs(z_u - z) < 0.01)
end