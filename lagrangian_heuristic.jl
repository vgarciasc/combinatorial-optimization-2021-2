include("plotting_helper.jl")
include("lagrangian_subproblem.jl")

function lagrangian_heuristic(n, d, α, x_s, y_s, w_s; heuristic=1)
    if α > 0 && heuristic != 1
        return lagrangian_heuristic_outliers(n, d, α, x_s, y_s, w_s)
    end

    N = 1:n
    x, y, w = copy(x_s), copy(y_s), copy(w_s)
    P = [i for i in N if y[i] > 0]
    for i in N
        x[i, :] .= 0
        if w[i] == 0 # point not an outlier
            # point not in any cluster, or in many clusters
            # either way, label it to closest cluster only
            closest_p = P[sortperm([d[i, p] for p in P])[1]]
            x[i, closest_p] = 1
        end
    end
    z = sum(d[i, j] * x[i, j] for i in N, j in N)

    x, y, w, z
end

function lagrangian_heuristic_outliers(n, d, α, x_s, y_s, w_s)
    N = 1:n
    x, y, w = copy(x_s), copy(y_s), copy(w_s)
    P = [i for i in N if y[i] > 0]
    C = floor(Int, α * n)
    dists = []
    for i in N
        x[i, :] .= 0
        # point not in any cluster, or in many clusters
        # either way, label it to closest cluster only
        closest_cluster = P[sortperm([d[i, p] for p in P])[1]]
        x[i, closest_cluster] = 1
        push!(dists, (i, closest_cluster, d[i, closest_cluster]))
    end
    w .= 0
    for (i, j, dist) in sort(dists, by=a->-a[3])[1:C]
        x[i, :] .= 0
        w[i] = 1
    end

    z = sum(d[i, j] * x[i, j] for i in N, j in N)

    x, y, w, z
end

if abspath(PROGRAM_FILE) == @__FILE__
    dist(a, b) = sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2)
    K = 2
    # pts = [[rand() * 10, rand() * 10] for _ in 1:500]
    n = length(pts)
    N = 1:n
    d = [dist(pts[i], pts[j]) for i in N, j in N]
    u = rand(n) * 100
    α = 0.2

    x_u, y_u, w_u, z_u = lagrangian_subproblem(u, d, K, α)
    x, y, w, z = lagrangian_heuristic(n, d, α, x_u, y_u, w_u)
    x2, y2, w2, z2 = lagrangian_heuristic_outliers(n, d, α, x_u, y_u, w_u)

    print("lower bound: $z_u")
    print("upper bound 1: $z")
    print("upper bound 2: $z2")

    pyplot()
    p1 = plot_2d_solution_multiple_clusters(pts, x_u, y_u, w_u, K)
    p2 = plot_2d_solution(pts, x, y, w, K)
    p3 = plot_2d_solution(pts, x2, y2, w2, K)
    
    plot(p1, p2, p3, layout=(1, 3))
    plot!(size=(1200, 400))

    # guaranteeing constraints
    print(sum(y) == K)
    print(all([x[i,j] <= y[j] for i in N, j in N]))
    print(all([sum(x[i, :]) == 1 - w[i] for i in N]))
    
    print(sum(y2) == K)
    print(all([x2[i,j] <= y2[j] for i in N, j in N]))
    print(all([sum(x2[i, :]) == 1 - w2[i] for i in N]))
end