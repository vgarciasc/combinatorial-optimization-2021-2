using Plots

include("lagrangian_subproblem.jl")
include("lagrangian_heuristic.jl")
include("plotting_helper.jl")

function subgradient(n, d, K; ρ=0.2, MAX_ITER=1000, THRESHOLD=0.5)
    N = 1:n
    u = zeros(n)

    hist_u, hist_l = [], []
    x_upper, y_upper, z_upper = nothing, nothing, +Inf
    z_lower = -Inf

    for k in 1:MAX_ITER
        x_u, y_u, z_lower = lagrangian_subproblem(u, d, K)
        x_k, y_k, z_k = lagrangian_heuristic(n, d, x_u, y_u)

        if z_k < z_upper
            x_upper, y_upper = x_k, y_k
            z_upper = z_k
        end

        G = [(1 - sum(x_u[i, j] for j in N)) for i in N]
        T = ρ * (z_k - z_lower) / sum(g^2 for g in G)
        # T += (rand() - 0.5)
        u = u + T * G
        
        push!(hist_l, z_lower)
        push!(hist_u, z_upper)

        if z_upper - z_lower < THRESHOLD
            println("Solved after $k iterations!")
            println("\tu_$k: $([round(u[i]; digits=2) for i in N])")
            println("\tz_k (upper): $z_upper")
            println("\tz_u (lower): $z_lower")
            println("\tclusters: $([i for i in N if y_k[i] > 0])")
            break
        end
    end

    x_upper, y_upper, z_upper, z_lower, hist_u, hist_l
end

if abspath(PROGRAM_FILE) == @__FILE__
    dist(a, b) = sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2)

    # Gerador de instâncias:
    # pts = [round.(rand(2)*10, digits=2) for i in 1:10]

    # Instâncias fáceis:
    pts = [[1 1], [1 4], [1 5], [2 3], [2 5], [3 4], [4 2], [5 1], [6 2], [6 6]]
    # pts = [[0 5], [1 2], [3 3], [4 3], [7 8], [9 9], [4 2], [10 1], [5 4], [7 6]]

    # Instância dificil:
    # pts = [[5.44, 9.35], [2.02, 8.39], [2.45, 7.81], [3.77, 3.31], [9.88, 6.72], [3.54, 4.03], [2.41, 4.27], [0.56, 3.46], [6.07, 8.86], [0.38, 6.35]]
    
    n = length(pts)
    d = [dist(pts[i], pts[j]) for i in 1:n, j in 1:n]
    K = 3; ρ = 0.5
    
    x, y, z_upper, z_lower, hist_u, hist_l = subgradient(n, d, K, ρ=ρ, MAX_ITER=500, THRESHOLD=0.5)

    p1 = plot(1:length(hist_u), hcat(hist_u, hist_l), label=["upper bound" "lower bound"], legend=:bottomright)
    p2 = plot_2d_solution(pts, x, y, K)
    plot(p1, p2, layout=(2, 1), title="solution: $(round(z_upper, digits=2)), ρ=$ρ")

    print(pts)
end