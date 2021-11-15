using Plots

include("lagrangian_subproblem.jl")
include("lagrangian_heuristic.jl")
include("plotting_helper.jl")

function subgradient(n, d, K; 
                     x_start=nothing, y_start=nothing, w_start=nothing, z_upper_start=Inf, z_lower_start=-Inf, u_start=zeros(n),
                     α=0, ϵ=0.1, ρ_min = 0.001, MAX_ITER=1000, THRESHOLD=0.5, θ = 1, heuristic = 1, verbose = false)
    N = 1:n
    hist_u, hist_l = [], []
    is_solved = false
    x_upper, y_upper, w_upper, z_upper = x_start, y_start, w_start, z_upper_start
    z_lower = z_lower_start
    u = u_start
    ρ = 2    
    
    improve = 0
    G_k = zeros(n)
    counter = 0

    for k in 1:MAX_ITER
        counter += 1
        if k == MAX_ITER
            println("Stopping after $k iterations: $k == MAX_ITER!")
            break
        end

        x_u, y_u, w_u, z_lower = lagrangian_subproblem(u, d, K, α)
        x_k, y_k, w_k, z_k = lagrangian_heuristic(n, d, α, x_u, y_u, w_u, heuristic=heuristic)

        if z_k < z_upper
            x_upper, y_upper, w_upper = x_k, y_k, w_k
            z_upper = z_k
            improve = 0
        else
            improve += 1

            if improve >= MAX_ITER/20
                ρ /= 2
                improve = 0

                if ρ < ρ_min
                    println("Stopping after $k iterations: ρ < ρ_min => $ρ < $ρ_min !")
                    break
                end
            end
        end

        G_k = θ * [(1 - w_u[i] - sum(x_u[i, j] for j in N)) for i in N]
        T = ρ * (z_k - z_lower) / sum(g^2 for g in G_k)
        u = u + (1 + ϵ) * T * G_k

        if verbose
            println("====================== ITERATION $k ===================")
            println("\tz_lower: $z_lower")
            println("\tz_k: $z_k")
            println("\tG_$k: $G_k")
            println("\tT_$k: $T")
            println("\tu_$k: $u")
        end

        push!(hist_l, z_lower)
        push!(hist_u, z_upper)

        if z_upper - z_lower < THRESHOLD
            is_solved = true
            println("Stopping after $k iterations: solved!")
            break
        end
    end

    println("\tu_$(counter): $u")
    println("\toptimality gap: $((z_upper - z_lower)/z_upper)")
    println("\tz_k (upper): $z_upper")
    println("\tz_u (lower): $z_lower")
    println("\tclusters: $([i for i in N if y_upper[i] > 0])")
    is_solved, x_upper, y_upper, w_upper, z_upper, z_lower, u, hist_u, hist_l
end

if abspath(PROGRAM_FILE) == @__FILE__
    dist(a, b) = sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2)

    # Gerador de instâncias:
    pts = [round.(rand(2)*10, digits=2) for i in 1:200]

    # Instâncias fáceis:
    # pts = [[1 1], [1 4], [1 5], [2 3], [2 5], [3 4], [4 2], [5 1], [6 2], [6 6]]
    # pts = [[0 5], [1 2], [3 3], [4 3], [7 8], [9 9], [4 2], [10 1], [5 4], [7 6]]

    # Instância dificil:
    # pts = [[5.44, 9.35], [2.02, 8.39], [2.45, 7.81], [3.77, 3.31], [9.88, 6.72], [3.54, 4.03], [2.41, 4.27], [0.56, 3.46], [6.07, 8.86], [0.38, 6.35]]

    n = length(pts)
    d = [dist(pts[i], pts[j]) for i in 1:n, j in 1:n]
    K = 2

    @time is_solved, x, y, w, z_upper, z_lower, hist_u, hist_l = subgradient(n, d, K, α=0, ϵ=0.1, ρ_min=0.001, MAX_ITER=10000, THRESHOLD=0.5, θ=1)

    pyplot()
    p1 = plot(1:length(hist_u), hcat(hist_u, hist_l), label=["upper bound" "lower bound"], ylim=(0, maximum(hist_u)), legend=:bottomright)
    p2 = plot_2d_solution(pts, x, y, w, K)
    plot(p1, p2, layout=(2, 1), title="solution: $(round(z_upper, digits=2))")

    print(pts)
end