using Statistics, LinearAlgebra, StatsBase
using DataFrames, CSV, Random
using Plots

include("lagrangian_subgradient.jl")

# input: a (n, p) matrix X, where p >= 1 is the number of attributes
function calc_distance_matrix(X)
    n = size(X, 1)
    d = zeros(n, n)
    X = hcat(pts...)'

    dt = fit(UnitRangeTransform, X, dims=1)
    StatsBase.transform(dt, X)

    for i in 1:n
        for j in 1:n
            d[i, j] = norm(X[i, :] - X[j, :])
        end
    end
    d, dt
end

function read_points_from_file(filename)
    f = open(filename)
    lines = readlines(f)
    header, pts = lines[1], lines[2:end]

    header = split(header, " ")
    if header[end] == "label"
        _, K = parse.(Int, header[1:end-1])
    else
        _, K = parse.(Int, header[1:end])
    end
    
    pts = [split(pt, " ") for pt in pts]
    pts = [parse.(Float64, pt) for pt in pts]
    
    labels = []
    if header[end] == "label"
        labels = [pt[end] for pt in pts]
        pts = [pt[1:end-1] for pt in pts]
    end

    K, pts, labels
end

function check_solvability_percentage(n, iter)
    solved_count = 0

    for _ in 1:iter
        pts = [round.(rand(2)*10, digits=2) for i in 1:n]
        d = [norm(pts[i] - pts[j]) for i in 1:n, j in 1:n]

        is_solved, x, y, z_upper, z_lower, _, hist_u, hist_l = subgradient(n, d, K, ϵ=0.1, ρ_min=0.0001, MAX_ITER=2000, THRESHOLD=0.5)

        if is_solved
            solved_count += 1
        else
            println(pts)
        end
    end

    println("percentage of problems solved: $(solved_count / iter)")
end

function run_basic_clustering(filename)
    K, pts = read_points_from_file(filename)
    n = length(pts)
    @time d, _ = calc_distance_matrix(pts)

    @time is_solved, x, y, w, z_upper, z_lower, _, hist_u, hist_l = subgradient(n, d, K, α=0, ϵ=0.1, ρ_min=0.0001, MAX_ITER=20000, THRESHOLD=0.5)

    scatter([pt[1] for pt in pts], [pt[2] for pt in pts], label="", msw=0, color=:lightgrey, markersize=5, xlabel="x1", ylabel="x2")
    p1 = plot(1:length(hist_u), hcat(hist_u, hist_l), 
              label=["upper bound" "lower bound"], 
              legend=:topright, 
              xlabel="iterations",
              ylabel="value",
              ylim=(0, maximum(hist_u))
            #   ylim=(quantile(hist_l, 0.2), maximum(hist_u))
              )
    p2 = plot_2d_solution(pts, x, y, w, K)
    plot(p1, p2, layout=(2, 1), title="solution: $(round(z_upper, digits=2))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    filename = "data/clustering/dim032.txt"
    K, pts = read_points_from_file(filename)
    n = length(pts)
    
    @time d, _ = calc_distance_matrix(pts)

    @time is_solved, x, y, w, z_upper, z_lower, u, hist_u, hist_l = subgradient(n, d, K, α=0.2, ϵ=0, 
        ρ_min=0.0001, MAX_ITER=10000, THRESHOLD=0.5)

    # @time is_solved, x, y, w, z_upper, z_lower, u, hist_u, hist_l = subgradient(n, d, K, α=0.2, ϵ=0.1, ρ_min=0.0001, MAX_ITER=10000, THRESHOLD=0.5) 
    # println("Optimality gap: ($z_upper - $z_lower) / $z_upper = $((z_upper - z_lower) / z_upper)")

    # x_s, y_s, w_s, z_upper = lagrangian_heuristic(n, d, 0, x, y, zeros(n))

    # @time is_solved, x_k, y_k, w_k, z_upper, z_lower, u, hist_u, hist_l = subgradient(n, d, K, α=0.2, ϵ=0, ρ_min=0.0001, MAX_ITER=5000, THRESHOLD=0.5,
    #     x_start = x_s, y_start = y_s, w_start = w_s, z_upper_start = z_upper, z_lower_start = z_lower, u_start = u)
    # println("Optimality gap: ($z_upper - $z_lower) / $z_upper = $((z_upper - z_lower) / z_upper)")

    optimality_gap = (z_upper - z_lower) / z_upper

    pyplot()
    scatter([pt[1] for pt in pts], [pt[2] for pt in pts], label="", msw=0, color=:lightgrey, markersize=5, xlabel="x1", ylabel="x2", title="")
    p1 = plot(1:length(hist_u), hcat(hist_u, hist_l), 
              label=["upper bound" "lower bound"], 
              legend=:topright, 
              xlabel="iterations",
              ylabel="value",
            #   ylim=(0, maximum(hist_u))
              ylim=(quantile(hist_l, 0.2), maximum(hist_u))
              )
    p2 = plot_2d_solution(pts, x, y, w, K)
    plot(p1, p2, layout=(1, 2), title="solution: $(round(z_upper, digits=2)), gap: $(round(optimality_gap, digits=2))")
    plot!(size=(1200, 400))
end