include("lagrangian_subgradient.jl")

function check_solvability_percentage(n, iter)
    dist(a, b) = sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2)
    solved_count = 0

    for _ in 1:iter
        pts = [round.(rand(2)*10, digits=2) for i in 1:n]
        d = [dist(pts[i], pts[j]) for i in 1:n, j in 1:n]
        K = 3

        is_solved, x, y, z_upper, z_lower, hist_u, hist_l = subgradient(n, d, K, ϵ=0.1, ρ_min=0.0001, MAX_ITER=2000, THRESHOLD=0.5)

        if is_solved
            solved_count += 1
        else
            println(pts)
        end
    end

    println("percentage of problems solved: $(solved_count / iter)")
end

check_solvability_percentage(100, 100)