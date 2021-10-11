include("lagrangian_subproblem.jl")

function lagrangian_heuristic(n, d, x_s, y_s)
    N = 1:n
    x, y = copy(x_s), copy(y_s)
    C = [i for i in N if y[i] > 0]
    for i in N
        if sum(x[i, j] for j in N) != 1
            # point not in any cluster, or in many clusters
            # either way, label it to closest cluster only
            closest_c = C[sortperm([d[i, c] for c in C])[1]]
            x[i, :] .= 0
            x[i, closest_c] = 1
        end
    end
    z = sum(d[i, j] * x[i, j] for i in N, j in N)
    return x, y, z
end

if abspath(PROGRAM_FILE) == @__FILE__
    dist(a, b) = sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2)

    K = 3
    pts = [[1 1], [1 4], [1 5], [2 3], [2 5], [3 4], [4 2], [5 1], [6 2], [6 6]]
    n = length(pts)
    N = 1:n
    d = [dist(pts[i], pts[j]) for i in N, j in N]
    u = rand(-10:10, 10)

    x_u, y_u, z_u = lagrangian_subproblem(u, d, K)
    x, y, z = lagrangian_heuristic(n, d, x_u, y_u)

    print("lower bound: $z_u")
    print("upper bound: $z")

    # guaranteeing constraints
    print(sum(y) == K)
    print(all([x[i,j] <= y[j] for i in N, j in N]))
    print(all([sum(x[i, :]) == 1 for i in N]))
end