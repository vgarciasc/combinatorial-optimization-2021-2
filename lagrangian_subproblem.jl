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
    x = [lagrangian_c[i, j] <= 0 && y[j] == 1 for i in 1:M, j in 1:N]

    lower_bound = sum(lagrangian_c .* x) + sum(u)

    return x, y, lower_bound
end