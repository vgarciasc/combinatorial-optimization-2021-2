using Plots

function plot_2d_solution(pts, x, y, K)
    N = 1:length(y)
    C = [i for i in N if y[i] > 0]

    p = Plots.plot(title="")
    for k in 1:K
        scatter!([pts[i][1] for i in N if value(x[i, C[k]]) != 0],
                 [pts[i][2] for i in N if value(x[i, C[k]]) != 0],
                 label="",
                #  label="cluster $k",
                 msw=0,
                 xlabel="x1", ylabel="x2",
                 markersize=6)
    end
    scatter!([pts[c][1] for c in C],
             [pts[c][2] for c in C],
             color="black",
             label="",
             markersize=6)
    
    p
end