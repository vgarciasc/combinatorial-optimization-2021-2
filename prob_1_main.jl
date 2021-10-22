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

function shuffle_data(pts, labels)
    order = shuffle(1:length(pts))
    pts[order, :], labels[order, :]
end

function check_solvability_percentage(n, iter)
    solved_count = 0

    for _ in 1:iter
        pts = [round.(rand(2)*10, digits=2) for i in 1:n]
        d = [norm(pts[i] - pts[j]) for i in 1:n, j in 1:n]

        is_solved, x, y, z_upper, z_lower, hist_u, hist_l = subgradient(n, d, K, ϵ=0.1, ρ_min=0.0001, MAX_ITER=2000, THRESHOLD=0.5)

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

    @time is_solved, x, y, z_upper, z_lower, hist_u, hist_l = subgradient(n, d, K, ϵ=0.1, ρ_min=0.0001, MAX_ITER=10000, THRESHOLD=0.5)
    println("Optimality gap: $((z_upper - z_lower) / z_upper)")

    pyplot()
    scatter([pt[1] for pt in pts], [pt[2] for pt in pts], label="", msw=0, color=:lightgrey, markersize=5, xlabel="x1", ylabel="x2")
    p1 = plot(1:length(hist_u), hcat(hist_u, hist_l), 
              label=["upper bound" "lower bound"], 
              legend=:topright, 
              xlabel="iterations",
              ylabel="value",
            #   ylim=(0, maximum(hist_u))
              ylim=(quantile(hist_l, 0.2), maximum(hist_u))
              )
    p2 = plot_2d_solution(pts, x, y, K)
    # plot(p1, p2, layout=(2, 1), title="solution: $(round(z_upper, digits=2))")
end

function run_clustering_for_classification()
    filename = "data/clustering/movement_libras.txt"
    K, pts, labels = read_points_from_file(filename)
    pts, labels = shuffle_data(pts, labels)
    possible_labels = unique(labels)
    
    fold_qt = 10
    fold_size = trunc(Int, length(pts) / fold_qt)

    accuracies = []

    for k in 1:fold_qt
        data = hcat(hcat(pts...)', hcat(labels...)')

        test = data[max(((k - 1) * fold_size), 1) : (k * fold_size), :]
        training_1 = k == 1       ? [] : data[begin : ((k - 1) * fold_size), :]
        training_2 = k == fold_qt ? [] : data[(k * fold_size) : end, :]
        training = isempty(training_1) ? training_2 : (isempty(training_2) ? training_1 : vcat(training_1, training_2))
        
        training_pts, training_lbs = training[:, 1:end-1], training[:, end]
        test_pts, test_lbs = test[:, 1:end-1], test[:, end]

        n = size(training_pts, 1)
        d, dt = calc_distance_matrix(training_pts)

        @time is_solved, x, y, z_upper, z_lower, hist_u, hist_l = subgradient(n, d, K, ϵ=0.1, ρ_min=0.0001, MAX_ITER=10000, THRESHOLD=0.5)
        
        cluster_centers = [i for i in 1:n if y[i] == 1]
        pts_per_cluster = [[(i, "_", training_lbs[i]) for i in 1:n if x[i, j] == 1] for j in cluster_centers]
        
        classification_matrix = [[length(filter(a -> a[3] == label, pts_per_cluster[j])) / length(pts_per_cluster[j]) for label in 1:K] for j in 1:K]
        classification_matrix = hcat(classification_matrix...)'
        cluster_label_pairings = []

        for _ in 1:K
            assigned_clusters = [cl for (cl, lb) in cluster_label_pairings]
            assigned_labels = [lb for (cl, lb) in cluster_label_pairings]

            best_pairing = convert(Tuple, argmax(classification_matrix))
            counter = 0
            while best_pairing[1] in assigned_clusters || best_pairing[2] in assigned_labels
                i, j = best_pairing
                classification_matrix[i, j] = -1
                
                best_pairing = convert(Tuple, argmax(classification_matrix))
                counter += 1
                if counter > 500
                    println("SOCORRO!")
                    break
                end
            end

            i, j = best_pairing
            push!(cluster_label_pairings, best_pairing)
            sort!(cluster_label_pairings, by = a -> a[1])
            println("added $best_pairing to $cluster_label_pairings")
            
            classification_matrix[i, :] = zeros(K)
            classification_matrix[:, j] = zeros(K)
        end
        
        predicted_clusters = []
        for cluster in 1:K
            assigned_label = [assigned_label for (cl, assigned_label) in cluster_label_pairings if cl == cluster][1]
            push!(predicted_clusters, [(idx, pt, label, (label == assigned_label ? 1 : 0)) for (idx, pt, label) in pts_per_cluster[cluster]])
        end

        correct = 0
        test_pts = StatsBase.transform(dt, test_pts)
        for i in 1:size(test_pts, 1)
            pt, lb = test_pts[i, :], test_lbs[i]

            best_dist = 10000
            assigned_cluster = -1
            for k in 1:K
                cdist = norm(StatsBase.transform(dt, pt') - StatsBase.transform(dt, training_pts[cluster_centers[k], :]'))
                if cdist < best_dist
                    assigned_cluster = k
                    best_dist = cdist
                end
            end
            
            _, assigned_label = filter(a -> a[1] == assigned_cluster, cluster_label_pairings)[1]
            if lb == assigned_label
                correct += 1
            end

            cluster_label_pairings
        end
        
        accuracy = sum([sum([correct for (_,_,_,correct) in predicted_clusters[k]]) for k in 1:K]) / n
        # accuracy = correct / size(test_pts, 1)
        push!(accuracies, accuracy)
    end

    println("mean accuracy: $(mean(accuracies))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_basic_clustering("data/clustering/dim032.txt")
end