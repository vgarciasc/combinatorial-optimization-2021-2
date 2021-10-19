using LightGraphs, GraphRecipes, Plots

function read_graph_from_txt(filename)
    f = open(filename);
    lines = readlines(f)
    header, edges = lines[1], lines[2:end]

    header = split(header, " ")
    n, m = parse.(Int, header[3:4])

    edges = [split(edge, " ")[2:end] for edge in edges]
    edges = [parse.(Int, edge) for edge in edges]

    m, n, edges
end

function neighbors(G, v)
    _, _, edges = G
    [e[1] for e in edges if e[2] == v] ∪ [e[2] for e in edges if e[1] == v] 
end

function bron_kerbosch_1(G, R, P, X, cliques=[])
    if isempty(P) && isempty(X)
        cliques = cliques ∪ [R]
    end
    
    for v in P
        N_v = neighbors(G, v)
        cliques = cliques ∪ bron_kerbosch_1(G, R ∪ [v], P ∩ N_v, X ∩ N_v)
        P = setdiff(P, [v])
        X = X ∪ [v]
    end

    cliques
end

function bron_kerbosch_2(G, R, P, X, cliques=[])
    if isempty(P) && isempty(X)
        cliques = cliques ∪ [R]
    end
    
    u = (P ∪ X)[begin]

    for v in setdiff(P, neighbors(G, u))
        N_v = neighbors(G, v)
        cliques = cliques ∪ bron_kerbosch_1(G, R ∪ [v], P ∩ N_v, X ∩ N_v)
        P = setdiff(P, [v])
        X = X ∪ [v]
    end

    cliques
end

function plot_graph(G, stable_set=[])
    m, n, edges = G
    
    A = zeros(n, n)
    for edge in edges
        i, j = edge
        A[i, j] = 1
        A[j, i] = 1
    end

    graphplot(SimpleGraph(A),
              markercolor = :lightgrey,
              names = [(i in stable_set ? "$i*" : i) for i in 1:n],
              markersize = 0.2,
              fontsize = 10,
              curves = false,
              nodeshape = [(i in stable_set ? :rect : :circle) for i in 1:n])
end

if abspath(PROGRAM_FILE) == @__FILE__
    G = read_graph_from_txt("data/1dc_15.txt")
    m, n, edges = G
    @time cliques = bron_kerbosch_2(G, [], 1:n, [])
    plot_graph(G)
end