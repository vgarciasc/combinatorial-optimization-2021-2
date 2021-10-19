function read_graph_from_txt(filename)
    f = open(filename);
    lines = readlines(f)
    header, edges = lines[1], lines[2:end]

    header = split(header, " ")
    n, m = parse.(Int, header[3:4])

    edges = [split(edge, " ")[2:end] for edge in edges]
    edges = [parse.(Int, edge) for edge in edges]

    n, m, edges
end