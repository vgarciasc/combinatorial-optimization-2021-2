function sp_cut_callback(cb_data)
    # println("Cut")

    x_vals = callback_value.(Ref(cb_data), x)

    R = []
    X = []
    P = [i for i in 1:length(x_vals) if x_vals[i] > 0]

    for ex in explored
        if P ⊆ ex
            # println("Exp: $(explored)")
            return
        end
    end
    append!(explored, [P])
    # println(x_vals)
        
    # println("P: $(P)")
    cliques = bron_kerbosch_2(G, R, P, X)
    for R in cliques
        if length(R) >= 3      
            # println("Clique: $(Array{Int}(R))")
            con = @build_constraint(sum(x[Array{Int}(R)]) <= 1)
            # if sum(x_vals[Array{Int}(R)]) < 1
            #     println("Sum does not overcome 0")
            # end
            # println("Adding $(con)")
            MOI.submit(model, MOI.UserCut(cb_data), con)
        end
    end
end
