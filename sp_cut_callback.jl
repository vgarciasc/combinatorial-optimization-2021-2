function sp_cut_callback(cb_data)

    x_vals = callback_value.(Ref(cb_data), x)

    for c in cliques
        if length(c) >= 3 && sum(x_vals[Array{Int}(c)]) > 1
            con = @build_constraint(sum(x[Array{Int}(c)]) <= 1)
            MOI.submit(model, MOI.UserCut(cb_data), con)
        end
    end
end
