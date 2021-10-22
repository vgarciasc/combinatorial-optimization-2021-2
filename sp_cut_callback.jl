function sp_cut_callback(cb_data)
    global flag_all_explored
    if flag_all_explored
        return
    end

    x_vals = callback_value.(Ref(cb_data), x)

    for c in cliques
        if length(c) >= 3 && sum(x_vals[Array{Int}(c)]) > 1
            con = @build_constraint(sum(x[Array{Int}(c)]) <= 1)
            MOI.submit(model, MOI.UserCut(cb_data), con)
        end
    end

    flag_all_explored = true
end

function sp_cut_callback_one_cut(cb_data)
    global flag_all_explored
    if flag_all_explored
        return
    end

    x_vals = callback_value.(Ref(cb_data), x)

    for xv in 1:length(x_vals)
        if x_vals[xv] > 0 && xv ∉ explored
            res = bron_kerbosch_2_with_adj_list_maximal(adj_list, empty1, adj_include[xv], empty2)
            union!(explored, xv)
            if length(res) >= 3 && sum(x_vals[res]) > 1
                con = @build_constraint(sum(x[res]) <= 1)
                MOI.submit(model, MOI.UserCut(cb_data), con)
                return
            end
        end
    end

    if isempty(setdiff(explored, 1:n))
        flag_all_explored = true
    end
end
