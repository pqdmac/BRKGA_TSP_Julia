

# TO DO List

function decoder(population, dist)

    pop = deepcopy(population)
    len_pop = length(pop)
    pop_sol = TspData[]
    
    for i in 1:len_pop
        route, fo_total =  Cheapest_Insertion(dist, pop[i])       
        aux_TspData = TspData(pop[i], route, fo_total)
        push!(pop_sol, aux_TspData)
    end

    return pop_sol
end

function Cheapest_Insertion(dist, individual)

    NR = sortperm(individual)
    NR_aux = deepcopy(NR)

    len_initial_route = 2
    route = NR[1:len_initial_route]

    fo = 0
    size = length(route)

    for i in 1:size
        vi = route[i]
        if i == size
            j = 1
        else
            j = i + 1
        end       
        vj = route[j]
        fo += dist[vi, vj]
    end

    deleteat!(NR_aux, 1:len_initial_route)

    min_vk, min_j = 0, 0
    z_min = Inf
    len_NR_aux = length(NR_aux)

    for k in 1:len_NR_aux
        vk = NR_aux[k]
        z_min = Inf
        size_route = length(route)      # Tamanho da rota considerada
        for i in 1:size_route
            vi = route[i]
            if i == size_route
                j = 1
            else
                j = i + 1
            end
            vj = route[j]          
            z = dist[vi, vk] + dist[vk, vj] - dist[vi, vj]
            if z < z_min
                z_min = z
                min_j = j
                min_vk = vk
            end
        end
        insert!(route, min_j, min_vk)
        fo += z_min
    end

    fo = round(fo, digits=2)

    return route, fo

end