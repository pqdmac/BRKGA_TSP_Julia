
using DelimitedFiles

# TO DO List

function ReadData(datafilename)

    data = DelimitedFiles.readdlm(datafilename, ' ')

    Coord_X = hcat(data[7:end,2]...)       # Coordinate X
    Coord_Y = hcat(data[7:end,3]...)       # Coordinate Y
    N = length(Coord_X)                    # Number of nodes
    
    distance_matrix = [sqrt((Coord_X[i] - Coord_X[j])^2 + 
                        (Coord_Y[i] - Coord_Y[j])^2) for i in 1:N, j in 1:N]

    return Coord_X, Coord_Y, distance_matrix, N
end

function ReadBKS(datafilename)
    
    data = DelimitedFiles.readdlm(datafilename, '\t')

    path_instances = hcat(data[1:end,1]...)
    BKS_instances = hcat(data[1:end,2]...)

    return path_instances, BKS_instances
end

function check_route(route, dist, N)
    cost = 0
    for i in 1:N
        vi = route[i]
        if i == N
            j = 1
        else
            j = i + 1
        end
        vj = route[j]
        cost += dist[vi, vj]
    end
    cost = round(cost, digits=2)   
    println("\nCost verified: $cost")
    return
end

function Write_Data(probname, s0, fo0, rpd0, t0, tt)
    len_route = length(s0)
    arquivo = open("Results_per_instance/$probname.txt", "w")
    for i in 1:len_route

        write(arquivo, string(s0[i]))

        if i == len_route
            write(arquivo, "\n")
        else
            write(arquivo, " ")
        end

    end

    write(arquivo, "\nBest FO\n")
    write(arquivo, string(fo0))
    write(arquivo, "\n\nBest RPD\n")
    write(arquivo, string(rpd0))
    write(arquivo, "\n\nTime to best FO\n")
    write(arquivo, string(t0))
    write(arquivo, "\n\nTotal Time\n")
    write(arquivo, string(tt))
    close(arquivo)

    return
end



