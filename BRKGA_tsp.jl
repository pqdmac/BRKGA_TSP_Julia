
using CSV
using DataFrames
using StatsBase

# TO DO List

include("ReadFile.jl")
include("Decoder.jl")

ti = time()

mutable struct TspData
    rk:: Vector{Float64}
    sol:: Vector{Int64}
    fo:: Float64
end

function crossover(pop_decoded, size_pop, int_pe, int_pm, rhoe, dist, N)

    id_eliteParent = rand(1:int_pe)
    eliteParent = pop_decoded[id_eliteParent]
    id_non_eliteParent = rand(1:size_pop)               # Consider the entire population
    non_eliteParent = pop_decoded[id_non_eliteParent]
    rand_aux = round.(rand(N), digits = 3)
    rk_aux = Float64[]

    for i in 1:N
        if rand_aux[i] < rhoe
            push!(rk_aux, eliteParent.rk[i])
        else
            push!(rk_aux, non_eliteParent.rk[i])
        end
    end

    route, fo =  Cheapest_Insertion(dist, rk_aux)
    return rk_aux, route, fo
end

function generate_offspring(pop_decoded, size_pop, int_pe, int_pm, rhoe, dist, N)
    
    pop_gen = TspData[]
    for i in (int_pe + 1 ):(size_pop - int_pm)

        rk_gen, route_gen, fo_gen = crossover(pop_decoded, size_pop, int_pe, int_pm, rhoe, dist, N)
        aux_TspData = TspData(rk_gen, route_gen, fo_gen)
        push!(pop_gen, aux_TspData)
    end

    return pop_gen
end

function generate_mutants(int_pm, N, dist)

    pop_m = [round.(rand(N); digits = 3) for i in 1:int_pm]
    pop_m_decoded = decoder(pop_m, dist)

    return pop_m_decoded
end

function run_BRKGA(pe, pm, rhoe, size_pop, N, max_Iter, dist)

    time2best = 0
    int_pe = Int(round(pe * size_pop))
    int_pm = Int(round(pm * size_pop))

    # Create vectors of random keys
    pop = [round.(rand(N); digits = 3) for i in 1:size_pop]
    pop_decoded = decoder(pop, dist)
    sort!(pop_decoded, by = v -> v.fo)

    best_fo = pop_decoded[1].fo        
    best_sol = pop_decoded[1].sol
    # check_route(best_sol, dist, N)

    Iter = 1
    while Iter < max_Iter

        # println("\nIteration: $Iter")
        new_pop = pop_decoded[1:int_pe]
        pop_mut = generate_mutants(int_pm, N, dist)
        append!(new_pop, pop_mut)
        pop_offspring = generate_offspring(pop_decoded, size_pop, int_pe, int_pm, rhoe, dist, N)
        append!(new_pop, pop_offspring)
        sort!(new_pop, by = v -> v.fo)
        fo = new_pop[1].fo

        if fo < best_fo
            time2best = time()
            best_fo = fo
            # println("\n##########################\nIteration: $Iter")
            # println("Partial_FO: $best_fo")
            best_sol = new_pop[1].sol
        end
        pop_decoded = deepcopy(new_pop)
        Iter += 1
    end

    return best_fo, best_sol, time2best
end


global BRKGA_results = DataFrame(Instancia = String[], BKS = Float64[], 
    Best_FO = Float64[], Mean_FO = Float64[], Worse_FO = Float64[], Best_RPD = Float64[],
    Mean_RPD = Float64[], Time_to_best_FO = Float64[], Mean_Time_to_best_FO = Float64[],
    Run_time_best_FO = Float64[], Mean_run_time = Float64[]
)

# datafile = "Instances\\list_40_instances_plus_BKS.txt"
datafile = "Instances\\list_20_instances_plus_BKS.txt"

inst, BKS = ReadBKS(datafile)

nr_files = length(BKS)
nr_runs = 3

println("Running...")
    
for i in 1:nr_files

    if i == 5
        break
    end

    fos = Float64[]
    rpds = Float64[]
    best_times = Float64[]
    times = Float64[]

    solutions = []

    instance_results = DataFrame(BKS = Float64[], FO_BRKGA = Float64[],  RPD_BRKGA = Float64[],
                                    Time_to_best_FO = Float64[], Total_time = Float64[])

    instance = inst[i]

    probname = instance[15:end-4]
    println("\nInstance:\t $probname")

    size_pop = 600      # population size
    pe = 0.20           # elite fraction
    pm = 0.05           # mutants fraction
    rhoe = 0.7          # probability to inherit genes from elite
    max_Iter = 600      # number of generations

    for r in 1:nr_runs

        initial_time_BRKGA = time()

        X, Y, dist, N = ReadData(instance)

        fo_BRKGA, sol_BRKGA, time2best = run_BRKGA(pe, pm, rhoe, size_pop, N, max_Iter, dist)

        final_time_BRKGA = time()
        time_BRKGA = round(final_time_BRKGA - initial_time_BRKGA, digits=2)
        rpd = round(100 * (fo_BRKGA - BKS[i])/BKS[i], digits = 2)
        time2best = round(time2best - initial_time_BRKGA, digits=2)

        push!(instance_results, (BKS[i], fo_BRKGA, rpd, time2best, time_BRKGA))

        push!(solutions, sol_BRKGA)
        push!(fos, fo_BRKGA)
        push!(rpds, rpd)
        push!(best_times, time2best)
        push!(times, time_BRKGA)

    end

    println(instance_results)

    # CSV.write("Results_per_instance/$probname.csv", instance_results)

    best_fo = minimum(fos)
    mean_fo = round(mean(fos), digits=2)
    worse_fo = maximum(fos)

    # Searching for an element in a 1D ARRAY 
    sch = best_fo  
    positionArray = indexin( sch, fos ) 

    best_sol = solutions[positionArray[1]]

    time2bestFO = best_times[positionArray[1]]
    mean_time2bestFO = round(mean(best_times), digits=2)

    run_time_best_FO = times[positionArray[1]]
    mean_run_time = round(mean(times), digits=2)

    total_time_best_fo = times[positionArray[1]]

    best_rpd = minimum(rpds)
    mean_rpd = round(mean(rpds), digits=2)


    # Write_Data(probname, best_sol, best_fo, best_rpd, time2bestFO, total_time_best_fo)


    push!(BRKGA_results, (probname, BKS[i], best_fo, mean_fo, worse_fo, 
                        best_rpd, mean_rpd, time2bestFO, mean_time2bestFO,
                        run_time_best_FO, mean_run_time))
    # CSV.write("BRKGA_results.csv", BRKGA_results)

end

println("\n")
println(BRKGA_results)
println("\n")
# CSV.write("BRKGA_results.csv", BRKGA_results)

tf = time()
general_time = round(tf - ti, digits = 2)
println("\nTotal_run_time: $general_time")


# # ############################




