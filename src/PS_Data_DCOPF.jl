"""
PS_Data_DCOPF.jl

This module extracts power system data for Direct Current Optimal Power Flow (DCOPF) analysis,
including line data, load data, generator data, bus-generator mappings, and generator status.

# Dependencies
- `PowerSystems`: For power system data handling.
- `DataFrames`, `CSV`: For data input and processing.
- `LinearAlgebra`: For matrix operations.
- `Printf`: For logging.

# Exported Functions
- `get_PS_data`: Main function to extract power system data from PTDF and ROP files.
"""
module PS_Data_DCOPF

using PowerSystems
using DataFrames
using LinearAlgebra
using CSV
using Printf

export get_PS_data

"""
    get_PS_data(PTDF_csvfile::String, rop_file::String, PS::System) -> Tuple{Matrix{Float64}, Matrix{Float32}, Matrix{Float32}, Matrix{Int64}, Vector{Int64}}

Extract power system data for DCOPF analysis, including line data, load data, generator data,
bus-generator mappings, and generator status.

# Arguments
- `PTDF_csvfile::String`: Path to the CSV file containing the PTDF matrix.
- `rop_file::String`: Path to the Resource Optimization Problem (ROP) file.
- `PS::System`: PowerSystems system object.

# Returns
- `_line_data::Matrix{Float64}`: Matrix containing line data [Line#, Max_Limit, Min_Limit, PTDFs].
- `Acc_Load::Matrix{Float32}`: Matrix containing load data [Bus#, Total Load].
- `_Gen_data::Matrix{Float32}`: Matrix containing generator data [Gen#, Bus#, P_max, P_min, cost_intervals, interval_power, interval_cost...].
- `_Bus_Gen::Matrix{Int64}`: Matrix mapping buses to generators (1 for active, 0 for inactive).
- `_GEN_STAT::Vector{Int64}`: Vector of generator status (1 for active, 0 for inactive).

# Notes
- The PTDF matrix is loaded from the CSV file and appended to the line data.
- Generator cost data is extracted from the ROP file using the `rop_parser` function.

# Examples
```julia
julia> using PS_Data_DCOPF
julia> line_data, load_data, gen_data, bus_gen, gen_stat = get_PS_data("ptdf.csv", "case.rop", system)
```
"""
function get_PS_data(PTDF_csvfile, rop_file, PS)
    _no_of_buses = length(get_bus_numbers(PS))
    _no_of_lines = length(get_components(ACBranch, PS))
    _PTDF_file = CSV.File(joinpath(PTDF_csvfile); header=true)
    _PTDF = DataFrame(_PTDF_file)
    _PTDF = Matrix(_PTDF)
    _line_data = zeros(Float64, _no_of_lines, _no_of_buses + 3)
    in = 0
    for q in get_components(ACBranch, PS)
        in += 1
        _line_data[in, 1] = in
        _line_data[in, 2] = 100 * get_rating(q)
        _line_data[in, 3] = -100 .* get_rating(q)
    end
    _line_data[:, 4:_no_of_buses+3] = _PTDF

    _P_Load = []
    _P_Load_Bus = []
    for q in get_components(StandardLoad, PS)
        push!(_P_Load, 100 * get_constant_active_power(q))
        push!(_P_Load_Bus, get_number(get_bus(q)))
    end
    Acc_Load = zeros(Float32, _no_of_buses, 2)
    M = (size(_P_Load))[1]
    for i in 1:M
        for j in 1:_no_of_buses
            if _P_Load_Bus[i] == j
                Acc_Load[j, 2] = Acc_Load[j,] + _P_Load[i]
            end
        end
    end

    Gen_dispatch, APower_dispatch, Cost_Table = rop_parser(rop_file, PS)
    _no_of_gens = (size(Gen_dispatch))[1]
    _Gen_data = zeros(Float32, _no_of_gens, 17)
    _Gen_data[:, 1] = Gen_dispatch[:, 4]
    _Gen_data[:, 2] = Gen_dispatch[:, 1]
    _Gen_data[:, 3] = APower_dispatch[:, 2]
    _Gen_data[:, 4] = APower_dispatch[:, 3]
    _Gen_data[:, 5:17] = Cost_Table[:, 2:14]

    _f = 0
    gen_bus_num = []
    gen_num = []
    gen_status = []
    for _jj in get_components(ThermalStandard, PS)
        _f = length(string(get_number(get_bus(_jj))))
        push!(gen_bus_num, get_number(get_bus(_jj)))
        push!(gen_num, parse(Int64, get_name(_jj)[11+_f+1]))
        if get_status(_jj) == true
            status = 1
        else
            status = -1
        end
        push!(gen_status, status)
    end
    final = []
    final = hcat(gen_bus_num, gen_num)
    final = hcat(final, gen_status)
    max_gen = findmax(final[:, 2])
    bus_gen = zeros(Int64, _no_of_buses, max_gen[1])

    for er in 1:(size(final))[1]
        bus_gen[final[er, 1], final[er, 2]] = final[er, 3]
    end
    gen = 0
    _Bus_Gen = zeros(Int64, _no_of_buses, _no_of_gens)
    for ia in 1:_no_of_buses
        for _ja in 1:(size(bus_gen))[2]
            if (bus_gen[ia, _ja] == 1)
                gen += 1
                _Bus_Gen[ia, gen] = 1
            elseif (bus_gen[ia, _ja] == -1)
                gen += 1
                _Bus_Gen[ia, gen] = 0
            end
        end
    end
    _GEN_STAT = zeros(Int64, _no_of_gens,)
    for temp_gen_stat in 1:_no_of_gens
        _GEN_STAT[temp_gen_stat] = sum(_Bus_Gen[:, temp_gen_stat])
    end

    return _line_data, Acc_Load, _Gen_data, _Bus_Gen, _GEN_STAT
end

"""
    rop_parser(rop_file::String, PSYS::System) -> Tuple{Matrix{Float32}, Matrix{Float32}, Matrix{Float32}}

Parse a Resource Optimization Problem (ROP) file to extract generator dispatch, active power dispatch,
and piece-wise linear cost data.

# Arguments
- `rop_file::String`: Path to the ROP file.
- `PSYS::System`: PowerSystems system object.

# Returns
- `Gen_dispatch::Matrix{Float32}`: Generator dispatch data [Gen#, Bus#, Dispatch, Status].
- `APower_dispatch::Matrix{Float32}`: Active power dispatch data [Gen#, P, P_min, P_max, Ramp, Reserve, Cost].
- `Cost_Table::Matrix{Float32}`: Piece-wise linear cost data [Gen#, Segments, Power1, Cost1, ...].
"""
function rop_parser(rop_file, PSYS)
    file = CSV.File(rop_file; header=false)
    r_file = (size(file))[1]
    k = 0
    GN_NO = 0
    for GN in get_components(ThermalStandard, PSYS)
        GN_NO += 1
    end
    @info "Number of generators: $GN_NO"
    Gen_dispatch = zeros(Float32, GN_NO, 4)
    APower_dispatch = zeros(Float32, GN_NO, 7)
    Cost_Table = zeros(Float32, 2320, 14)
    for i in 1:r_file
        if file[i,][1] == "0 / END OF ADJUSTABLE BUS LOAD TABLES BEGIN GENERATOR DISPATCH DATA"
            j = 1
            for j in i+1:r_file
                if file[j,][1] == "0 / END OF GENERATOR DISPATCH DATA BEGIN ACTIVE POWER DISPATCH TABLES"
                    break
                end
                Gen_dispatch[j-i, 1] = parse(Int, file[j,][1])
                Gen_dispatch[j-i, 2] = parse(Int, file[j,][2])
                Gen_dispatch[j-i, 3] = file[j,][3]
                Gen_dispatch[j-i, 4] = file[j,][4]
            end
        end

        if file[i,][1] == "0 / END OF GENERATOR DISPATCH DATA BEGIN ACTIVE POWER DISPATCH TABLES"
            j = 1
            for j in i+1:r_file
                if file[j,][1] == "0 / END OF ACTIVE POWER DISPATCH TABLES BEGIN GENERATION RESERVE DATA"
                    break
                end
                APower_dispatch[j-i, 1] = parse(Int, file[j,][1])
                APower_dispatch[j-i, 2] = parse(Float32, file[j,][2])
                APower_dispatch[j-i, 3] = file[j,][3]
                APower_dispatch[j-i, 4] = file[j,][4]
                APower_dispatch[j-i, 5] = file[j,][5]
                APower_dispatch[j-i, 6] = file[j,][6]
                APower_dispatch[j-i, 7] = file[j,][7]
            end
        end

        if file[i,][1] == "0 / END OF ADJUSTABLE BRANCH REACTANCE DATA BEGIN PIECE-WISE LINEAR COST TABLES"
            j = i + 1
            global m = 1
            while j <= r_file
                if file[j,][1] == "0 / END OF PIECE-WISE LINEAR COST TABLES BEGIN PIECEWISE QUADRATIC COST TABLES"
                    break
                end
                Cost_Table[m, 1] = parse(Float32, file[j,][1])
                Cost_Table[m, 2] = file[j,][3]
                if coalesce(file[j,][3] == 6)
                    Cost_Table[m, 3] = parse(Float32, file[j+1,][1])
                    Cost_Table[m, 4] = parse(Float32, file[j+1,][2])
                    Cost_Table[m, 5] = parse(Float32, file[j+2,][1])
                    Cost_Table[m, 6] = parse(Float32, file[j+2,][2])
                    Cost_Table[m, 7] = parse(Float32, file[j+3,][1])
                    Cost_Table[m, 8] = parse(Float32, file[j+3,][2])
                    Cost_Table[m, 9] = parse(Float32, file[j+4,][1])
                    Cost_Table[m, 10] = parse(Float32, file[j+4,][2])
                    Cost_Table[m, 11] = parse(Float32, file[j+5,][1])
                    Cost_Table[m, 12] = parse(Float32, file[j+5,][2])
                    Cost_Table[m, 13] = parse(Float32, file[j+6,][1])
                    Cost_Table[m, 14] = parse(Float32, file[j+6,][2])
                    j = j + 7
                elseif coalesce(file[j,][3] == 2)
                    Cost_Table[m, 3] = parse(Float32, file[j+1,][1])
                    Cost_Table[m, 4] = parse(Float32, file[j+1,][2])
                    Cost_Table[m, 5] = parse(Float32, file[j+2,][1])
                    Cost_Table[m, 6] = parse(Float32, file[j+2,][2])
                    j = j + 3
                end
                m += 1
            end
        end
    end
    Cost_Table = Cost_Table[1:m-1, :]

    return Gen_dispatch, APower_dispatch, Cost_Table
end

end # module PS_Data_DCOPF