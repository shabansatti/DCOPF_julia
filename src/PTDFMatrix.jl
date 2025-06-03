"""
PTDFMatrix.jl

This module provides functionality to compute the Power Transfer Distribution Factor (PTDF) matrix
for a power system, which describes the sensitivity of power flows on lines to injections at buses.
The PTDF matrix is calculated using the PowerSystems package and written to a CSV file.

# Dependencies
- `PowerSystems`: For power system data handling.
- `DataFrames`, `CSV`, `Tables`: For data manipulation and output.
- `LinearAlgebra`, `SparseArrays`: For matrix operations.
- `Printf`: For logging.

# Exported Functions
- `get_PTDF_matrix`: Main function to compute the PTDF matrix and save it to a CSV file.
"""
module PTDFMatrix

using PowerSystems
using DataFrames
using LinearAlgebra
using SparseArrays
using CSV
using Printf

export get_PTDF_matrix

"""
    get_PTDF_matrix(file_dir::String, out_file::String)

Compute the Power Transfer Distribution Factor (PTDF) matrix for a power system and write it to a CSV file.

# Arguments
- `file_dir::String`: Path to the power system data file.
- `out_file::String`: Path to the output CSV file for the PTDF matrix.

# Returns
- `PTDF::Matrix{Float64}`: The PTDF matrix, including the reference bus.
- `power_system::System`: The PowerSystems system object.

# Notes
- The PTDF matrix is computed by excluding the reference bus during intermediate calculations
  and reinserting it as a column of zeros in the final matrix.
- The output CSV file contains the PTDF matrix with bus numbers as headers.

# Examples
```julia
julia> using PTDFMatrix
julia> PTDF, system = get_PTDF_matrix("system_data", "output.csv")
```
"""
function get_PTDF_matrix(file_dir::String, out_file::String)
    power_system = System(joinpath(file_dir); bus_name_formatter=x -> "$(x["name"])-$(x["index"])")
    ref_bus, n_bus, n_line = get_system_params(power_system)
    H_t, B_hat_t = calculate_H_B_matrices(power_system, ref_bus, n_bus, n_line)
    PTDF_no_ref = calculate_PTDF(B_hat_t, H_t, n_bus, n_line)
    PTDF = add_refbus_in_PTDF(PTDF_no_ref, ref_bus, n_bus, n_line)
    output_csv_file(power_system, PTDF, out_file)
    return PTDF, power_system
end

"""
    get_system_params(system::System) -> Tuple{Int, Int, Int}

Extract system parameters from a PowerSystems system object.

# Arguments
- `system::System`: The PowerSystems system object.

# Returns
- `ref_bus::Int`: Reference bus number.
- `n_bus::Int`: Total number of buses.
- `n_line::Int`: Total number of lines.
"""
function get_system_params(system::System)
    ref_bus = 0
    for b in get_components(ACBus, system)
        if string(get_bustype(b)) == "REF"
            ref_bus = get_number(b)
            break
        end
    end
    @info "Reference Bus: $ref_bus"

    n_bus = length(get_bus_numbers(system))
    n_line = length(get_components(ACBranch, system))
    return ref_bus, n_bus, n_line
end

"""
    calculate_H_B_matrices(system::System, ref_bus::Int, n_bus::Int, n_line::Int) -> Tuple{SparseMatrixCSC, SparseMatrixCSC}

Calculate the incidence (H_t) and weighted admittance (B_hat_t) matrices for PTDF computation.

# Arguments
- `system::System`: The PowerSystems system object.
- `ref_bus::Int`: Reference bus number.
- `n_bus::Int`: Number of buses.
- `n_line::Int`: Number of lines.

# Returns
- `H_t::SparseMatrixCSC`: Transpose of the weighted incidence matrix.
- `B_hat_t::SparseMatrixCSC`: Weighted admittance matrix.
"""
function calculate_H_B_matrices(system::System, ref_bus::Int, n_bus::Int, n_line::Int)
    lines = collect(1:n_line)
    from_bus = zeros(Int64, n_line)
    to_bus = zeros(Int64, n_line)
    x_inv = zeros(Float16, n_line)
    z = 0

    for q in get_components(ACBranch, system)
        z += 1
        from_bus[z] = get_number(get_from(get_arc(q)))
        to_bus[z] = get_number(get_to(get_arc(q)))
        x_inv[z] = round(-1 / get_x(q), digits=4)
    end

    X = zeros(Int64, n_line, 3)
    X = hcat(hcat(hcat(lines, from_bus), to_bus))
    Y = zeros(Float16, n_line, 2)
    for i in 1:n_line
        if X[i, 2] > X[i, 3]
            Y[i, 1] = -1
            Y[i, 2] = 1
        else
            Y[i, 1] = 1
            Y[i, 2] = -1
        end
    end
    L = vcat(X[:, 1], X[:, 1])
    M = vcat(X[:, 2], X[:, 3])
    N = vcat(Y[:, 1], Y[:, 2])

    A = sparse(L, M, N)
    A = A[:, 1:end .!= ref_bus]
    A_t = sparse(M, L, N)
    A_t = A_t[1:end .!= ref_bus, :]

    B = sparse(lines, lines, x_inv)
    H = B * A
    H_t = A_t * B
    B_hat = A_t * H
    B_hat_t = H_t * A
    return H_t, B_hat_t
end

"""
    calculate_PTDF(B_hat_t::SparseMatrixCSC, H_t::SparseMatrixCSC, n_bus::Int, n_line::Int) -> Matrix{Float64}

Compute the PTDF matrix excluding the reference bus.

# Arguments
- `B_hat_t::SparseMatrixCSC`: Weighted admittance matrix.
- `H_t::SparseMatrixCSC`: Transpose of the weighted incidence matrix.
- `n_bus::Int`: Number of buses.
- `n_line::Int`: Number of lines.

# Returns
- `PTDF::Matrix{Float64}`: PTDF matrix without the reference bus.
"""
function calculate_PTDF(B_hat_t::SparseMatrixCSC, H_t::SparseMatrixCSC, n_bus::Int, n_line::Int)
    F_fac = lu(B_hat_t)
    L = F_fac.L
    U = F_fac.U
    p = F_fac.p
    Rs = F_fac.Rs
    q = F_fac.q
    Rs = Rs[p, :]

    PTDF = zeros(Float64, n_line, n_bus - 1)
    for i in 1:n_line
        temp_3 = H_t[:, i]
        H_temp = Rs .* temp_3[p]
        temp_1 = L \ H_temp
        temp_2 = U \ temp_1
        PTDF[i, :] = temp_2[invperm(q)]
    end
    return PTDF
end

"""
    add_refbus_in_PTDF(PTDF_no_ref::Matrix{Float64}, ref_bus::Int, n_bus::Int, n_line::Int) -> Matrix{Float64}

Reinsert the reference bus column (zeros) into the PTDF matrix.

# Arguments
- `PTDF_no_ref::Matrix{Float64}`: PTDF matrix without the reference bus.
- `ref_bus::Int`: Reference bus number.
- `n_bus::Int`: Number of buses.
- `n_line::Int`: Number of lines.

# Returns
- `PTDF::Matrix{Float64}`: PTDF matrix with the reference bus column included.
"""
function add_refbus_in_PTDF(PTDF_no_ref::Matrix{Float64}, ref_bus::Int, n_bus::Int, n_line::Int)
    temp = zeros(Float16, n_line, 1)
    PTDF1 = PTDF_no_ref[:, 1:ref_bus-1]
    PTDF2 = PTDF_no_ref[:, ref_bus:n_bus-1]
    PTDF = hcat(PTDF1, temp)
    PTDF = hcat(PTDF, PTDF2)
    PTDF = round.(PTDF, digits=5)
    return PTDF
end

"""
    output_csv_file(system::System, PTDF::Matrix{Float64}, out_file::String)

Write the PTDF matrix to a CSV file with bus numbers as headers.

# Arguments
- `system::System`: The PowerSystems system object.
- `PTDF::Matrix{Float64}`: The PTDF matrix.
- `out_file::String`: Path to the output CSV file.
"""
function output_csv_file(system::System, PTDF::Matrix{Float64}, out_file::String)
    buses = get_bus_numbers(system)
    CSV.write(joinpath(out_file), Tables.table(PTDF); header=string.(buses), writeheader=true)
end

end # module PTDFMatrix