"""
DCOPF.jl

This script implements a Direct Current Optimal Power Flow (DCOPF) model using the PTDF matrix
and power system data. It minimizes the generation cost subject to power balance, line flow,
and generation constraints.

# Dependencies
- `PTDFMatrix.jl`: Module for computing the Power Transfer Distribution Factor (PTDF) matrix.
- `PS_Data_DCOPF.jl`: Module for extracting power system data.
- `JuMP`, `PiecewiseLinearOpt`, `HiGHS`: Optimization packages for modeling and solving.
- `CSV`, `Tables`: For data input/output.

# Inputs
- `raw_file_dir::String`: Path to the raw power system data file.
- `PTDF_csv_dir::String`: Path to the output CSV file for the PTDF matrix.
- `rop_file_dir::String`: Path to the Resource Optimization Problem (ROP) file.

# Outputs
- A JuMP model is solved to determine optimal generator dispatch (`Pg`).
- The PTDF matrix and system data are retrieved and used in the optimization.

# Examples
```julia
julia> include("DCOPF.jl")
# Runs the DCOPF model and optimizes generator dispatch
```
"""

include("PTDFMatrix.jl")
include("PS_Data_DCOPF.jl")

using .PTDFMatrix
using .PS_Data_DCOPF
using JuMP
using PiecewiseLinearOpt
using HiGHS
using CSV
using Tables

# Input file paths
raw_file_dir = "grid/case.raw"
PTDF_csv_dir = "grid/case.CSV"
rop_file_dir = "grid/case.rop"

# Retrieve PTDF matrix and power system
PS_PTDF, P_system = get_PTDF_matrix(raw_file_dir, PTDF_csv_dir)

# Extract power system data
Line_Data, Load_Data, Gen_Data, Bus_Gen, Gen_Stat = get_PS_data(PTDF_csv_dir, rop_file_dir, P_system)
@info "Data dimensions:" Line=size(Line_Data) Load=size(Load_Data) Generators=size(Gen_Data) Bus_Gen=size(Bus_Gen)

# Prepare optimization parameters
PL_max = Line_Data[:, 2]
PL_min = Line_Data[:, 3]
Pg_max = Gen_Data[:, 3] .* Gen_Stat
Pg_min = Gen_Data[:, 4] .* Gen_Stat
PWL_Cost = Gen_Data[:, 5:17] .* Gen_Stat
NG = size(Pg_max, 1)
NB = size(Bus_Gen, 1)

# Define and solve the DCOPF model
model = Model()

@variable(model, Pg[1:NG])
@constraint(model, c1, 1.5 * PL_min .<= PS_PTDF * (Bus_Gen * Pg - Load_Data[:, 2]) .<= 1.5 * PL_max)
@constraint(model, c2, Pg_min .<= Pg .<= Pg_max)
@constraint(model, c3, sum(Pg[G] * Gen_Stat[G] for G in 1:NG) == sum(Load_Data[B, 2] for B in 1:NB))

# Define piecewise linear cost function for each generator
obj_func = []
for gen_i in 1:NG
    x = Float64[]
    y = Float64[]
    if Pg_max[gen_i] == 0
        push!(x, 0)
        push!(y, 0)
    else
        for ind in 1:2:Int(PWL_Cost[gen_i, 1]) * 2
            push!(x, PWL_Cost[gen_i, ind + 1])
            push!(y, PWL_Cost[gen_i, ind + 2])
        end
    end
    push!(obj_func, piecewiselinear(model, Pg[gen_i], x, y))
end

@objective(model, Min, sum(obj_func[i] for i in 1:NG))

# Set optimizer and solve
set_optimizer(model, HiGHS.Optimizer)
optimize!(model)

