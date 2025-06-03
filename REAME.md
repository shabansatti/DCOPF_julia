# DCOPF

A Julia package for Direct Current Optimal Power Flow (DCOPF) analysis, computing optimal generator dispatch using Power Transfer Distribution Factor (PTDF) matrices and power system data.

---

## Overview

The `DCOPF` package implements a DCOPF model to minimize generation costs subject to power balance, line flow, and generation constraints. It uses the PTDF matrix to model line flows and integrates with the [`PowerSystems.jl`](https://github.com/NREL-Sienna/PowerSystems.jl) package for power system data handling.

The package consists of three main modules:

- **`PTDFMatrix.jl`**: Computes the PTDF matrix, describing the sensitivity of line flows to bus injections.
- **`PS_Data_DCOPF.jl`**: Extracts line, load, generator, and bus-generator mapping data from PTDF and Resource Optimization Problem (ROP) files.
- **`DCOPF.jl`**: Defines and solves the DCOPF optimization problem using the JuMP framework with the HiGHS optimizer.

---

## Installation

1. **Clone the repository**:
   ```bash
   git clone https://github.com/shabansatti/DCOPF.git
   cd DCOPF
   ```

2. **Install Julia**:  
   Download and install Julia (version 1.6 or higher) from the official website:  
   [https://julialang.org/downloads/](https://julialang.org/downloads/)

3. **Activate and install dependencies**:  
   In the Julia REPL, navigate to the project directory and run:
   ```julia
   using Pkg
   Pkg.activate(".")
   Pkg.instantiate()
   ```

---

## Usage

### Prepare Input Files

- A raw power system data file (e.g., `case.raw`)
- A Resource Optimization Problem (ROP) file (e.g., `case.rop`)

Place these files in the project directory or update the file paths in `src/DCOPF.jl`.

### Run the DCOPF Model

In the Julia REPL:
```julia
include("src/DCOPF.jl")
```

This script will:
- Compute the PTDF matrix and save it to a specified CSV file
- Extract power system data (line, load, generator, and bus-generator mappings)
- Define and solve the DCOPF model to determine the optimal generator dispatch

---

## Example

```julia
julia> include("src/DCOPF.jl")
# Output: Logs data dimensions and solves the DCOPF model
```

---
