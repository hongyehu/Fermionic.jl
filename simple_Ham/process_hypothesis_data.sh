#!/bin/bash
julia --project hypothesis_measurement_process.jl --L 40 --type 1 --distance 5 --Delta 0.2 --measurements 50
julia --project hypothesis_measurement_process.jl --L 40 --type 1 --distance 5 --Delta 0.2 --measurements 100
julia --project hypothesis_measurement_process.jl --L 40 --type 1 --distance 5 --Delta 0.2 --measurements 150
julia --project hypothesis_measurement_process.jl --L 40 --type 1 --distance 5 --Delta 0.2 --measurements 200
julia --project hypothesis_measurement_process.jl --L 40 --type 1 --distance 5 --Delta 0.2 --measurements 300
julia --project hypothesis_measurement_process.jl --L 40 --type 1 --distance 5 --Delta 0.2 --measurements 400
julia --project hypothesis_measurement_process.jl --L 40 --type 1 --distance 5 --Delta 0.2 --measurements 500
julia --project hypothesis_measurement_process.jl --L 40 --type 1 --distance 5 --Delta 0.2 --measurements 600
julia --project hypothesis_measurement_process.jl --L 40 --type 1 --distance 5 --Delta 0.2 --measurements 1000


