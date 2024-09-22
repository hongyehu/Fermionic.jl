#!/bin/bash
julia --project hypothesis_measurement.jl --L 40 --type 0 --distance 5 --Delta 0.1 --samples 10000000 --left 1 --right 1
julia --project hypothesis_measurement.jl --L 40 --type 0 --distance 5 --Delta 0.1 --samples 10000000 --left 1 --right 2
julia --project hypothesis_measurement.jl --L 40 --type 0 --distance 5 --Delta 0.1 --samples 10000000 --left 1 --right 3
julia --project hypothesis_measurement.jl --L 40 --type 0 --distance 5 --Delta 0.1 --samples 10000000 --left 1 --right 4
julia --project hypothesis_measurement.jl --L 40 --type 0 --distance 5 --Delta 0.1 --samples 10000000 --left 2 --right 1
julia --project hypothesis_measurement.jl --L 40 --type 0 --distance 5 --Delta 0.1 --samples 10000000 --left 2 --right 2
julia --project hypothesis_measurement.jl --L 40 --type 0 --distance 5 --Delta 0.1 --samples 10000000 --left 2 --right 3
julia --project hypothesis_measurement.jl --L 40 --type 0 --distance 5 --Delta 0.1 --samples 10000000 --left 2 --right 4
julia --project hypothesis_measurement.jl --L 40 --type 0 --distance 5 --Delta 0.1 --samples 10000000 --left 3 --right 1
julia --project hypothesis_measurement.jl --L 40 --type 0 --distance 5 --Delta 0.1 --samples 10000000 --left 3 --right 2
julia --project hypothesis_measurement.jl --L 40 --type 0 --distance 5 --Delta 0.1 --samples 10000000 --left 3 --right 3
julia --project hypothesis_measurement.jl --L 40 --type 0 --distance 5 --Delta 0.1 --samples 10000000 --left 3 --right 4
julia --project hypothesis_measurement.jl --L 40 --type 0 --distance 5 --Delta 0.1 --samples 10000000 --left 4 --right 1
julia --project hypothesis_measurement.jl --L 40 --type 0 --distance 5 --Delta 0.1 --samples 10000000 --left 4 --right 2
julia --project hypothesis_measurement.jl --L 40 --type 0 --distance 5 --Delta 0.1 --samples 10000000 --left 4 --right 3
julia --project hypothesis_measurement.jl --L 40 --type 0 --distance 5 --Delta 0.1 --samples 10000000 --left 4 --right 4
