#!/bin/bash
# julia --project snr_measurement.jl --L 56 --type 1 --distance 7 --Delta 0.1&
# julia --project snr_measurement.jl --L 56 --type 1 --distance 7 --Delta 0.2&
# wait
# julia --project snr_measurement.jl --L 56 --type 1 --distance 7 --Delta 0.3&
# julia --project snr_measurement.jl --L 56 --type 1 --distance 7 --Delta 0.4&
# wait
julia --project snr_measurement.jl --L 56 --type 1 --distance 7 --Delta 0.5&
julia --project snr_measurement.jl --L 56 --type 1 --distance 7 --Delta 0.6&
wait
julia --project snr_measurement.jl --L 56 --type 1 --distance 7 --Delta 0.7&
julia --project snr_measurement.jl --L 56 --type 1 --distance 7 --Delta 0.8&
wait
julia --project snr_measurement.jl --L 56 --type 1 --distance 7 --Delta 0.9&
julia --project snr_measurement.jl --L 56 --type 1 --distance 7 --Delta 1.0&
wait
julia --project snr_measurement.jl --L 56 --type 1 --distance 7 --Delta 1.2&
julia --project snr_measurement.jl --L 56 --type 1 --distance 7 --Delta 1.4&
wait
julia --project snr_measurement.jl --L 56 --type 1 --distance 7 --Delta 1.6&
julia --project snr_measurement.jl --L 56 --type 1 --distance 7 --Delta 1.8&
wait
julia --project snr_measurement.jl --L 56 --type 1 --distance 7 --Delta 2.0&