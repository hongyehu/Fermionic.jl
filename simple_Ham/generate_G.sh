#!/bin/bash
julia --project generate_G.jl --L 40 --type 1 --Delta 0.2&
julia --project generate_G.jl --L 40 --type 0 --Delta 0.2&