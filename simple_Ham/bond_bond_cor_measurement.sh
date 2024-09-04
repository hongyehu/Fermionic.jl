#!/bin/bash

# Loop through values of super_x, super_y, and index
for super_x in {1..6}; do
  for super_y in {1..6}; do
    for index in {1..8}; do
      # Exclude the case when super_x == super_y == index == 1
      if [ $super_x -eq 1 ] && [ $super_y -eq 1 ] && [ $index -eq 1 ]; then
        continue
      fi
      # Run the Julia script in series (one after the other)
      julia --project bond_bond_cor_measurement.jl --L 24 --type 1 --super_x $super_x --super_y $super_y --index $index --samples 500
    done
  done
done
