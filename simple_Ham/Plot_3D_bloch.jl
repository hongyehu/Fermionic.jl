using LinearAlgebra
using SparseArrays
using Fermionic
using ProgressMeter
using Statistics
using ArgParse
using JSON3
using JLD2
using PythonCall
include("./analytical_2d_bond.jl")
import numpy as np
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

# # Create a sphere
# phi, theta = np.mgrid[0.0:2.0 * np.pi:100j, 0.0:np.pi:50j]
# x = np.sin(theta) * np.cos(phi)
# y = np.sin(theta) * np.sin(phi)
# z = np.cos(theta)

# # Create a 3D plot
# fig = plt.figure(figsize=(8, 6))
# ax = fig.add_subplot(111, projection='3d')

# # Plot the surface of the sphere
# ax.plot_surface(x, y, z, color='w', alpha=0.3, rstride=5, cstride=5, linewidth=0.5, edgecolor='black')

# # Now, let's add some colored dots on the sphere at specific angles
# n_dots = 20  # Number of dots to place on the sphere

# # These angles determine the placement of the dots
# theta = np.linspace(0, np.pi, n_dots)
# phi = np.linspace(0, 2 * np.pi, n_dots)

# # Convert spherical coordinates to Cartesian coordinates for the dots
# x_dots = np.sin(theta) * np.cos(phi)
# y_dots = np.sin(theta) * np.sin(phi)
# z_dots = np.cos(theta)

# # Plot the dots
# ax.scatter(x_dots, y_dots, z_dots, color='blue')

# # Setting the labels
# ax.set_xlabel('X axis')
# ax.set_ylabel('Y axis')
# ax.set_zlabel('Z axis')

# Show the plot