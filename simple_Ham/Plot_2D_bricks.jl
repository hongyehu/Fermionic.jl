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

plt = pyimport("matplotlib.pyplot")
cm = pyimport("matplotlib.cm")
patches = pyimport("matplotlib.patches")
# mcolors = pyimport("matplotlib.colors")

"""
A function that generate xy location, width, height of the FancyBboxPatch
    based on super lattice 2D location and index within super lattice
"""
function generate_patch_info(super_x,super_y,index)
    tmp = 0.5
    if index == 1
        x = -tmp/2
        y = -tmp/2
        width = tmp
        height = 1.0+tmp
    elseif index == 2
        x = -tmp/2+1.0
        y = -tmp/2
        width = tmp
        height = 1.0+tmp
    elseif index == 3
        x = 2.0-tmp/2
        y = -tmp/2
        width = 1.0+tmp
        height = tmp
    elseif index == 4
        x = 2.0-tmp/2
        y = -tmp/2+1.0
        width = 1.0+tmp
        height = tmp
    elseif index == 5
        x = -tmp/2
        y = 2.0-tmp/2
        width = 1.0+tmp
        height = tmp
    elseif index == 6
        x = -tmp/2
        y = 3.0-tmp/2
        width = 1.0+tmp
        height = tmp
    elseif index == 7
        x = 2.0-tmp/2
        y = 2.0-tmp/2
        width = tmp
        height = 1.0+tmp
    elseif index == 8
        x = 3.0-tmp/2
        y = 2.0-tmp/2
        width = tmp
        height = 1.0+tmp
    end
    return x+4.0*(super_x-1),y+4.0*(super_y-1),width,height
end
rectangles = []
for super_x in 1:6
    for super_y in 1:6
        for index in 1:8
            x,y,width,height = generate_patch_info(super_x,super_y,index)
            push!(rectangles,Dict("xy"=>(x,y),"width"=>width,"height"=>height))
        end
    end
end
L=24
type = :d
G_file_name = "/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/simple_direct_measurement_data/G/L$(L)_type1_Delta0.3_mu0.5.jld2"
@load G_file_name G

bond_bond_correlation(G, 1,4,6)

scalars = [] # Example scalar values
for super_x in 1:6
    for super_y in 1:6
        for index in 1:8
            println("super x:$(super_x) super y:$(super_y) index:$(index)")
            push!(scalars,bond_bond_correlation(G, super_x,super_y,index))
        end
    end
end

scalars
rescaled_scalars = real.(scalars)
begin
    plt.plot(rescaled_scalars)
    plt.show()
end
max_value = maximum(abs.(rescaled_scalars))
rescaled_scalars = (rescaled_scalars.+max_value)./(2*max_value)
### Test on plot the batches
begin
    # Create a figure and axis
    fig, ax = plt.subplots(figsize = (3,3))
    boxstyle = "round,pad=0.02,rounding_size=0.1" 
    colormap = cm.bwr
    colors = colormap(rescaled_scalars)  # Map scalar values to colors

    # Add rectangles to the plot
    for i in 1:size(rectangles)[1]
        rect = rectangles[i]
        ax.add_patch(
            patches.FancyBboxPatch(
                rect["xy"], rect["width"], rect["height"],
                edgecolor="black",
                boxstyle=boxstyle,
                facecolor=colors[i-1],
                linewidth=0.5,
                alpha=0.5  # Transparency
            )
        )
    end
    # Set limits, labels, etc.
    ax.set_xlim(-1, 4*6)
    ax.set_ylim(-1, 4*6)
    ax.set_aspect("equal", "box")
    ax.axis("off")

    # Show plot
    plt.show()
end