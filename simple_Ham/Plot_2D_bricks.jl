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
"""
Real Data
"""
scalars = [] # Example scalar values
for super_x in 1:6
    for super_y in 1:6
        for ori_index in 1:8
            shift_super_x = mod(super_x-3,6)+1
            shift_super_y = mod(super_y-3,6)+1
            if ori_index == 1
                index = 8
            elseif ori_index == 2
                index = 7
            elseif ori_index == 3
                index = 6
            elseif ori_index == 4
                index = 5
            elseif ori_index == 5
                index = 4
            elseif ori_index == 6
                index = 3
            elseif ori_index == 7
                index = 2
            elseif ori_index == 8
                index = 1
            end
            if shift_super_x == 1 && shift_super_y == 1 && index == 1
                push!(scalars,0.0)
            else
                println("super x:$(super_x) super y:$(super_y) index:$(index)")
                file_name = "/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/simple_direct_measurement_data/bond_bond_data/L24_type1_superx$(shift_super_x)_supery$(shift_super_y)_index$(index)_Delta0.3_mu0.5_samples500.json"
                open(file_name, "r") do io
                    data = JSON3.read(io)
                    tmp = data["mean"]
                    push!(scalars,tmp)
                end
            end
        end
    end
end


scalars

"""
Synthetic Data 
"""
# L=24
# type = :d
# G_file_name = "/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/simple_direct_measurement_data/G/L$(L)_type1_Delta0.3_mu0.5.jld2"
# @load G_file_name G

# bond_bond_correlation(G, 1,4,6)

# scalars = [] # Example scalar values
# for super_x in 1:6
#     for super_y in 1:6
#         for index in 1:8
#             println("super x:$(super_x) super y:$(super_y) index:$(index)")
#             push!(scalars,bond_bond_correlation(G, super_x,super_y,index))
#         end
#     end
# end
# scalars
# rescaled_scalars = real.(scalars)
rectangles
"""
Plot 2D bond-bond-bricks
"""
rescaled_scalars = real.(scalars)
for i in 1:size(scalars)[1]
    if real(scalars[i])==0.0
        println(i)
    end
end
rectangles[120]
begin
    plt.plot(rescaled_scalars,"o")
    plt.show()
end
max_value = maximum(abs.(rescaled_scalars))
rescaled_scalars = (rescaled_scalars.+max_value)./(2*max_value)
begin
    plt.plot(rescaled_scalars,"o")
    plt.show()
end

### Test on plot the batches
begin
    # Create a figure and axis
    fig, ax = plt.subplots(figsize = (3,3))
    boxstyle = "round,pad=0.02,rounding_size=0.1" 
    colormap = cm.seismic
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
    ax.add_patch(
        patches.FancyBboxPatch(
            (10.75,9.75), 0.5, 1.5,
            edgecolor="green",
            boxstyle=boxstyle,
            facecolor="whitesmoke",  # No fill
            hatch="//////",  # Diagonal stripes (you can use different patterns)
            linewidth=1.0,     # Optional: adjust the line width as needed
            alpha=1.0         # Optional: full opacity
        )
    )
    # Set limits, labels, etc.
    ax.set_xlim(-1, 4*6)
    ax.set_ylim(-1, 4*6)
    ax.set_aspect("equal", "box")
    # Create the color bar using the colormap and normalization
    sm = plt.cm.ScalarMappable(cmap=colormap)
    # sm.set_array([])  # Required for colorbar to work
    # cbar = fig.colorbar(sm, ax=ax)
    # cbar.set_label('Scalar Value')  # Optional: label for the color bar
    ax.axis("off")
    fig.savefig("./bond_bond_correlation_2.svg",bbox_inches = "tight")
    # Show plot
    plt.show()
end

Normalize = pyimport("matplotlib.colors").Normalize
np = pyimport("numpy")
begin
    # Example data array
    data_values = np.linspace(-max_value,max_value,100)  # Your data array
    # Normalize the data based on its min and max
    norm = Normalize(vmin=np.min(data_values), vmax=np.max(data_values))
    # Create a colormap (you can use any colormap)
    colormap = cm.seismic
    # Create the ScalarMappable object to map the data to colors
    sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)
    # Create the color bar
    fig, ax = plt.subplots(figsize=(2, 5))  # Adjust the size as needed
    cbar = fig.colorbar(sm, ax=ax, orientation="vertical")  # Vertical color bar
    # cbar.set_label("Data Value")  # Add a label for the color bar
    # Show the plot with the color bar
    plt.axis("off")
    plt.savefig("./colorbar_2.svg",bbox_inches = "tight")
    plt.show()
end