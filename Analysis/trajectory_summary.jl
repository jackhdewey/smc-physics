# Data analysis
# Concatenates trajectory plot images into a grid for side-by-side comparison

using Plots
using FileIO
using Images

include("../Utilities/fileio.jl")

dir1 = string("Analysis/Plots/Modelv5/Sphere/ObsVar1/Test/")
dir2 = string("Analysis/Plots/Modelv5/Sphere/ObsVar075/Test/")
dir3 = string("Analysis/Plots/Modelv5/Sphere/PosVar05/Test/")
dir4 = string("Analysis/Plots/Modelv5/Sphere/ObsVar025/Test/")

files1 = readdir(dir1)
sort!(files1, lt=png_particle_order)
files2 = readdir(dir2)
sort!(files2, lt=png_particle_order)
files3 = readdir(dir3)
sort!(files3, lt=png_particle_order)
files4 = readdir(dir4)
sort!(files4, lt=png_particle_order)

for i in eachindex(files1)
    plt1 = load(string(dir1, files1[i]))
    plt2 = load(string(dir2, files2[i]))
    plt3 = load(string(dir3, files3[i]))
    plt4 = load(string(dir4, files4[i]))
    display(mosaicview(plt1, plt3, plt2, plt4; nrow=2))
end