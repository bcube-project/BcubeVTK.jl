using Test
using BcubeVTK
using DelimitedFiles
using Bcube
using StaticArrays
import BcubeVTK: write_vtk, write_vtk_lagrange
using WriteVTK # TODO : get rid of it
using SHA

"""
Custom way to "include" a file to print infos.
"""
function custom_include(path)
    filename = split(path, "/")[end]
    print("Running test file " * filename * "...")
    include(path)
    println("done.")
end

# This dir will be removed at the end of the tests
tempdir = mktempdir()

# Reading sha1 checksums
f = readdlm(joinpath(@__DIR__, "checksums.sha1"), String)
fname2sum = Dict(r[2] => r[1] for r in eachrow(f))

@testset "BcubeVTK.jl" begin
    custom_include("./test_vtk.jl")
end
