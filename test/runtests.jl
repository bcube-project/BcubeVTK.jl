using Test
using BcubeVTK
using Bcube
using StaticArrays
using BcubeVTK: write_vtk_lagrange, read_file, write_file
using SHA
using ReferenceTests
using LinearAlgebra

include("xml_parser.jl")
include("utils.jl")

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
tempdir = mktempdir(; cleanup = true)
@show tempdir

@testset "BcubeVTK.jl" begin
    custom_include("./test_vtk.jl")
end
