# BcubeVTK.jl
Implementation of [`Bcube`](https://github.com/bcube-project/Bcube.jl) IO interface for VTK format. Checkout the relative `Bcube` [documentation](https://bcube-project.github.io/Bcube.jl/stable/api/io/io_interface/) for more infos.

For now, only the `write_file` interface is implemented.

## Basic usage
```julia
using Bcube
using BcubeVTK

mesh = rectangle_mesh(10, 20)
U = TrialFESpace(FunctionSpace(:Lagrange, 1), mesh)
u = FEFunction(U)
projection_l2!(u, PhysicalFunction(x -> sum(x)), CellDomain(mesh))

write_file("output.pvd", mesh, Dict("u" => u, "grad_u" => ∇(u)))
```

## Limitations

* The `write_file` has not been tested for hexahedral elements of order > 1
