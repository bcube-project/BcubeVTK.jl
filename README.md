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

## Keyword arguments

Besides the keyword arguments defined by the Bcube
[`write_file`](https://bcube-project.github.io/Bcube.jl/stable/api/io/io_interface/)
interface, `BcubeVTK` exposes a couple of VTK-specific ones. Any other keyword
argument is transparently forwarded to
[`WriteVTK.jl`](https://github.com/JuliaVTK/WriteVTK.jl) (more precisely to
`vtk_grid`, and to `paraview_collection` / `vtk_save` which are called
internally).

* `collection_append::Bool = false`: when `true`, the current snapshot is
  appended to an existing `.pvd` collection (i.e. `append = true` is passed to
  `paraview_collection`) instead of (over)writing a brand new one. This is handy
  when several time steps of a time-dependent simulation are written from
  independent `write_file` calls and you want them all gathered in the same
  `.pvd` file. Note that the keyword is named `collection_append` (and not
  `append`) precisely so that an `append` kwarg can still be forwarded to
  WriteVTK's `vtk_grid`.
* `pad_vectors_to_3d::Bool = true`: VTK only recognises a data array as a
  *vector* field when it has exactly 3 components. On a 2D mesh, vector fields
  naturally have only 2 components, so when this option is enabled each 2D
  vector is padded with a zero third component, allowing ParaView to interpret
  it as a proper vector field (and to use e.g. glyph or streamline filters).
  Set it to `false` to keep the original number of components.

All remaining keyword arguments are passed through to `WriteVTK.jl`. For
example, pass `ascii = true` (together with `append = false`) to write an ASCII
file instead of a binary/compressed one.

## Limitations

* The `write_file` has not been tested for hexahedral elements of order > 1
