"""
    write_vtk(basename::String, it::Int,time::Real, mesh::AbstractMesh{topoDim,spaceDim}, vars::Dict{String,Tuple{V,L}}; append=false) where{topoDim,spaceDim,V,L<:WriteVTK.AbstractFieldData}

Write a set of variables on the mesh nodes or cell centers to a VTK file.

# Example
```julia
mesh = basic_mesh()
u = rand(ncells(mesh))
v = rand(nnodes(mesh))
dict_vars = Dict( "u" => (u, VTKCellData()),  "v" => (v, VTKPointData()) )
write_vtk("output", 0, 0.0, mesh, dict_vars)
```
"""
function write_vtk(
    basename::String,
    it::Int,
    time::Real,
    mesh::AbstractMesh{topoDim, spaceDim},
    vars::Dict{String, Tuple{V, L}};
    append = false,
) where {topoDim, spaceDim, V, L <: WriteVTK.AbstractFieldData}
    pvd = paraview_collection(basename; append = append)

    # Create coordinates arrays
    vtknodes = reshape(
        [get_coords(n)[idim] for n in get_nodes(mesh) for idim in 1:spaceDim],
        spaceDim,
        nnodes(mesh),
    )

    # Connectivity
    c2n = connectivities_indices(mesh, :c2n)
    # Create cell array
    vtkcells =
        [MeshCell(vtk_entity(cells(mesh)[icell]), c2n[icell]) for icell in 1:ncells(mesh)]

    # Define mesh for vtk
    new_name = _build_fname_with_iterations(basename, it)
    vtkfile = vtk_grid(new_name, vtknodes, vtkcells)

    for (varname, (value, loc)) in vars
        vtkfile[varname, loc] = value
    end

    pvd[float(time)] = vtkfile
    vtk_save(pvd) # also triggers `vtk_save(vtkfile)`
end

function write_vtk_bnd_discontinuous(
    basename::String,
    it::Int,
    time::Real,
    domain::BoundaryFaceDomain,
    vars::Dict{String, Tuple{V, L}},
    degree::Int;
    append = false,
) where {V, L <: WriteVTK.AbstractFieldData}
    pvd = paraview_collection(basename; append = append)

    mesh = get_mesh(domain)
    sdim = spacedim(mesh)

    # Connectivities
    c2n = connectivities_indices(mesh, :c2n)
    f2n = connectivities_indices(mesh, :f2n)
    f2c = connectivities_indices(mesh, :f2c)

    # Cell and face types
    celltypes = cells(mesh)

    bndfaces = Bcube.get_cache(domain)

    fs = FunctionSpace(:Lagrange, max(1, degree)) # here, we implicitly impose that the mesh is composed of Lagrange elements only

    a = map(bndfaces) do iface
        icell = f2c[iface][1]
        sideᵢ = Bcube.cell_side(celltypes[icell], c2n[icell], f2n[iface])
        localfacedofs =
            Bcube.idof_by_face_with_bounds(fs, shape(celltypes[icell]))[sideᵢ]
        ξ = get_coords(fs, shape(celltypes[icell]))[localfacedofs]
        xdofs = map(
            _ξ -> Bcube.mapping(celltypes[icell], get_nodes(mesh, c2n[icell]), _ξ),
            ξ,
        )
        ftype =
            Bcube.entity(Bcube.face_shapes(shape(celltypes[icell]), sideᵢ), Val(degree))
        ftype, Bcube.rawcat(xdofs)
    end
    ftypes = getindex.(a, 1)
    vtknodes = reshape(Bcube.rawcat(getindex.(a, 2)), sdim, :)

    # Create elements array
    vtkcells = MeshCell[]
    count = 0
    for ftype in ftypes
        _nnode = get_ndofs(fs, shape(ftype))
        push!(vtkcells, MeshCell(vtk_entity(ftype), collect((count + 1):(count + _nnode))))
        count += _nnode
    end

    # Define mesh for vtk
    new_name = _build_fname_with_iterations(basename, it)
    vtkfile = vtk_grid(new_name, vtknodes, vtkcells)

    for (varname, (value, loc)) in vars
        vtkfile[varname, loc] = value
    end

    pvd[float(time)] = vtkfile
    vtk_save(pvd)
end

"""
    write_vtk(basename::String, mesh::AbstractMesh{topoDim,spaceDim}) where{topoDim,spaceDim}

Write the mesh to a VTK file.

# Example
```julia
write_vtk("output", basic_mesh())
```
"""
function write_vtk(
    basename::String,
    mesh::AbstractMesh{topoDim, spaceDim},
) where {topoDim, spaceDim}
    dict_vars = Dict{String, Tuple{Any, WriteVTK.AbstractFieldData}}()
    write_vtk(basename, 1, 0.0, mesh, dict_vars)
end

"""
    vtk_entity(t::AbstractEntityType)

Convert an `AbstractEntityType` into a `VTKCellType`. To find the correspondance, browse the `WriteVTK`
package AND check the Doxygen (for numbering) : https://vtk.org/doc/nightly/html/classvtkTriQuadraticHexahedron.html
"""
function vtk_entity(t::AbstractEntityType)
    error("Entity type $t doesn't have a VTK correspondance")
end

vtk_entity(::Node_t) = VTKCellTypes.VTK_VERTEX
vtk_entity(::Bcube.Bar2_t) = VTKCellTypes.VTK_LINE
vtk_entity(::Bcube.Bar3_t) = VTKCellTypes.VTK_LAGRANGE_CURVE
vtk_entity(::Bcube.Bar4_t) = VTKCellTypes.VTK_LAGRANGE_CURVE
#vtk_entity(::Bar5_t)    = error("undefined")
vtk_entity(::Bcube.Tri3_t) = VTKCellTypes.VTK_TRIANGLE
vtk_entity(::Bcube.Tri6_t) = VTKCellTypes.VTK_LAGRANGE_TRIANGLE #VTK_QUADRATIC_TRIANGLE
vtk_entity(::Bcube.Tri9_t) = error("undefined")
vtk_entity(::Bcube.Tri10_t) = VTKCellTypes.VTK_LAGRANGE_TRIANGLE
#vtk_entity(::Tri12_t)   = error("undefined")
vtk_entity(::Bcube.Quad4_t) = VTKCellTypes.VTK_QUAD
vtk_entity(::Bcube.Quad8_t) = VTKCellTypes.VTK_QUADRATIC_QUAD
vtk_entity(::Bcube.Quad9_t) = VTKCellTypes.VTK_BIQUADRATIC_QUAD
vtk_entity(::Bcube.Quad16_t) = VTKCellTypes.VTK_LAGRANGE_QUADRILATERAL
vtk_entity(::Bcube.Tetra4_t) = VTKCellTypes.VTK_TETRA
vtk_entity(::Bcube.Tetra10_t) = VTKCellTypes.VTK_QUADRATIC_TETRA
vtk_entity(::Bcube.Penta6_t) = VTKCellTypes.VTK_WEDGE
vtk_entity(::Bcube.Hexa8_t) = VTKCellTypes.VTK_HEXAHEDRON
vtk_entity(::Bcube.Pyra5_t) = VTKCellTypes.VTK_PYRAMID
#vtk_entity(::Hexa27_t) = VTK_TRIQUADRATIC_HEXAHEDRON # NEED TO CHECK NODE NUMBERING : https://vtk.org/doc/nightly/html/classvtkTriQuadraticHexahedron.html

vtk_entity(::Line, ::Val{Degree}) where {Degree} = VTKCellTypes.VTK_LAGRANGE_CURVE
vtk_entity(::Square, ::Val{Degree}) where {Degree} = VTKCellTypes.VTK_LAGRANGE_QUADRILATERAL
vtk_entity(::Triangle, ::Val{Degree}) where {Degree} = VTKCellTypes.VTK_LAGRANGE_TRIANGLE

get_vtk_name(c::VTKCellType) = Val(Symbol(c.vtk_name))
const VTK_LAGRANGE_QUADRILATERAL = get_vtk_name(VTKCellTypes.VTK_LAGRANGE_QUADRILATERAL)
const VTK_LAGRANGE_TRIANGLE = get_vtk_name(VTKCellTypes.VTK_LAGRANGE_TRIANGLE)

"""
Return the node numbering of the node designated by its position in the x and y direction.

See https://www.kitware.com/modeling-arbitrary-order-lagrange-finite-elements-in-the-visualization-toolkit/.
"""
function _point_index_from_IJK(t::Val{:VTK_LAGRANGE_QUADRILATERAL}, degree, i, j)
    1 + _point_index_from_IJK_0based(t, degree, i - 1, j - 1)
end

# see : https://github.com/Kitware/VTK/blob/675adbc0feeb3f62730ecacb2af87917af124543/Filters/Sources/vtkCellTypeSource.cxx
#        https://github.com/Kitware/VTK/blob/265ca48a79a36538c95622c237da11133608bbe5/Common/DataModel/vtkLagrangeQuadrilateral.cxx#L558
function _point_index_from_IJK_0based(::Val{:VTK_LAGRANGE_QUADRILATERAL}, degree, i, j)
    # 0-based algo
    ni = degree
    nj = degree
    ibnd = ((i == 0) || (i == ni))
    jbnd = ((j == 0) || (j == nj))
    # How many boundaries do we lie on at once?
    nbnd = (ibnd > 0 ? 1 : 0) + (jbnd > 0 ? 1 : 0)

    if nbnd == 2 # Vertex DOF
        return (i > 0 ? (j > 0 ? 2 : 1) : (j > 0 ? 3 : 0))
    end

    offset = 4
    if nbnd == 1 # Edge DOF
        ibnd == 0 && (return (i - 1) + (j > 0 ? degree - 1 + degree - 1 : 0) + offset)
        jbnd == 0 &&
            (return (j - 1) + (i > 0 ? degree - 1 : 2 * (degree - 1) + degree - 1) + offset)
    end

    offset = offset + 2 * (degree - 1 + degree - 1)
    # Face DOF
    return (offset + (i - 1) + (degree - 1) * (j - 1))
end

function _vtk_lagrange_node_index_vtk_to_bcube(shape::AbstractShape, degree)
    return 1:length(get_coords(FunctionSpace(Lagrange(), degree), shape))
end

"""
Bcube node numbering -> VTK node numbering (in a cell)
"""
function _vtk_lagrange_node_index_bcube_to_vtk(shape::Union{Square, Cube}, degree)
    n = _get_num_nodes_per_dim(
        QuadratureRule(shape, Quadrature(QuadratureUniform(), Val(degree))),
    )
    IJK = CartesianIndices(ntuple(i -> 1:n[i], length(n)))

    vtk_cell_name = get_vtk_name(vtk_entity(shape, Val(degree)))

    # We loop over all the nodes designated by (i,j,k) and find the corresponding
    # VTK index for each node.
    index = map(vec(IJK)) do ijk
        _point_index_from_IJK(vtk_cell_name, degree, Tuple(ijk)...)
    end
    return index
end

"""
VTK node numbering (in a cell) -> Bcube node numbering
"""
function _vtk_lagrange_node_index_vtk_to_bcube(shape::Union{Square, Cube}, degree)
    return invperm(_vtk_lagrange_node_index_bcube_to_vtk(shape, degree))
end

function _vtk_coords_from_lagrange(shape::AbstractShape, degree)
    return get_coords(FunctionSpace(Lagrange(), degree), shape)
end

"""
Coordinates of the nodes in the VTK cell, ordered as expected by VTK.
"""
function _vtk_coords_from_lagrange(shape::Union{Square, Cube}, degree)
    fs = FunctionSpace(Lagrange(), degree)
    c = get_coords(fs, shape)
    index = _vtk_lagrange_node_index_vtk_to_bcube(shape, degree)
    return c[index]
end

"""
    write_vtk_lagrange(
        basename::String,
        vars::Dict{String, F},
        mesh::AbstractMesh,
        it::Integer = -1,
        time::Real = 0.0;
        mesh_degree::Integer = 1,
        discontinuous::Bool = true,
        functionSpaceType::AbstractFunctionSpaceType = Lagrange(),
        collection_append::Bool = false,
        vtk_kwargs...,
    ) where {F <: AbstractLazy}

    write_vtk_lagrange(
        basename::String,
        vars::Dict{String, F},
        mesh::AbstractMesh,
        U_export::AbstractFESpace,
        it::Integer = -1,
        time::Real = 0.0;
        collection_append::Bool = false,
        vtk_kwargs...,
    ) where {F <: AbstractLazy}

Write the provided FEFunction/MeshCellData/CellFunction on a mesh of order `mesh_degree` with the nodes
defined by the `functionSpaceType` (only `Lagrange{:Uniform}` supported for now). The boolean `discontinuous`
indicate if the node values should be discontinuous or not.

`vars` is a dictionnary of variable name => Union{FEFunction,MeshCellData,CellFunction} to write.

# Dev notes
- in order to write an ASCII file, you must pass both `ascii = true` and `append = false`
- `collection_append` is not named `append` to enable passing correct `kwargs` to `vtk_grid`
- remove (once fully validated) : `write_vtk_discontinuous`
"""
function write_vtk_lagrange(
    basename::String,
    vars::Dict{String, F},
    mesh::AbstractMesh,
    it::Integer = -1,
    time::Real = 0.0;
    mesh_degree::Integer = 1,
    discontinuous::Bool = true,
    functionSpaceType::AbstractFunctionSpaceType = Lagrange(),
    collection_append::Bool = false,
    vtk_kwargs...,
) where {F <: AbstractLazy}
    U_export = TrialFESpace(
        FunctionSpace(functionSpaceType, mesh_degree),
        mesh;
        isContinuous = !discontinuous,
    )
    write_vtk_lagrange(
        basename,
        vars,
        mesh,
        U_export,
        it,
        time;
        collection_append,
        vtk_kwargs...,
    )
end

function write_vtk_lagrange(
    basename::String,
    vars::Dict{String, F},
    mesh::AbstractMesh,
    U_export::AbstractFESpace,
    it::Integer = -1,
    time::Real = 0.0;
    collection_append = false,
    vtk_kwargs...,
) where {F <: AbstractLazy}
    # FE space stuff
    fs_export = get_function_space(U_export)
    @assert get_type(fs_export) <: Lagrange "Only FunctionSpace of type Lagrange are supported for now"
    degree_export = get_degree(fs_export)
    dhl_export = Bcube._get_dhl(U_export)
    nd = get_ndofs(dhl_export)

    # extract all `MeshCellData` in `vars` as this type of variable
    # will not be interpolated to nodes and will be written with
    # the `VTKCellData` attribute.
    vars_cell = filter(((k, v),) -> v isa MeshData{<:CellData}, vars)
    vars_point = filter(((k, v),) -> k ∉ keys(vars_cell), vars)

    # Get ncomps and type of each `point` variable
    type_dim = map(var -> get_return_type_and_codim(var, mesh), values(vars_point))

    # VTK stuff
    coords_vtk = zeros(spacedim(mesh), nd)
    values_vtk = map(((_t, _d),) -> zeros(_t, _d..., nd), type_dim)
    cells_vtk = MeshCell[]
    sizehint!(cells_vtk, ncells(mesh))
    nodeweigth_vtk = zeros(nd)

    # Loop over mesh cells
    for cinfo in DomainIterator(CellDomain(mesh))
        # Cell infos
        icell = cellindex(cinfo)
        ctype = celltype(cinfo)
        _shape = shape(ctype)

        # Get Lagrange dofs/nodes coordinates in ref space (Bcube order)
        coords_bcube = get_coords(fs_export, _shape)

        # Get VTK <-> BCUBE dofs/nodes mapping
        v2b = _vtk_lagrange_node_index_vtk_to_bcube(_shape, degree_export)

        # Add nodes coordinates to coords_vtk
        for (iloc, ξ) in enumerate(coords_bcube)
            # Global number of the node
            iglob = get_dof(dhl_export, icell, 1, iloc)

            # Create CellPoint (in reference domain)
            cpoint = CellPoint(ξ, cinfo, ReferenceDomain())

            # Map it to the physical domain and assign to VTK
            coords_vtk[:, iglob] .= get_coords(change_domain(cpoint, PhysicalDomain()))
            nodeweigth_vtk[iglob] += 1.0

            # Evaluate all vars on this node
            for (ivar, var) in enumerate(values(vars_point))
                _var = Bcube.materialize(var, cinfo)
                values_vtk[ivar][:, iglob] .+= Bcube.materialize(_var, cpoint)
            end
        end

        # Create the VTK cells with the correct ordering
        iglobs = get_dof(dhl_export, icell)
        push!(
            cells_vtk,
            MeshCell(vtk_entity(_shape, Val(degree_export)), view(iglobs, v2b)),
        )
    end

    # Averaging contribution at nodes :
    # For discontinuous variables written as a continuous
    # field, node averaging is necessary because each adjacent
    # cells gives different interpolated values at nodes, then
    # a averaging step is used to define a unique node value.
    # This step has no effect on node values for :
    # - discontinuous variables written as a discontinous
    #   fields as nodes are duplicated in this case and
    #   there is one cell contribution for each duplicated
    #   node (nodeweigth_vtk[i] is equal to 1, ∀i).
    # - continous variables written as a continous fields
    #   because interpolated values are all equals and
    #   averaging gives the same value.
    for val in values_vtk
        for i in eachindex(nodeweigth_vtk)
            val[:, i] .= val[:, i] ./ nodeweigth_vtk[i]
        end
    end

    # Write VTK file
    pvd = paraview_collection(basename; append = collection_append)
    new_name = _build_fname_with_iterations(basename, it)
    vtkfile = vtk_grid(new_name, coords_vtk, cells_vtk; vtk_kwargs...)

    for (varname, value) in zip(keys(vars_point), values_vtk)
        vtkfile[varname, VTKPointData()] = value
    end
    for (varname, value) in zip(keys(vars_cell), values(vars_cell))
        vtkfile[varname, VTKCellData()] = get_values(value)
    end

    pvd[float(time)] = vtkfile
    vtk_save(pvd)
end

"""
Append the number of iteration (if positive) to the basename
"""
function _build_fname_with_iterations(basename::String, it::Integer)
    it >= 0 ? @sprintf("%s_%08i", basename, it) : basename
end

function Bcube.write_file(
    ::VTKIoHandler,
    filepath::String,
    mesh::AbstractMesh,
    U_export::AbstractFESpace,
    data = nothing,
    it::Integer = -1,
    time::Real = 0.0;
    collection_append = false,
    kwargs...,
)

    # Remove extension from filename
    basename = first(splitext(filepath))

    # Just write the mesh if `data` is `nothing`
    if isnothing(data)
        write_vtk(basename, mesh)
        return
    end

    # We don't use FlowSolution names in VTK, so we flatten everything
    _data = data
    if valtype(data) <: Dict
        _keys = map(d -> collect(keys(d)), values(data))
        _keys = vcat(_keys...)
        _values = map(d -> collect(values(d)), values(data))
        _values = vcat(_values...)
        _data = Dict(_keys .=> _values)
    end

    # Write !
    write_vtk_lagrange(
        basename,
        _data,
        mesh,
        U_export,
        it,
        time;
        collection_append,
        kwargs...,
    )
end
