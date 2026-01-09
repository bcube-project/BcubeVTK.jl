module BcubeVTK
using Bcube
using WriteVTK
using Printf

import Bcube:
    AbstractEntityType,
    AbstractMesh,
    AbstractShape,
    AbstractLazy,
    AbstractFunctionSpaceType,
    AbstractFESpace,
    Line,
    Square,
    Triangle,
    Cube,
    Tetra,
    Prism,
    connectivities_indices,
    get_function_space,
    get_type,
    Lagrange,
    CellData,
    MeshData,
    get_return_type_and_codim,
    DomainIterator,
    cellindex,
    celltype,
    shape,
    _get_num_nodes_per_dim,
    get_dof,
    CellPoint,
    ReferenceDomain,
    PhysicalDomain,
    change_domain,
    cells,
    Node_t

include("common.jl")
include("write.jl")
# include("read.jl") # to be implemented with ReadVTK
end
