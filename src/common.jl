struct VTKIoHandler <: Bcube.AbstractIoHandler end

Bcube._filename_to_handler(::Union{Val{:pvd}, Val{:vtk}, Val{:vtu}}) = VTKIoHandler()