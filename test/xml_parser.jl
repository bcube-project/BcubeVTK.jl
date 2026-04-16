using XML

"""
    parse_vtk_xml(filepath::AbstractString)

Parse a VTK XML file (.vtu) and return the solution as a dictionary of dictionaries.

Returns a Dict with keys:
- "Points": matrix of point coordinates (n_dim_spa × n_points)
- "Cells": Dict with "connectivity", "offsets", "types"
- "PointData": Dict of point data arrays
- "CellData": Dict of cell data arrays (if present)
"""
function parse_vtk_xml(filepath::AbstractString)
    content = read(filepath, String)
    doc = XML.parse(content, XML.LazyNode)

    result = Dict{String, Any}()

    # Navigate to Piece element
    # Structure: VTKFile -> UnstructuredGrid -> Piece
    # Skip XML declaration node (first child with tag=nothing)
    vtkfile = nothing
    for child in doc.children
        if child.tag == "VTKFile"
            vtkfile = child
            break
        end
    end
    if vtkfile === nothing
        error("VTKFile element not found")
    end

    unstructured_grid = find_child_by_tag(vtkfile, "UnstructuredGrid")
    if unstructured_grid === nothing
        error("UnstructuredGrid element not found")
    end

    piece = find_child_by_tag(unstructured_grid, "Piece")
    if piece === nothing
        error("Piece element not found")
    end

    # Parse Points
    points_section = find_child_by_tag(piece, "Points")
    points_array = find_child_by_tag(points_section, "DataArray")
    n_dims = parse(Int, get_attr(points_array, "NumberOfComponents"))
    points_data = parse.(Float64, split(get_text_content(points_array)))
    n_points = div(length(points_data), n_dims)
    points_matrix = reshape(points_data, n_dims, n_points)
    result["Points"] = points_matrix

    # Parse Cells
    cells_section = find_child_by_tag(piece, "Cells")
    names = ("connectivity", "offsets", "types")
    cells_dict = Dict(
        name =>
            parse.(Int, split(get_text_content(find_child_by_name(cells_section, name)))) for name in names
    )
    result["Cells"] = cells_dict

    # Parse PointData
    point_data_section = find_child_by_tag(piece, "PointData")
    if point_data_section !== nothing && !isempty(point_data_section.children)
        name_and_values = parse_data_section(point_data_section)
        result["PointData"] = Dict(name => array for (name, array) in name_and_values)
    end

    # Parse CellData if present
    cell_data_section = find_child_by_tag(piece, "CellData")
    if cell_data_section !== nothing && !isempty(cell_data_section.children)
        name_and_values = parse_data_section(point_data_section)
        result["CellData"] = Dict(name => array for (name, array) in name_and_values)
    end

    return result
end

# Helper function to get text content from a node
function get_text_content(node)
    if node.value !== nothing
        return String(node.value)
    end
    # If the node has children, get the text from the first child
    if !isempty(node.children)
        return String(node.children[1].value)
    end
    return ""
end

# Helper function to get attribute value
function get_attr(node, attr_name)
    attrs = attributes(node)
    if haskey(attrs, attr_name)
        return attrs[attr_name]
    end
    return nothing
end

# Helper function to find child by tag name
function find_child_by_tag(node, tag_name)
    for child in node.children
        if child.tag == tag_name
            return child
        end
    end
    return nothing
end

function find_child_by_name(node, name)
    for child in node.children
        if get_attr(child, "Name") == name
            return child
        end
    end
    return nothing
end

function parse_data_section(data_section)
    dArrays = filter(child -> child.tag == "DataArray", data_section.children)
    name_and_values = map(dArrays) do dArray
        name = get_attr(dArray, "Name")
        n_components = parse(Int, get_attr(dArray, "NumberOfComponents"))
        data_content = get_text_content(dArray)
        data_values = parse.(Float64, split(data_content))
        n_points_data = length(data_values) ÷ n_components
        vals = if n_components == 1
            data_values
        else
            reshape(data_values, n_components, n_points_data)
        end
        return name, vals
    end
    return name_and_values
end