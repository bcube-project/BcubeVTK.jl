@testset "vtk" begin
    @testset "write_vtk" begin
        mesh = rectangle_mesh(10, 10; xmax = 2π, ymax = 2π)
        val_sca = PhysicalFunction(x -> cos(x[1]) * sin(x[2]))
        val_vec = PhysicalFunction(x -> SA[cos(x[1]), sin(x[2])])
        meshdata = MeshCellData(Bcube.get_cell_centers(mesh))
        basename = "write_vtk_rectangle"
        write_file(
            joinpath(tempdir, basename * ".pvd"),
            mesh,
            Dict("u" => val_sca, "v" => val_vec, "meshdata" => meshdata),
            1,
            0.0;
            ascii = true,
            compress = false,
            append = false,
        )
        fname = BcubeVTK._build_fname_with_iterations(basename, 1) * ".vtu"
        test_ref(joinpath(tempdir, fname))

        basename = "write_vtk_mesh"
        write_file(
            joinpath(tempdir, basename * ".pvd"),
            basic_mesh();
            ascii = true,
            compress = false,
            append = false,
        )
        fname = basename * ".vtu"
        test_ref(joinpath(tempdir, fname))
    end

    @testset "write_vtk_lagrange" begin
        mesh = rectangle_mesh(6, 7; xmin = -1.0, xmax = 1.0, ymin = -1.0, ymax = 1.0)
        u = FEFunction(TrialFESpace(FunctionSpace(:Lagrange, 4), mesh))
        projection_l2!(u, PhysicalFunction(x -> x[1]^2 + x[2]^2), mesh)

        vars = Dict("u" => u, "grad_u" => ∇(u))

        # bmxam: for some obscur reason, order 3 and 5 lead to different sha1sum
        # when running in standard mode or in test mode...
        for mesh_degree in (1, 2, 4)
            basename = "write_vtk_lagrange_deg$(mesh_degree)"
            write_vtk_lagrange(
                joinpath(tempdir, basename),
                vars,
                mesh;
                mesh_degree,
                discontinuous = false,
                vtkversion = v"1.0",
                ascii = true,
                compress = false,
                append = false,
            )

            # Check
            fname = basename * ".vtu"
            test_ref(joinpath(tempdir, fname))
        end

        # add var MeshCellData :
        quad = Quadrature(4)
        dΩ = Measure(CellDomain(mesh), quad)
        vars["umean"] = cell_mean(u, dΩ)
        basename = "write_vtk_lagrange_deg4_with_mean"
        write_vtk_lagrange(
            joinpath(tempdir, basename),
            vars,
            mesh;
            mesh_degree = 4,
            discontinuous = false,
            vtkversion = v"1.0",
            ascii = true,
            compress = false,
            append = false,
        )

        # Check
        fname = basename * ".vtu"
        test_ref(joinpath(tempdir, fname))

        basename = "write_vtk_lagrange_deg4_dg_with_mean"
        write_vtk_lagrange(
            joinpath(tempdir, basename),
            vars,
            mesh;
            mesh_degree = 4,
            discontinuous = true,
            vtkversion = v"1.0",
            ascii = true,
            compress = false,
            append = false,
        )

        # Check
        fname = basename * ".vtu"
        test_ref(joinpath(tempdir, fname))
    end

    @testset "write_file hexa" begin
        mesh = one_cell_mesh(:hexa)
        U = TrialFESpace(FunctionSpace(:Lagrange, 1), mesh)
        u = FEFunction(U, collect(1:nnodes(mesh)))
        basename = "write_vtk_lagrange_hexa_deg1.vtu"
        write_file(
            joinpath(tempdir, basename),
            mesh,
            Dict("u" => u);
            vtkversion = v"1.0",
            ascii = true,
            compress = false,
            append = false,
        )

        # Check
        fname = basename
        test_ref(joinpath(tempdir, fname))
    end

    @testset "write_file tetra" begin
        mesh = one_cell_mesh(:tetra)
        U = TrialFESpace(FunctionSpace(:Lagrange, 1), mesh)
        u = FEFunction(U, collect(1:nnodes(mesh)))
        basename = "write_vtk_lagrange_tetra_deg1.vtu"
        write_file(
            joinpath(tempdir, basename),
            mesh,
            Dict("u" => u);
            vtkversion = v"1.0",
            ascii = true,
            compress = false,
            append = false,
        )

        # Check
        fname = basename
        test_ref(joinpath(tempdir, fname))
    end

    @testset "misc" begin
        mesh = one_cell_mesh(:line)
        f = PhysicalFunction(x -> 1.0)
        d = Dict(
            "Gas" => Dict("Temperature" => f, "u" => f),
            "Droplet" => Dict("v" => f, "Temperature" => f),
        )
        basename = "write_vtk_dict_of_dict.vtu"
        write_file(
            joinpath(tempdir, basename),
            mesh,
            d;
            vtkversion = v"1.0",
            ascii = true,
            compress = false,
            append = false,
        )

        # Check
        fname = basename
        test_ref(joinpath(tempdir, fname))
        # Reading VTK is not supported for now
        # result = read_file(tmppath)
        # @assert "Gas_Temperature" ∈ keys(result.d)
        # @assert "Droplet_Temperature" ∈ keys(result.d)
        # @assert "Temperature" ∉ keys(result.d)
        # @assert "u" ∈ keys(result.d)
        # @assert "v" ∈ keys(result.d)
    end
end
