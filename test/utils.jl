#------------------------------------------------#
#  Functions for tests based on ReferenceTest.jl #
#------------------------------------------------#

"""
    refpath(filename::String)

Return the full path of the reference file named `filename`.
It assumes that the file is located in the `test/references/`
directory of the current project.
"""
refpath(filename::String) = joinpath(@__DIR__, "references/", filename)

function compare_arrays(a, b; rtol, atol)
    # We mainly rely on `isapprox` (when everything's fine)
    # and only perform additionnal checks to explain failure if needed
    check_failed = any((!isapprox).(a, b; atol, rtol))

    if check_failed
        println("Comparison failure info: atol=$atol, rtol=$rtol")
        # We perform the check on all entries, without stopping on the first failure
        # (except for Inf/NaN)
        for (x, y) in zip(a, b)
            # Code below adapted from "isapprox"
            if !isfinite(x) || !isfinite(y)
                println("isfinite check failed")
                break
            end
            if isnan(x) || isnan(y)
                println("isnan check failed")
                break
            end

            x′, y′ = promote(x, y) # to avoid integer overflow
            n = norm(x - y)
            r = rtol * max(norm(x′), norm(y′))
            if n > max(atol, r)
                println(
                    "x=$x, y=$y, norm(x-y)=$n, rtol*max(norm(x),norm(y))=$(rtol*max(norm(x),norm(y)))",
                )
            end
        end
        return false
    else
        return true
    end
end

"""
    compare(
        d1::Dict{String, <:AbstractArray},
        d2::Dict{String, <:AbstractArray},
        rtol,
        atol,
    )

Return `true` if dictionaries `d1` and `d2` have the
same keys and if these keys are associated with arrays
of equal values, with a relative tolerance `rtol` and
an absolute tolerance `atol`.
"""
function compare(
    d1::Dict{String, <:AbstractArray},
    d2::Dict{String, <:AbstractArray};
    atol,
    rtol,
)
    # From the official doc for `Dict`:
    # -
    # "When the values are stored internally in a hash table,
    # as is the case for Dict, the order in which they are returned may vary".
    # -
    # Then `keys(d1) == keys(d2)` is not recommended
    # and we do the following two tests instead :
    length(keys(d1)) != length(keys(d2)) && return false
    any(x -> !haskey(d1, x), keys(d2)) && return false

    return all(keys(d1)) do key
        _x = d1[key]
        _y = d2[key]
        return compare_arrays(_x, _y; rtol, atol)
    end
end
compare(; atol::Real = 1.0e-12, rtol::Real = 1.0e-12) = (a, b) -> compare(a, b; atol, rtol)

"""
    compare_recursive(
        d1::Dict,
        d2::Dict,
        rtol,
        atol,
    )

Return `true` if dictionaries `d1` and `d2` have the same structure and content.
Performs recursive comparison: nested dictionaries are compared recursively,
and arrays are compared using element-wise approximate equality.
"""
function compare_recursive(d1::Dict, d2::Dict; rtol, atol)
    # Check if dictionaries have the same keys
    length(keys(d1)) != length(keys(d2)) && return false
    any(x -> !haskey(d1, x), keys(d2)) && return false

    return all(keys(d1)) do key
        v1 = d1[key]
        v2 = d2[key]

        if v1 isa Dict && v2 isa Dict
            # If both values are dicts, compare recursively
            return compare_recursive(v1, v2; rtol, atol)

        elseif v1 isa AbstractArray && v2 isa AbstractArray
            # If both values are arrays, use approximate comparison
            return compare_arrays(v1, v2; rtol, atol)
        else
            error("Type not supported")
        end
    end
end
function compare_recursive(; atol::Real = 1.0e-12, rtol::Real = 1.0e-12)
    (a, b) -> compare_recursive(a, b; atol, rtol)
end

function test_ref(filepath::String; rtol = eps(), atol = eps())
    filename = basename(filepath)
    d1 = parse_vtk_xml(refpath(filename))
    d2 = parse_vtk_xml(filepath)
    @test compare_recursive(d1, d2; rtol, atol)
end