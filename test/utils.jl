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

    for key in keys(d1)
        _x = d1[key]
        _y = d2[key]

        # We mainly rely on `isapprox` (when everything's fine)
        # and only perform additionnal checks to explain failure if needed
        check_failed = any((!isapprox).(d1[key], d2[key]; atol, rtol))

        if check_failed
            println("Comparison failure info: atol=$atol, rtol=$rtol")
            # We perform the check on all entries, without stopping on the first failure
            # (except for Inf/NaN)
            for (x, y) in zip(_x, _y)
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
        end
    end
    return true
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

        # If both values are dicts, compare recursively
        if v1 isa Dict && v2 isa Dict
            return compare_recursive(v1, v2; rtol, atol)

            # If both values are arrays, use approximate comparison
        elseif v1 isa AbstractArray && v2 isa AbstractArray
            return compare(v1, v2; rtol, atol)
        else
            error("Type not supported")
        end
    end
end
function compare_recursive(; atol::Real = 1.0e-12, rtol::Real = 1.0e-12)
    (a, b) -> compare_recursive(a, b; atol, rtol)
end

"""
    test_ref(filename_ref::String, data, comp::Function = compare())

Test that the values `data` with reference `filename_ref` (stored in `./test/references`)
using equality test strategy given by `comp`.

By default, `comp=compare()` assumes `data` is a subtype of
`Dict{String, <:AbstractArray}`.

If `data` is not a valid subtype for the defaut test strategy, one could :
- defined the method `_as_dict(data)` to convert `data` to a valid subtype.
- provide another test function with the keyword argument `comp`.
The function must have the same signature as `Base.isequal` function.
"""
function test_ref(filename_ref::String, data, comp::Function = compare())
    println("Testing against $(filename_ref)")
    @test_reference refpath(filename_ref) _as_dict(data) by = comp
end
_as_dict(a::Dict) = a
_as_dict(a::AbstractArray) = Dict("array" => a)
_as_dict(a::AbstractSparseMatrix) = Dict(zip("sparse_" .* ("I", "J", "V"), findnz(a)))