using Random
using SparseArrays
using LinearAlgebra

function example_sparse_matrix()
    # Make sure we get the same result each time we run this function
    Random.seed!(42)

    # Type `?sprand` on the REPL to get the documentation for this function
    return Matrix(5I + sprand(5,5,0.5))
end
