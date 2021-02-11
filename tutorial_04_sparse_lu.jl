using Random
using SparseArrays
using LinearAlgebra

function path_theorems()
    Random.seed!(1)
    # ^ Resets the random number generator. This ensures we get the same random
    # matrix each time we run this function. You can create additional exercise
    # material by changing the above seed.

    A = Matrix(5I + sprand(5,5,0.25))

    println("Sparsity pattern of A:")
    display(2 .* (A .!= 0))
    println()

    println("Sparsity pattern of A^2:")
    display((A .!= 0) .+ (A^2 .!= 0))
    println()

    println("Sparsity pattern of inv(A):")
    display((A .!= 0) .+ (inv(A) .!= 0))
    println()

    println("Sparsity pattern of LU factorisation:")
    L,U = lu(A, Val(false))
    display((A .!= 0) .+ (L+U .!= 0))
    println()
    println("(2 = original nonzero, 1 = fill-in)")
end
