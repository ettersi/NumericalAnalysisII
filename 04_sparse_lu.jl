using PyPlot
using LinearAlgebra
using SparseArrays
using AMD

using Printf

function upleft_arrow_matrix(n)
    A = spdiagm(0 => fill(n,n))
    A[2:end, 1] .= A[1, 2:end] .= 1.0
    return A
end

function downright_arrow_matrix(n)
    A = spdiagm(0 => fill(n,n))
    A[1:end-1, end] .= A[end, 1:end-1] .= 1.0
    return A
end

function lu_benchmark()
    n = 2_000
    println("Runtime of LU factorisation:")
    for matrix in (upleft_arrow_matrix, downright_arrow_matrix)
        A = matrix(n)
        time = @elapsed(ldlt(A, perm=1:n))
        # LDLt is the symmetric version of the LU factorisation.
        # Ignore the `perm = 1:n` argument.

        @printf("%22s: %6f seconds\n", matrix, time)
    end
end

function lu_structures()
    n = 5
    for matrix in (upleft_arrow_matrix, downright_arrow_matrix)
        A = matrix(n)
        L,U = lu(Matrix(A), Val(false))  # Ignore the `Val(false)` argument
        println(matrix, ":")
        println("L:")
        display(L .!= 0)
        println("U:")
        display(U .!= 0)
        println()
    end
end
