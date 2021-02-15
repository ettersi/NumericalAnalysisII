using PyPlot
using LinearAlgebra
using SparseArrays
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
    n = 1000
    println("Runtime of LU factorisation:")
    for matrix in (
        upleft_arrow_matrix,
        downright_arrow_matrix,
    )
        A = matrix(n)
        time = @elapsed(ldlt(A, perm=1:n))
        # LDLt is the symmetric version of the LU factorisation.
        # Ignore the `perm = 1:n` argument for now.

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

function matrix()
    p = [4,2,1,3]
    A = [
        1 0 0 1
        1 1 0 0
        1 0 1 0
        0 0 0 1
    ]
    println("A:")
    display(A)
    println()
    println("P*A*P':")
    display(A[p,p])
    println()
end

function permutation()
    p = [4,2,1,3]
    q = [3,2,4,1]
    @show p[q]
end



using Combinatorics

function check_fillin()
    # Define matrix
    Random.seed!(42)
    A = [
           9   rand() rand()    0
        rand()    9      0   rand()
        rand()    0      9   rand()
           0   rand() rand()    9
    ]

    # Loop over all permutations
    for p in permutations(1:4)
        for q in permutations(1:4)
            # Skip if top-left entry is 0
            if A[p[1],q[1]] == 0; continue; end

            # Check that the number of fill-in entries is exactly 2
            L,U = lu(A[p,q], Val(false))
            n_fill = sum((L+U) .!= 0) - sum(A .!= 0)
            if n_fill != 2
                error("Permutations p = ", p, " and q = ", q, " lead to ", n_fill, " fill-in entries")
            end
        end
    end
    println("All permutations lead to exactly two fill-in entries")
end



laplacian_1d(n) = (n+1)^2*Tridiagonal(
    fill( 1.0,n-1), # subdiagonal
    fill(-2.0,n),   # diagonal
    fill( 1.0,n-1)  # superdiagonal
)

function laplacian_2d(n)
    Δ = sparse(laplacian_1d(n))
    Id = sparse(I,n,n)
    return kron(Δ,Id) + kron(Id,Δ)
end

"""
    nested_dissection(n) -> p

Compute the nested dissection permutation for `laplacian_2d(n)`.

Note that `n` must be of the form `n = 2^k - 1`.

For convenience, this function works with two-dimensional indices `i1,i2`
rather than the one-dimensional index `i = i1 + n*(i2-1)`. The output of this
function therefore needs to be converted using the
`matrix_to_vector_indices()` function provided below before it can be used to
permute a matrix.
"""
function nested_dissection(n)
    @assert isodd(n)
    n == 1 && return [(1,1)]

    nn = n÷2
    pp = nested_dissection(nn)
    return [
        [(i[1]     , i[2]     ) for i in pp];  # top left quadrant
        [(i[1]+nn+1, i[2]     ) for i in pp];  # bottom left quadrant
        [(i[1]     , i[2]+nn+1) for i in pp];  # top right quadrant
        [(i[1]+nn+1, i[2]+nn+1) for i in pp];  # bottom right quadrant
        [(nn+1,      i2) for i2 in 1:nn];      # left horizontal separator
        [(nn+1, nn+1+i2) for i2 in 1:nn];      # right horizontal separator
        [(i1, nn+1) for i1 in 1:n]             # vertical separator
    ]
end

matrix_to_vector_indices(n,p) = [i[1] + n*(i[2]-1) for i in p]

"""
    runtimes(n = 127)

Print the runtimes of LU factorisation applied to `P*laplacian_2d(n)*P'` for
`P in [I, nested_dissection]`.

Note that `n` must be of the form `n = 2^k - 1`.
"""
function runtimes(n = 127)
    Δ = -laplacian_2d(n)

    t = @elapsed( ldlt(Δ; perm = 1:n^2) )
    @printf("         Original: %5.3f seconds\n", t)

    p = matrix_to_vector_indices(n,nested_dissection(n))
    t = @elapsed( ldlt(Δ; perm = p) )
    @printf("Nested dissection: %5.3f seconds\n", t)

    t = @elapsed( ldlt(Δ) )
    @printf("          Default: %5.3f seconds\n", t)
end

"""
    sparsity_pattern(n, perm = "")

Plot the sparsity pattern of `L+U`, where `L`,`U` is the LU factorisation of
`Δ = laplacian_2d(n)`.

`perm in ["","nd"]` denotes the permutation to apply to `Δ` before
the factorisation.

Note that for `perm == "nd"`, `n` must be of the form `n = 2^k - 1`.
"""
function sparsity_pattern(n,perm = "")
    A = -laplacian_2d(n)
    if perm == "nd"
        p = matrix_to_vector_indices(n, nested_dissection(n))
        A = A[p,p]
    end

    L,U = lu(Matrix(A), Val(false)) # `Val(false)` disables pivoting for stability

    clf()
    spy(L.+U, marker="s", c="r", ms=3e2*n^(-2))
    spy(  A , marker="s", c="k", ms=3e2*n^(-2))
    display(gcf())
end

