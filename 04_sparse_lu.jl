using PyPlot
using LinearAlgebra
using SparseArrays
using AMD

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

function runtimes()
    n = 2^8-1
    Δ = -laplacian_2d(n)

    # The below code uses `cholesky()` instead of `lu()`.
    # `cholesky()` is the symmetric version of `lu()`.
    # We have to use `cholesky()` here because only this functions
    # allows us to prescribe the permutation.
    # For sparse matrices, `lu()` will always determine a permutation
    # itself and there is no way for us to turn that off.

    t = @elapsed( cholesky(Δ; perm = 1:n^2) )
    println("         Original: ", round(t, sigdigits=3), " sec")

    p = matrix_to_vector_indices(n,nested_dissection(n))
    t = @elapsed( cholesky(Δ; perm = p) )
    println("Nested dissection: ", round(t, sigdigits=3), " sec")

    t = @elapsed( cholesky(Δ) )
    println("              AMD: ", round(t, sigdigits=3), " sec")
end

function sparsity(n,perm = "")
    A = -laplacian_2d(n)
    if perm == "nd"
        p = matrix_to_vector_indices(n, nested_dissection(n))
        A = A[p,p]
    end
    if perm == "amd"
        p = amd(A)
        A = A[p,p]
    end

    # We use dense `lu()` here because sparse `cholesky()` does not
    # provide a convenient way to extract the factors.
    # The `Val(false)` argument turns off pivoting in the factorisation.
    L,U = lu(Matrix(A), Val(false))

    clf()
    spy(L.+U, marker="s", c="r", ms=3e2*n^(-2))
    spy(  A , marker="s", c="k", ms=3e2*n^(-2))
end
