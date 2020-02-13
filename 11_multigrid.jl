using LinearAlgebra
using SparseArrays
using PyPlot

laplacian(n) = (n+1)^2 * sparse(Tridiagonal(
    fill( 1.0,n-1),
    fill(-2.0,n),
    fill( 1.0,n-1)
))

function pmatrix(n)
    # I recommend you do *not* try to understand the code for this function.
    # Instead, have a look at `Matrix(pmatrix(3))` to see what this function
    # does.
    p = Vector(1:3:3*n+1)
    i = vec( (-1:1) .+ (2:2:2n)' )
    v = repeat([0.5,1,0.5],n)
    return SparseMatrixCSC(2n+1,n,p,i,v)
end

const n_jac = 1

function relaxed_jacobi_step(x,A,b)
    D = Diagonal(A)
    for i = 1:n_jac
        x += 2/3*(D\(b - A*x))
    end
    return x
end

function coarse_grid_correction(x,A,b)
    @assert all(size(A) .== length(b))
    @assert isodd(length(b))
    n = length(b)÷2  # Integer division: 7÷3 -> 2

    P = pmatrix(n)
    return x + P*(Tridiagonal(P'*A*P)\(P'*(b-A*x)))
end

function twogrid_step(x,A,b)
    x = relaxed_jacobi_step( x, A, b )
    x = coarse_grid_correction( x, A, b )
    return x
end

function coarse_grid_correction_recursive(x,A,b)
    @assert all(size(A) .== length(b))
    @assert isodd(length(b))
    n = length(b)÷2  # Integer division: 7÷2 -> 3

    if n == 1
        return coarse_grid_correction(x,A,b)
    end

    P = pmatrix(n)
    return x + P*multigrid_step(zeros(n), P'*A*P, P'*(b-A*x))
end

function multigrid_step(x,A,b)
    x = relaxed_jacobi_step(x,A,b)
    x = coarse_grid_correction_recursive(x,A,b)
    return x
end

function convergence_history(A,b, step, kmax)
    hist = zeros(kmax)
    x = copy(b)
    for k = 1:kmax
        x = step(x,A,b)
        hist[k] = norm(A*x-b)
    end
    return hist
end

function plot_convergence()
    n = 2^5-1
    kmax = 30
    A = -laplacian(n)
    b = ones(n)

    clf()
    kk = [5,kmax]
    semilogy(kk, 1e1*(1/3^n_jac).^kk, "k--")
    for (name,step) in (
            ("relaxed Jacobi", relaxed_jacobi_step),
            ("two-grid", twogrid_step),
            ("multigrid", multigrid_step),
        )
        hist = convergence_history(A,b, step, kmax)
        semilogy(1:kmax, hist./hist[1], label=name)
    end
    ylim(clamp.(ylim(), 1e-18, Inf))
    legend(loc="best", frameon=false)
    display(gcf())
end
