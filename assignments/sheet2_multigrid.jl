using LinearAlgebra
using SparseArrays
using PyPlot


const n_jacobi_steps = 2
const jacobi_step_length = #= TODO =#       # See line 41
const jacobi_convergence_rate = #= TODO =#  # See line 103


laplacian_1d(n) = (n+1)^2 * sparse(Tridiagonal(
    fill( 1.0,n-1),
    fill(-2.0,n),
    fill( 1.0,n-1)
))

function laplacian_2d(n)
    Δ = sparse(laplacian_1d(n))
    Id = sparse(I,n,n)
    return kron(Δ,Id) + kron(Id,Δ)
end

function pmatrix_1d(n)
    # I recommend you do *not* try to understand the code for this function.
    # Instead, have a look at `Matrix(pmatrix(3))` to see what this function
    # does.
    p = Vector(1:3:3*n+1)
    i = vec( (-1:1) .+ (2:2:2n)' )
    v = repeat([0.5,1,0.5],n)
    return SparseMatrixCSC(2n+1,n,p,i,v)
end

function pmatrix_2d(n)
    P = pmatrix_1d(n)
    return kron(P,P)
end

function relaxed_jacobi_step(x,A,b)
    D = Diagonal(Vector(diag(A)))
    for i = 1:n_jacobi_steps
        x += jacobi_step_length*(D\(b - A*x))
    end
    return x
end

function coarse_grid_correction(x,A,b)
    @assert length(x) == length(x)
    @assert all(size(A) .== length(b))
    @assert isodd(length(b))
    @assert isinteger(sqrt(length(b)))

    n = Int(sqrt(length(b)))÷2
    P = pmatrix_2d(n)
    return x + P*(Matrix(P'*A*P)\(P'*(b-A*x)))
end

function twogrid_step(x,A,b)
    x = relaxed_jacobi_step( x, A, b )
    x = coarse_grid_correction( x, A, b )
    return x
end

function coarse_grid_correction_recursive(x,A,b)
    @assert length(x) == length(x)
    @assert all(size(A) .== length(b))
    @assert isodd(length(b))
    @assert isinteger(sqrt(length(b)))

    n = Int(sqrt(length(b)))÷2

    if n == 1
        return coarse_grid_correction(x,A,b)
    end

    P = pmatrix_2d(n)
    return x + P*multigrid_step(zeros(n^2), P'*A*P, P'*(b-A*x))
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
    kmax = 40
    A = -laplacian_2d(n)
    b = rand(n^2)

    clf()
    kk = [5,35]
    semilogy(kk, 1e1*(jacobi_convergence_rate^n_jacobi_steps).^kk, "k--")
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
