using LinearAlgebra
using PyPlot

laplacian(n) = (n+1)^2*Tridiagonal(
    fill( 1.0,n-1), # subdiagonal
    fill(-2.0,n),   # diagonal
    fill( 1.0,n-1)  # superdiagonal
)

function jacobi_step(x0,b)
    n = length(x0)
    x1 = zeros(n)

    x1[1] = 0.5*(b[1]/(n+1)^2 + x0[2])
    for i = 2:n-1
        x1[i] = 0.5*(b[i]/(n+1)^2 + x0[i-1] + x0[i+1])
    end
    x1[n] = 0.5*(b[n]/(n+1)^2 + x0[n-1])

    return x1
end

function gauss_seidel_step(x,b)
    n = length(x)

    x[1] = 0.5*(b[1]/(n+1)^2 + x[2])
    for i = 2:n-1
        x[i] = 0.5*(b[i]/(n+1)^2 + x[i-1] + x[i+1])
    end
    x[n] = 0.5*(b[n]/(n+1)^2 + x[n-1])

    return x
end

function convergence_history(b, step, kmax)
    hist = zeros(kmax)
    A = -laplacian(length(b))
    x = copy(b)
    for k = 1:kmax
        x = step(x,b)
        hist[k] = norm(A*x-b)
    end
    return hist
end

function plot_convergence()
    n = 127
    kmax = 10_000
    b = ones(n)

    clf()
    λmax = cos(π / (n+1))
    semilogy(λmax.^(1:kmax), "k--")
    for (name,step) in (
            ("Jacobi", jacobi_step),
            ("Gauss-Seidel",gauss_seidel_step),
            # ("Multigrid",multigrid_step),
        )
        hist = convergence_history(b, step, kmax)
        semilogy(1:kmax, hist./hist[1], label=name)
    end
    legend(loc="best", frameon=false)
    display(gcf())
end
