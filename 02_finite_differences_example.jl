using PyPlot
using LinearAlgebra

function laplacian(n,a)
    x = [(i+1/2)/(n+1) for i = 0:n]
    ax = a.(x)
    return (n+1)^2*Tridiagonal(
        ax[2:end-1],
        -ax[1:end-1] .- ax[2:end],
        ax[2:end-1]
    )
end

function solve_poisson(a,f, n)
    x = LinRange(0,1,n+2)[2:end-1]
    Δ = laplacian(n,a)
    b = f.(x)
    return x, -Δ\b
end

function example()
    a = x-> x < 0.5 ? 0.1 : 1.0
    f = x-> 1
    x,u = solve_poisson(a,f,1000)

    clf()
    plot(x,u)
end
