using PyPlot
using LinearAlgebra

function laplacian(x)
    n = length(x)
    x = [0;x;1]
    d1 = 1.0./(x[2:end] .- x[1:end-1])
    d2 = 2.0./(x[3:end] .- x[1:end-2])
    return Tridiagonal(
        @.( d2[2:end] * d1[2:end-1] ),
        @.( -d2 * (d1[1:end-1] + d1[2:end]) ),
        @.( d2[1:end-1] * d1[2:end-1] )
    )
end

using Test
function test_laplacian()
    @testset ("laplacian") begin
        @test laplacian([0.1,0.2]) ≈ [
            -200     100
            200/9   -25
        ]
    end
end

function plot_solution()
    # Problem parameters
    n = 10
    p = 2  # Power for grid biasing
    f = x -> 0.25 * x^(-3/2)
    u = x -> sqrt(x) - x

    clf()

    # Plot reference solution
    xx = LinRange(0,1,1000)
    plot(xx, u.(xx), "k-", label="exact solution")

    # Plot finite difference solutions
    for (grid,x) = (
        ("uniform", LinRange(0,1,n+2)[2:end-1]),
        ("adaptive", LinRange(0,1,n+2)[2:end-1].^p),
    )
        Δ = laplacian(x)
        ũ = -Δ\f.(x)
        plot([0;x;1], [0;ũ;0], "-o", label="$grid grid")
    end

    # Add finishing touches to the plot
    xlabel(L"x")
    legend()
    display(gcf())
end

function convergence()
    # Problem parameters
    n = 2 .^ (1:14)
    p = 2  # Power for grid biasing
    f = x -> 0.25 * x^(-3/2)
    u = x -> sqrt(x) - x

    clf()

    # Plot reference lines
    loglog(n, n.^-(1/2), "k:", label=L"O(n^{-1/2})")
    loglog(n, n.^-2, "k--", label=L"O(n^{-2})")

    # Plot the convergence both for uniform and adaptive grids
    for (grid,xfun) = (
        ("uniform", n->LinRange(0,1,n+2)[2:end-1]),
        ("adaptive", n->LinRange(0,1,n+2)[2:end-1].^p),
    )
        errors = [begin
            x = xfun(n)
            Δ = laplacian(x)
            ũ = -Δ\f.(x)
            norm(ũ .- u.(x), Inf)
        end for n in n]
        loglog(n, errors, label="$grid grid")
    end

    # Add finishing touches to the plot
    xlabel(L"n")
    legend(frameon=false)
    display(gcf())
end