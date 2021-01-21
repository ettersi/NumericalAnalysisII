#######################
# Bisection vs Newton

using BenchmarkTools
using PyPlot
using Roots

f(x) = 4 - 3x + 2x^2 - x^3
df(x) = -3 + 4x - 3x^2

function plot_function()
    x = LinRange(1,2,1000)
    clf()
    plot(x[[1,end]], [0,0], "k-", lw=0.5)
    plot(x, f.(x))
    display(gcf())
end

function find_root()
    x_bisect = find_zero(f, (1.0,2.0), Roots.Bisection())
    x_Newton = find_zero((f,df), 1.0, Roots.Newton())

    println("Bisection result: ", x_bisect)
    println("   Newton result: ", x_bisect)
end



#######################
# Complex square roots

using LinearAlgebra
using Printf
using Test

function square_root(w; print_error=false)
    f = x -> [
        x[1]^2 - x[2]^2 - real(w),
            2*x[1]*x[2] - imag(w)
    ]
    df = x-> [
        2x[1]   -2x[2]
        2x[2]    2x[1]
    ]

    x = [real(w),imag(w)]
    for k = 1:20
        if print_error
            @printf("error(k = %1.d) = %.2e\n", k, norm(f(x)))
        end

        if norm(f(x)) < 10*eps()*norm(x)
            return x[1] + x[2]*im
        end
        x -= df(x) \ f(x)
    end
    error("Newton's method did not converge. Final iterate is x = $(x[1] + x[2]*im).")
end

function test_square_root()
    @testset "square_root" begin
        @test square_root(1.0) == 1.0
        @test square_root(2.0) ≈ sqrt(2.0)
        @test square_root(1.0im) ≈ sqrt(1.0im)
        @test square_root(2.0im) ≈ sqrt(2.0im)
        @test square_root(1.0+1.0im) ≈ sqrt(1.0+1.0im)
        @test_throws Exception square_root(-1.0)
    end
end

function quadratic_convergence()
    square_root(2im, print_error=true)
    #=
    Output:
        error(k = 1) = 4.47e+00
        error(k = 2) = 1.25e+00
        error(k = 3) = 3.12e-01
        error(k = 4) = 1.28e-02  # ≈ 1e-1^2
        error(k = 5) = 2.05e-05  # ≈ 1e-2^2
        error(k = 6) = 5.24e-11  # ≈ 1e-5^2
        error(k = 7) = 0.00e+00
    =#
end
