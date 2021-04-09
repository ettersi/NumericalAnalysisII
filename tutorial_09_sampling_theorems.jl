using PyPlot
using Statistics
using Printf
using BenchmarkTools

sin_pdf(x) = π/2*sin(π*x)
quad_pdf(x) = 6*x*(1-x)

rand_sin(n) = [rand_sin() for i = 1:n]
function rand_sin()
    # TODO: Your code here!
    return NaN
end

rand_quad(n) = [rand_quad() for i = 1:n]
function rand_quad()
    # TODO: Your code here!
    return NaN
end

function plot_pdfs()
    x = LinRange(0,1,1000)
    clf()
    plot(x, quad_pdf.(x), label=L"p(x)")
    plot(x, sin_pdf.(x), label=L"q(x)")
    xlabel("x")
    legend()
    display(gcf())
end

function plot_pdf_ratio()
    x = LinRange(0,1,1000)[2:end-1]
    clf()
    plot(x, quad_pdf.(x)./sin_pdf.(x), label=L"p(x) / q(x)")
    xlabel("x")
    legend()
    display(gcf())
end

function histogram(rand_fun)
    n = 1_000_000
    if rand_fun == rand_sin
        pdf = sin_pdf
    elseif rand_fun == rand_quad
        pdf = quad_pdf
    else
        error("Invalid argument rand_fun = $rand_fun")
    end

    clf()
    hist(rand_fun(n); bins = 100, density = true, label="empirical PDF")
    x = LinRange(0,1,1000)
    plot(x, pdf.(x), "k-", label="theoretical PDF")
    xlabel(L"x")
    legend()
    display(gcf())
end

function monte_carlo()
    N = 1000
    X = rand_quad(N)
    Y = rand_sin(N)

    println("Monte Carlo estimate for E[X]")
    @printf("       Direct sampling: %.3f\n", NaN) # TODO: Your code here!
    @printf("   Importance sampling: %.3f\n", NaN) # TODO: Your code here!
end

function comparison()
    ###########
    # Variance
    N = 1_000_000
    X = rand_quad(N)
    Y = rand_sin(N)
    var_dir = var(X)
    var_imp = var(Y.*quad_pdf.(Y)./sin_pdf.(Y))

    println("Variance:")
    @printf("       Direct sampling: %.4f\n", var_dir)
    @printf("   Importance sampling: %.4f\n", var_imp)
    println()

    ##########
    # Runtime
    t_dir = @elapsed rand_quad()
    t_imp = @elapsed (Y = rand_sin(); Y * quad_pdf(Y) / sin_pdf(Y) )

    println("Runtime per sample:")
    @printf("       Direct sampling: %4.0f nanoseconds\n", 1e9*t_dir)
    @printf("   Importance sampling: %4.0f nanoseconds\n", 1e9*t_imp)
    println()

    #############
    # Comparison

    # TODO: Your code here!
end
