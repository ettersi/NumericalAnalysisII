using PyPlot
using Statistics
using Printf
using BenchmarkTools

sin_pdf(x) = π/2*sin(π*x)
quad_pdf(x) = 6*x*(1-x)

rand_sin(n) = [rand_sin() for i = 1:n]
rand_sin() = acos(1-2rand())/π

rand_quad(n) = [rand_quad() for i = 1:n]
function rand_quad()
    M = 12/π^2
    while true
        Y = rand_sin()
        if rand() <= quad_pdf(Y)/(M*sin_pdf(Y))
            return Y
        end
    end
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
    @printf("       Direct sampling: %.3f\n", mean(X))
    @printf("   Importance sampling: %.3f\n", mean(Y.*quad_pdf.(Y)./sin_pdf.(Y)))
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
    t_dir = @belapsed(rand_quad(), seconds=0.1)
    t_imp = @belapsed((Y = rand_sin(); Y * quad_pdf(Y) / sin_pdf(Y)), seconds=0.1)

    println("Runtime per sample:")
    @printf("       Direct sampling: %2.0f nanoseconds\n", 1e9*t_dir)
    @printf("   Importance sampling: %2.0f nanoseconds\n", 1e9*t_imp)
    println()

    #############
    # Comparison
    println("Comparison metric ( sqrt([runtime per sample]*[variance]) ):")
    @printf("       Direct sampling: %.3e\n", sqrt(t_dir*var_dir))
    @printf("   Importance sampling: %.3e\n", sqrt(t_imp*var_imp))
    println()

    # Justifaction for comparison metric:
    #   [Monte Carlo error]
    #     = sqrt( [variance] / [# samples] )
    #     = sqrt( [variance] / ( [runtime] / [runtime per sample] ) )
    #     = sqrt( [runtime per sample] * [variance] / [runtime] )

    # Example output:
    #   Variance:
    #          Direct sampling: 0.0500
    #      Importance sampling: 0.0535
    #
    #   Runtime per sample:
    #          Direct sampling: 61 nanoseconds
    #      Importance sampling: 40 nanoseconds
    #
    #   Comparison metric ( sqrt([runtime per sample]*[variance]) ):
    #          Direct sampling: 5.511e-05
    #      Importance sampling: 4.601e-05

    # Conclusions:
    #   - Direct sampling has a lower variance but larger runtime per sample than
    #     importance sampling.
    #   - All in all, importance sampling is about 5.5 / 4.6 ≈ 1.2 times more
    #     efficient than direct sampling.
end


# Answer to Q7:
# We know from the lecture that rejection sampling on average requires M proposals
# to generate a single sample of the target distribution. Under the given assumptions,
# the runtime-per-sample ratio between direct and importance sampling is hence exactly
# M = 12/π^2 ≈ 1.2, which is reasonably close to the empirically observed factor 1.5.
