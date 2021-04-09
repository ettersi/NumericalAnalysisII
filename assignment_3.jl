using PyPlot
using Statistics
using BenchmarkTools


# ##############################################################################
# Uniformly distributed points on the disk

function randdisk_rejection()
    # TODO: Your code here
end

function randdisk_transform()
    r = NaN   # TODO: Your code here
    phi = NaN   # TODO: Your code here
    s,c = sincos(phi)
    #     ^^^^^^ Evaluate both sin and cos at once. This turns out to be slightly
    #            faster than evaluating them individually.
    return ( r*c, r*s )
end

function plot_samples()
    n_samples = 1000

    clf()
    for (randdisk,label) in (
            (randdisk_rejection, "Rejection"),
            (randdisk_transform, "Inverse transform"),
        )
        x = zeros(2,n_samples)
        for i = 1:n_samples
            x[:,i] .= randdisk()
        end
        plot(x[1,:], x[2,:], "o", ms=3, label=label)
    end
    axis("square")
    display(gcf())
end

function performance_shootout()
    println("Runtime randdisk_rejection():")
    @btime randdisk_rejection()
    println("Runtime randdisk_transform():")
    @btime randdisk_transform()
end



# ##############################################################################
# Importance sampling for highly concentrated integrals

normal_pdf(m,s,x) = exp(-0.5*(x-m)^2/s^2)/(sqrt(2π)*s)

function uniform_sampling(f,N)
    # TODO: your code here
    return E,Var
end

function importance_sampling(f,N,m,s)
    # TODO: your code here
    return E,Var
end

function plot_histogram(estimator,f,N; label = "")
    # Perform one long run to get accurate estimates for E[X] and Var[X]
    E_X,Var_X = estimator(f,1_000_000)

    # Estimate the expectation E_E and variance Var_E of the Monte Carlo
    # estimator using the expectation E_X and variance Var_X of the underlying
    # random variable
    E_E = E_X
    Var_E = Var_X/N

    # Plot the estimator PDF predicted by the central limit theorem
    x = E_E .+ sqrt(Var_E) .* LinRange(-3,3,1000)
    plot(x, normal_pdf.(E_E, sqrt(Var_E), x), "k")

    # Generate a large number of Monte Carlo estimates and plot the resulting
    # empirical PDF for comparison
    E = [estimator(f,N)[1] for i = 1:10_000]
    hist(E, bins=20, density=true, label=label)
end

function sin_integral()
    f = x->π/2*sin(π*x)
    N = 100

    clf()
    plot_histogram(uniform_sampling,f,N, label="Uniform sampling")
    plot([1,1], 0.5.*ylim(),"-",lw=5, label="Exact integral")
    legend(loc="best", frameon=false)
    display(gcf())
end

function concentrated_integral()
    f = x->exp(-(20*(x-0.5))^4) / 0.0906401
    N = 100

    clf()
    plot_histogram(uniform_sampling,f,N, label="Uniform sampling")
    if (importance = false)
        plot_histogram((f,N)->importance_sampling(f,N,0.5,0.03), f,N, label="Importance sampling")
    end
    plot([1,1],0.5.*ylim(),"-",lw=5, label="Exact integral")
    legend(loc="best", frameon=false)
    display(gcf())
end

function plot_integrand()
    x = LinRange(0,1,1000)
    f = x->exp(-(20*(x-0.5))^4) / 0.0906401
    clf()
    plot(x, f.(x), label=L"f(x)")
    plot(x, normal_pdf.(0.5,0.03,x), label="Gaussian")
    xlabel(L"x")
    legend(loc="best")
    display(gcf())
end
