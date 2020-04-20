using PyPlot
using Statistics


# ########################################
# Uniformly distributed points on the disk
# ########################################

function randdisk_rejection()
    while true
        x = -1.0 .+ 2.0.*rand(2)
        if x[1]^2 + x[2]^2 <= 1
            return x
        end
    end
end

function randdisk_transformation()
    r = sqrt(rand())
    phi = 2*π*rand()
    return [
        r*cos(phi),
        r*sin(phi)
    ]
end

function plot_samples()
    n_samples = 1000

    clf()
    for (randdisk,label) in (
            (randdisk_rejection, "Rejection"),
            (randdisk_transformation, "Transformation"),
        )
        x = zeros(2,n_samples)
        for i = 1:n_samples
            x[:,i] = randdisk()
        end
        plot(x[1,:], x[2,:], "o", ms=3, label=label)
    end
    axis("equal")
    display(gcf())
end



# #####################################################
# Importance sampling for highly concentrated integrals
# #####################################################

normal_pdf(m,s,x) = exp(-0.5*(x-m)^2/s^2)/(sqrt(2π)*s)

function uniform_sampling(f,N)
    F = f.(rand(N))
    E_F = mean(F)
    Var_F = mean((F.-E_F).^2)
    return E_F,Var_F
end

function importance_sampling(f,N,m,s)
    X = m .+ s.*randn(N)
    F = f.(X) .* (0 .<= X .<= 1) ./ normal_pdf.(m,s,X)
    E_F = mean(F)
    Var_F = mean((F.-E_F).^2)
    return E_F,Var_F
end

"""
    plot_histogram(estimator,f,N,Nhist)

Create a histogram of `Nhist` realisations of `estimator` with `N` samples per
realisation.

Parameters:
  - `estimator`: either `uniform_sampling` or `importance_sampling`
  - `f`,`N`: as in `uniform_sampling` and `importance_sampling`.
  - `Nhist`: number of realisations of `estimator` to compute for the histogram.
"""
function plot_histogram(estimator,f,N,Nhist; label="")
    # Perform one long run to get accurate estimates for E[F] and Var[F]
    E_F,Var_F = estimator(f,N*Nhist)

    # Perform many short runs to sample the estimator
    E = [estimator(f,N)[1] for i = 1:Nhist]

    # Estimate the expectation and variance of the Monte Carlo estimator using
    # the expectation E_F and variance Var_F of the underlying random variable
    E_E = E_F
    Var_E = Var_F/N

    # Create the histogram of the Monte Carlo estimates
    hist(E, 20, density=true, label=label)

    # Compare the histogram against the normal distribution predicted by the
    # central limit theorem
    x = E_E .+ sqrt(Var_E) .* LinRange(-3,3,1000)
    plot(x, normal_pdf.(E_E, sqrt(Var_E), x), "k")
end

function sin_integral()
    f = x->π/2*sin(π*x)
    N = 100
    Nhist = 10_000

    clf()
    plot_histogram(uniform_sampling,f,N,Nhist, label="Uniform sampling")
    plot([1,1],0.5.*ylim(),"-",lw=5, label="Exact integral")
    legend(loc="best", frameon=false)
    display(gcf())
end

function concentrated_integral()
    f = x->exp(-(20*(x-0.5))^4) / 0.0906401
    N = 100
    Nhist = 10_000

    clf()
    plot_histogram(uniform_sampling,f,N,Nhist, label="Uniform sampling")
    plot_histogram((f,N)->importance_sampling(f,N,0.5,0.05),f,N,Nhist, label="Importance sampling")
    plot([1,1],0.5.*ylim(),"-",lw=5, label="Exact integral")
    legend(loc="best", frameon=false)
    display(gcf())
end

function plot_concentrated_integrand()
    x = LinRange(0,1,1000)
    f = x->exp(-(20*(x-0.5))^4) / 0.0906401
    clf()
    plot(x, f.(x), label=L"f(x)")
    plot(x, normal_pdf.(0.5,0.05,x), label=L"p_Y(x)")
    xlabel(L"x")
    legend(loc="best")
    display(gcf())
end
