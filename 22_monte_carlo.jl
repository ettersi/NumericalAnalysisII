using PyPlot
using SpecialFunctions
using StaticArrays
using Statistics


"""
    midpoint(f,d,n)

Compute the integral of `f` over `[0,1]^d` using the midpoint rule with `n`
quadrature points in each dimension.
"""
midpoint(f,d,n) = midpoint_nested(f,n,ntuple(k->n,d))
function midpoint_nested(f,n,nn)
    q = 0.0
    x = LinRange(0,1,2n+1)[2:2:end-1]
    for i in CartesianIndices(nn)
        q += f((ik->x[ik]).(i.I))
    end
    return q/n^length(nn)
end

"""
    monte_carlo(f,d,N)

Compute the integral of `f` over `[0,1]^d` using `N` uniformly distributed
samples.
"""
monte_carlo(f,d,N) = monte_carlo(f,Val(d),N)
function monte_carlo(f,::Val{d},N) where {d}
    # Passing `d` in the form `::Val{d}` allows the compiler to produce code
    # which is optimised for the particular dimension `d` that we are using.
    q = 0.0
    for i = 1:N
        q += f(@SVector rand(d))  #
    end
    return q/N
end

function convergence()
    d = (2,4,7,10,14)
    N = round.(Int, 10.0.^LinRange(0,6,51))
    f = x-> exp(-sum(x.^2))
    I = sqrt(π)*erf(1)/2
    ylims = [1e-8,4e0]

    clf()

    subplot(1,2,1)
    title("Midpoint")
    for (i,d) in enumerate(d)
        n = round.(Int, N.^(1/d))
        loglog(n.^d, [abs(I^d - midpoint(f,d,n))/I^d for n in n], label="d = $d")

        NN = (1e2,N[end])
        offset = (5e-2,1e-1,2e-1,6e-1,1e0)
        loglog(NN, offset[i].*inv.(NN).^(2/d), "k--");
    end
    xlabel("N")
    ylabel("relative error")
    legend(loc="best")
    ylim(ylims)

    subplot(1,2,2)
    title("Monte Carlo")
    for d in d
        loglog(N, [abs(I^d - monte_carlo(f,d,N))/I^d for N in N])
    end
    NN = (1e2,N[end])
    loglog(NN, 6e0.*sqrt.(inv.(NN)), "k--");
    xlabel("N")
    ylabel("relative error")
    ylim(ylims)
    display(gcf())
end
