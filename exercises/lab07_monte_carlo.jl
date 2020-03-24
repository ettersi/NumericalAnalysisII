using Random
using PyPlot

function compute_pi(N)
    # TODO
end

function convergence()
    N = round.(Int,10.0.^LinRange(0,6,51))
    q = compute_pi.(N)
    clf()
    loglog(N,N.^(-1/2), "k--")
    loglog(N, abs.(q.-π))
    xlabel(L"N")
    ylabel("error")
    display(gcf())
end
