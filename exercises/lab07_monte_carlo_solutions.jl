using Random
using PyPlot

function compute_pi(N)
    r = 0
    for i = 1:N
        r += ( rand()^2 + rand()^2 < 1 )
    end
    return 4r/N
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
