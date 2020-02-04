using LinearAlgebra
using PyPlot

# To install this package, type `] add IterativeSolvers` on the REPL.
using IterativeSolvers

function repeated_eigenvalues()
    clf()
    for λ in [
        [1,2,3,4,5,6],
        [1,1,3,4,5,6],
    ]
        A = Diagonal(λ)
        b = ones(length(λ))
        _,log = gmres(A,b, log=true)

        semilogy(
            1:log.iters,
            log[:resnorm],
            label="lambda = $λ",
            "o-"
        )
    end
    xlabel(L"k")
    ylabel(L"\|Ax_k - b\|_2")
    ylim([1e-3,3e0])
    legend(loc="best")
    display(gcf())
end

function alternating_eigenvalues()
    clf()
    for λ in [
        [1,-1,2,-2,3,-3,4,-4,5,-5,6,-6],
    ]
        A = Diagonal(λ)
        b = ones(length(λ))
        _,log = gmres(A,b, log=true)

        semilogy(
            1:log.iters,
            log[:resnorm],
            label="lambda = $λ",
            "o-"
        )
    end
    xlabel(L"k")
    ylabel(L"\|Ax_k - b\|_2")
    ylim([3e-1,5e0])
    legend(loc="best")
    display(gcf())
end
