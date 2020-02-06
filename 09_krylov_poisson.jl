using LinearAlgebra
using SparseArrays
using IterativeSolvers
using IncompleteLU


# Plotting packages
using PyPlot

ioff() # Disable interactive mode. It messes with small figsizes
rcdefaults()
rc("text", usetex=true)
rc("text.latex", preamble="\\usepackage{amsmath},\\usepackage{amssymb},\\usepackage{sfmath}")
rc("font", size=8)

laplacian_1d(n) = (n+1)^2*Tridiagonal(
    fill( 1.0,n-1), # subdiagonal
    fill(-2.0,n),   # diagonal
    fill( 1.0,n-1)  # superdiagonal
)

function laplacian_2d(n)
    Δ = sparse(laplacian_1d(n))
    Id = sparse(I,n,n)
    return kron(Δ,Id) + kron(Id,Δ)
end


function one_dimensional()
    n = 1000
    A = -laplacian_1d(n)
    b = ones(n)

    fig = figure(figsize=(4.2,1.6))
    subplot(1,2,1)
    for (name,method) in (("MinRes",minres), ("CG",cg))
        _,log = method(A,b, log=true, maxiter=n)
        semilogy(
            1:log.iters,
            log[:resnorm]/sqrt(n),
            label=name
        )
    end
    ylim([1e-5,1e2])
    xlabel(L"k")
    ylabel(L"\|r_k\|_2 / \|b\|_2")
    legend(loc="best", frameon=false)

    subplot(1,2,2)
    n = 10:1000
    niters = [cg(laplacian_1d(n),ones(n), log=true, tol=1e-8)[2].iters for n in n]
    loglog(n, niters, "C2")
    loglog(n, n, "k--", label=L"\mathcal{O}\bigl(n\bigr)")
    xlabel(L"n")
    ylabel(L"k")
    legend(loc="best", frameon=false)

    tight_layout()
    savefig("pics/one_dimensional_laplacian.pdf")
    close(fig)
end

function two_dimensional()
    n = 100
    A = -laplacian_2d(n)
    b = ones(n^2)
    kmax = 3n÷2

    fig = figure(figsize=(4.2,1.6))
    subplot(1,2,1)
    for (name,method) in (("MinRes",minres), ("CG",cg))
        _,log = method(A,b, log=true, maxiter=kmax)
        semilogy(
            1:log.iters,
            log[:resnorm]/n,
            label=name
        )
    end
    κ = 4*(n+1)^2/π^2
    semilogy([1,kmax], 16*((sqrt(κ) - 1) / (sqrt(κ) + 1)).^[1,kmax], "k--")
    xlabel(L"k")
    ylabel(L"\|r_k\|_2 / \|b\|_2")
    legend(loc="lower left", frameon=false)

    subplot(1,2,2)
    n = 10:100
    niters = [cg(laplacian_2d(n),ones(n^2), log=true, tol=1e-8)[2].iters for n in n]
    loglog(n, niters, "C2")
    loglog(n, 2.5n, "k--", label=L"\mathcal{O}\bigl(n\bigr)")
    legend(loc="best", frameon=false)
    xlabel(L"n")
    ylabel(L"k")

    tight_layout()
    savefig("pics/two_dimensional_laplacian.pdf")
    close(fig)
end


function two_dimensional_ilu()
    tol = 10.0.^.-(1:2)

    fig = figure(figsize=(4.2,1.6))
    subplot(1,2,1)
    n = 100
    A = -laplacian_2d(n)
    b = ones(n^2)
    for tol in tol
        _,log = cg(A,b, log=true, Pl = ilu(A; τ = 4*(n+1)^2*tol))
        semilogy(
            1:log.iters,
            log[:resnorm]/n,
            label="tol = $tol"
        )
    end
    xlabel(L"k")
    ylabel(L"\|r_k\|_2 / \|b\|_2")
    xlim([xlim()[1],150])
    legend(loc="best", frameon=false)



    n = 10:10:100
    f = zeros(length(n), length(tol))
    k = zeros(Int, length(n), length(tol))
    for (i,n) = enumerate(n)
        A = -laplacian_2d(n)
        b = ones(n^2)
        for (j,tol) = enumerate(tol)
            F = ilu(A; τ = 4*(n+1)^2*tol)
            _,log = cg(A,b, log=true, tol=1e-8, Pl = F)
            f[i,j] = nnz(F)/nnz(A)
            k[i,j] = log.iters
        end
    end

    subplot(1,2,2)
    for (j,tol) = enumerate(tol)
        loglog(n, k[:,j], "-o", ms=3)
    end
    loglog(n, 2n, "k--", label=L"\mathcal{O}\bigl(n\bigr)")
    legend(loc="best", frameon=false)
    xlabel(L"n")
    ylabel(L"k")

    tight_layout()
    savefig("pics/two_dimensional_laplacian_with_ilu.pdf")
    close(fig)
end
