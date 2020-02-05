using PyPlot
using LinearAlgebra
using IterativeSolvers

function gmres_convergence()
    N = 100
    κ = 10

    clf()
    for (label, A, line_style) in [
        (L"[1,\kappa]", Diagonal(LinRange(1,κ,N)), "-"),
        (L"[10,10\kappa]", Diagonal(LinRange(10,10*κ,N)), "-."),
        (L"[1,\kappa] \cup \{10^3 \,\kappa\}", Diagonal([LinRange(1,κ,N-1);  1e3*κ]), "-"),
        (L"[1,\kappa] \cup \{10^{-3} \,\kappa\}", Diagonal([LinRange(1,κ,N-1); 1e-3*κ]), "-")
    ]
        _,log = gmres(A,ones(N); log=true, restart = N)
        semilogy(
            1:log.iters, log[:resnorm],
            line_style,
            label=label
        )
    end
    idx = [5,28]
    semilogy(idx, 1e2*((sqrt(κ)-1)/(sqrt(κ)+1)).^idx, "k--")
    xlabel(L"k")
    ylabel(L"\|Ax_k - b\|_2")
    legend(loc="best")
    display(gcf())
end

function restarted_gmres_good()
    N = 100
    κ = 10

    clf()
    for k = (N,5)
        A = Diagonal(LinRange(1,κ,N))
        _,log = gmres(A,ones(N); log=true, restart = k)
        semilogy(
            1:log.iters, log[:resnorm],
            label="restart = $k"
        )
    end
    xlabel(L"k")
    ylabel(L"\|Ax_k - b\|_2")
    legend(loc="best")
    display(gcf())
end

function restarted_gmres_bad()
    clf()
    for k = (10,2,1)
        λ = 200 .+ (1:5)
        A = Diagonal([λ; .-λ])
        _,log = gmres(A,ones(10); log=true, restart = k)
        semilogy(
            1:log.iters, log[:resnorm],
            label="restart = $k"
        )
    end
    xlabel(L"k")
    ylabel(L"\|Ax_k - b\|_2")
    legend(loc="best")
    display(gcf())
end
