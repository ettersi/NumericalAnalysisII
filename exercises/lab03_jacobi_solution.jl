using PyPlot
using LinearAlgebra

function jacobi_step(x,A,b)
    D = Diagonal(A)
    return D\(b - (A-D)*x)
end

function gauss_seidel_step(x,A,b)
    L = tril(A,-1)
    D = Diagonal(A)
    U = triu(A,1)
    return (D+U)\(b - L*x)
end

function convergence_history(A,b, step, kmax)
    hist = zeros(kmax)
    x = zeros(length(b))
    for k = 1:kmax
        x = step(x,A,b)
        hist[k] = norm(A*x-b)
    end
    return hist
end

function plot_history(A,b)
    kmax = 50

    L = tril(A,-1)
    D = Diagonal(A)
    U = triu(A,1)

    ρjac = maximum(abs.(eigvals(D\(A-D))))
    ρgs = maximum(abs.(eigvals((D+U)\L)))

    clf()
    semilogy(1:kmax, ρjac.^(1:kmax), "k--", label="O($(round(ρjac, sigdigits=3))^k)")
    semilogy(1:kmax, ρgs.^(1:kmax), "k-.", label="O($(round(ρgs, sigdigits=3))^k)")
    for (name,step) in (
            ("Jacobi", jacobi_step),
            ("Gauss-Seidel",gauss_seidel_step)
        )
        hist = convergence_history(A,b, step, kmax)
        semilogy(1:kmax, hist./hist[1], label=name)
    end
    legend(loc="best", frameon=false)
    xlabel(L"k")
    ylabel(L"\|Ax_k - b\|_2")
    ylim(clamp.(ylim(),1e-16,Inf))
    display(gcf())
end

function task2()
    A = [
        1 2
        2 1
    ]
    plot_history(A, ones(2))
end

function task3_A1()
    A = [
        1   0.75  -0.5
        0    1     0.75
        1    0      1
    ]
    plot_history(A, ones(3))
end

function task3_A2()
    A = [
        1   0.75  0.5
        0    1     0.75
        1    0      1
    ]
    plot_history(A, ones(3))
end
