using PyPlot
using LinearAlgebra
using IterativeSolvers
using Random
using SparseArrays
using Printf

function gmres_example()
    # Define linear system
    n = 100
    A = Diagonal(LinRange(0.1,1.0,n))
    b = ones(n)

    # Call GMRES
    x,log = gmres(
        A,b;
        abstol = 1e-8,    # absolute tolerance
        reltol = 1e-2,    # relative tolerance
        maxiter = 1000,   # max number of iterations / max degree
        restart = n,      # ignore for now
        log = true,       # enable the `log` return value
    )

    # Show errors
    println()
    @printf("Number of iterations: %d\n", log.iters)
    println()
    @printf("Absolute residual: %.3e\n", norm(b - A*x))
    @printf("   Absolute error: %.3e\n", norm(A\b - x))
    println()
    @printf("Relative residual: %.3e\n", norm(b - A*x)/norm(b))
    @printf("   Relative error: %.3e\n", norm(A\b - x)/norm(A\b))
end



function convergence_subtleties()
    # Define matrices
    n = 50
    A1 = Diagonal(LinRange(0.1,1.0,n))
    A2 = spdiagm(
        0 => LinRange(0.1,1.0,n),
        1 => fill(0.55, n-1)
    )
    b1 = ones(n)
    b2 = [zeros(n÷2); ones(n÷2)]

    # Compute condition number of eigenvectors
    d1,V1 = eigen(Matrix(A1))
    d2,V2 = eigen(Matrix(A2))
    println();
    @printf("cond(V1,2) = %.1e\n", cond(V1,2))
    @printf("cond(V2,2) = %.1e\n", cond(V2,2))

    # Plot convergence histories and big O reference lines
    clf()

    _,hist = gmres(A1,b1, abstol=1e-6, restart=n, log=true)
    semilogy([norm(b1); hist[:resnorm]], "C0", label=L"A_1 x = b_1")

    nn = (5,24)
    κ = 10; r = (sqrt(κ)-1)/(sqrt(κ)+1)
    semilogy(nn, 2e1.*r.^nn, "C0--", label=L"O(\rho(10)^n)")

    _,hist = gmres(A1,b2, abstol=1e-6, restart=n, log=true)
    semilogy([norm(b2); hist[:resnorm]], "C1", label=L"A_1 x = b_2")

    nn = (1,8)
    κ = 2; r = (sqrt(κ)-1)/(sqrt(κ)+1)
    semilogy(nn, 2e-1.*r.^nn, "C1--", label=L"O(\rho(2)^n)")

    _,hist = gmres(A2,b1, abstol=1e-6, restart=n, log=true)
    semilogy([norm(b1); hist[:resnorm]], "C2", label=L"A_2 x = b1")

    ylim(clamp.(ylim(), 1e-7,Inf))
    xlabel(L"n")
    ylabel("GMRES residual")
    legend()
    display(gcf())
end



function restarted_gmres_good()
    n = 100
    A = Diagonal(LinRange(1,10,n))
    b = ones(n)

    clf()
    for (i,k) = enumerate((2,5,10,100))
        _,log = gmres(A,b; log=true, restart = k)
        semilogy(
            1:log.iters,
            log[:resnorm],
            "C$(i-1)-"
        )
        semilogy(
            1:k:log.iters,
            log[:resnorm][1:k:end],
            "C$(i-1)o",
            ms = 4
        )
        semilogy(
            [NaN], [NaN],
            "C$(i-1)-o",
            label="restart = $k",
            ms = 4
        )
    end
    xlabel(L"k")
    ylabel(L"\|Ax_k - b\|_2")
    legend(loc="best")
    display(gcf())
end

function restarted_gmres_bad()
    n = 100
    A = Diagonal([0.1; LinRange(1,10,n-1)])
    b = ones(n)

    clf()
    for (i,k) = enumerate((2,5,10,100))
        _,log = gmres(A,b; log=true, restart = k)
        semilogy(
            1:log.iters,
            log[:resnorm],
            "C$(i-1)-"
        )
        semilogy(
            1:k:log.iters,
            log[:resnorm][1:k:end],
            "C$(i-1)o",
            ms = 4
        )
        semilogy(
            [NaN], [NaN],
            "C$(i-1)-o",
            label="restart = $k",
            ms = 4
        )
    end
    xlabel(L"k")
    ylabel(L"\|Ax_k - b\|_2")
    legend(loc="best")
    display(gcf())
end



function gmres_vs_minres()
    n = 200
    Random.seed!(42)
    A = rand(n,n)
    A = A+A' + 12*I
    b = rand(n)

    clf()
    for (label, log) = (
        ("GMRES", gmres(A,b, log=true, restart=length(b))[2]),
        ("MinRes", minres(A,b, log=true)[2]),
    )
        semilogy(log[:resnorm], "-o", ms=2, label=label)
    end
    xlabel(L"m")
    ylabel(L"\|A \, p(A) \, b - b\|_2")
    legend(frameon=false)
    display(gcf())
end



function laplacian_1d(n)
    return (n+1)^2 * spdiagm(
        -1 => fill( 1.0,n-1),
         0 => fill(-2.0,n),
         1 => fill( 1.0,n-1)
    )
end

function laplacian_2d(n)
    Δ = laplacian_1d(n)
    Id = sparse(I,n,n)
    return kron(Id,Δ) + kron(Δ,Id)
end

function laplacian_3d(n)
    Δ = laplacian_1d(n)
    Id = sparse(I,n,n)
    return kron(Id,Id,Δ) + kron(Id,Δ,Id) + kron(Δ,Id,Id)
end


function cg_poisson_1d()
    clf()
    for n = (500,1000,1500,2000)
        Random.seed!(42)
        A = -laplacian_1d(n)
        b = rand(n)
        r = cg(A,b, log = true, reltol = eps())[2][:resnorm]
        semilogy(0:length(r)-1, r, label=latexstring("n = $n"))
    end
    xlabel(L"n")
    ylabel("CG residual norm")
    legend()
    display(gcf())
end

function cg_poisson_2d()
    clf()
    for n = (50,100,150,200)
        Random.seed!(42)
        A = -laplacian_2d(n)
        b = rand(n^2)
        r = cg(A,b, log = true, reltol = eps())[2][:resnorm]
        semilogy(0:length(r)-1, r, label=latexstring("n = $n"))
    end
    xlabel(L"n")
    ylabel("CG residual norm")
    legend()
    display(gcf())
end

function cg_poisson_3d()
    clf()
    for n = (10,20,30,40)
        Random.seed!(42)
        A = -laplacian_3d(n)
        b = rand(n^3)
        r = cg(A,b, log = true, reltol = eps())[2][:resnorm]
        κ = 4*(n+1)^2/π^2
        semilogy(0:length(r)-1, r, label=latexstring("n = $n"))
    end
    xlabel(L"n")
    ylabel("CG residual norm")
    legend()
    display(gcf())
end