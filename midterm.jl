using PyPlot
using LinearAlgebra
using SparseArrays
using IterativeSolvers

function laplacian(D)
    n = length(D)
    return (n+1)^2 .* spdiagm(
        -1 => D[2:end-1],
         0 => .-(D[1:end-1] .+ D[2:end]),
         1 => D[2:end-1],
    )
end

function example()
    # Specify the problem parameters
    n = 100
    D = x -> ifelse(x < 0.5, 0.1, 1.0)
    f = x -> 1.0

    # Assemble the linear system
    x = LinRange(0,1,n+2)[2:end-1]
    m = LinRange(0,1,2n+3)[2:2:end-1]
    Δ = laplacian(D.(m))
    u = -Δ\f.(x)

    # Plot the solution
    clf()
    plot([0;x;1], [0;u;0])
    xlabel(L"x")
    ylabel(L"u(x)")
    display(gcf())
end



function n_iterations()
    D = x -> ifelse(x < 0.5, 0.1, 1.0)
    f = x -> 1.0

    n = round.(Int, 10.0.^LinRange(0,4,10))
    error = [begin
        x = LinRange(0,1,n+2)[2:end-1]
        m = LinRange(0,1,2n+3)[2:2:end-1]
        Δ = laplacian(D.(m))

        _,hist = cg(-Δ,f.(x), log=true)
        hist.iters
    end for n in n]

    clf()
    loglog(n, error, "-o", label="# CG iterations")
    loglog(n, 1.5e0.*n, "k--", label=L"O(n)")
    xlabel(L"# grid points $n$")
    # ylabel("# CG iterations")
    legend()
    display(gcf())
end


###############################################################################
# WARNING: The following code uses advanced Julia features. You are not
# expected to understand the details of how this code works, but you are
# welcome to ask me questions in case you would like to know more.

struct FourierPreconditioner{DST,Diag,iDST}
    V::DST
    iD::Diag
    Vt::iDST
end

using FFTW

function FourierPreconditioner(n)
    # Eigenector matrix and its inverse == transpose
    u = zeros(n)
    V = FFTW.plan_r2r(u, FFTW.RODFT00)
    Vt = inv(V)

    # Inverse of eigenvalue matrix
    iD = Diagonal(@.( inv(2*(n+1)^2*(1 - cos(π*(1:n)/(n+1))))))

    # Assemble `FourierPreconditioner` object
    return FourierPreconditioner(V,iD,Vt)
end

function LinearAlgebra.ldiv!(v,P::FourierPreconditioner,u)
    # Extract matrices from `FourierPreconditioner` object
    V = P.V
    iD = P.iD
    Vt = P.Vt

    # Apply preconditioner to `u` and store result in `v`
    return v .= V*(iD*(Vt*u))
end

# End of advanced Julia code
###############################################################################

function fourier_preconditioning()
    fig = figure(figsize=(8,6))

    for (i,ε) in enumerate(10.0.^.-(0:3))
        # Specify the problem parameters
        n = 50
        D = x -> 1-4*(1-ε)*x*(1-x)
        f = x -> x^2 #1.0 #(x-0.5)^2

        # Assemble the linear system
        x = LinRange(0,1,n+2)[2:end-1]
        m = LinRange(0,1,2n+3)[2:2:end-1]
        A = -laplacian(D.(m))
        b = f.(x)

        # Plot convergence histories
        subplot(2,2,i)

        _,hist = cg(A,b, log=true)
        semilogy([norm(b); hist[:resnorm]], label="no preconditioning")

        _,hist = cg(A,b, Pl=FourierPreconditioner(n), log=true)
        semilogy([norm(b); hist[:resnorm]], label="Fourier preconditioning")

        title(latexstring("\\varepsilon = $ε"))
        xlabel("Number of iterations")
        ylabel("CG residual")
        ylim(1e-6,1e2)
        if i == 4; legend(loc="lower right"); end
    end
    tight_layout(pad=3.0)
    display(gcf())
    close(fig)
end
