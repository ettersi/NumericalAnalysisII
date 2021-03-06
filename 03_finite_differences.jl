using PyPlot
using LinearAlgebra
using SparseArrays

function laplacian_1d(n)
    # Use `Matrix(laplacian_1d(n))` to print the output of this function in a
    # readable form
    return (n+1)^2 * spdiagm(
        -1 => fill( 1.0,n-1), # subdiagonal
         0 => fill(-2.0,n),   # diagonal
         1 => fill( 1.0,n-1)  # superdiagonal
    )
end

function solve_poisson_1d(f)
    n = length(f)
    Δ = laplacian_1d(n)
    return -Δ\f
end

function example_1d()
    n = 4
    x = LinRange(0,1,n+2)[2:end-1]
    f = x -> π^2*sin(π*x)
    uref = x -> sin(π*x)
    u = solve_poisson_1d(f.(x))

    clf()
    xx = LinRange(0,1,1000)
    plot(xx,uref.(xx), "k-", label="exact solution")
    plot([0;x;1],[0;u;0], label="FD solution")
    legend(frameon=false)
    xlabel(L"x")
    ylabel(L"u(x)")
    display(gcf())
end

function convergence_1d()
    # Define problem and solution
    f = x -> π^2 * sin(π*x)
    u = x -> sin(π*x)

    # Compute errors
    n = 2 .^ (1:15)
    error = [begin
        x = LinRange(0,1,n+2)[2:end-1]
        ũ = solve_poisson_1d(f.(x))
        norm(ũ .- u.(x), 2)/sqrt(n+1)
    end for n in n]

    # Plot
    clf()
    loglog(n, error, label=L"\|u - u_n\|_{2,n}")
    loglog(n, n.^-2, "k--", label=L"O(n^{-2})")
    xlabel(L"n")
    legend(frameon=false)
    display(gcf())
end



function laplacian_2d(n)
    Δ = laplacian_1d(n)
    Id = sparse(I,n,n)     # n x n identity matrix
    return kron(Id,Δ) + kron(Δ,Id)
end

function solve_poisson_2d(f)
    @assert size(f,1) == size(f,2)
    n = size(f,1)
    Δ = laplacian_2d(n)
    return reshape(-Δ\vec(f), (n,n))
end

function example_2d()
    n = 100
    x = LinRange(0,1,n+2)[2:end-1]
    f = (x1,x2)->x1*x2
    u = solve_poisson_2d(f.(x,x'))

    clf()
    imshow(u, extent=(0,1,0,1), origin="bottom left")
    colorbar()
    display(gcf())
end

function convergence_2d()
    # Define problem and solution
    f = (x1,x2) -> 5*π^2 * sin(π*x1) * sin(2π*x2)
    u = (x1,x2) -> sin(π*x1) * sin(2π*x2)

    # Compute errors
    n = 2 .^ (1:9)
    error = [begin
        x = LinRange(0,1,n+2)[2:end-1]
        ũ = solve_poisson_2d(f.(x,x'))
        norm(ũ .- u.(x,x'), 2)/(n+1)
    end for n in n]

    # Plot
    clf()
    loglog(n, error, label=L"\|u - u_n\|_{2,n}")
    loglog(n, 2e0*n.^-2, "k--", label=L"O(n^{-2})")
    xlabel(L"n")
    legend(frameon=false)
    display(gcf())
end



function laplacian_3d(n)
    Δ = laplacian_1d(n)
    Id = sparse(I,n,n)    # n x n identity matrix
    return kron(Id,Id,Δ) + kron(Id,Δ,Id) + kron(Δ,Id,Id)
end

function solve_poisson_3d(f)
    @assert size(f,1) == size(f,2) == size(f,3)
    n = size(f,1)
    Δ = laplacian_3d(n)
    return reshape(-Δ\vec(f), (n,n,n))
end

function convergence_3d()
    # Define problem and solution
    f = (x1,x2,x3) -> 3*π^2 * sin(π*x1) * sin(π*x2) * sin(π*x3)
    u = (x1,x2,x3) -> sin(π*x1) * sin(π*x2) * sin(π*x3)

    # Compute errors
    n = 2 .^ (1:5)
    error = [begin
        x = LinRange(0,1,n+2)[2:end-1]
        ũ = solve_poisson_3d(f.(x,x',reshape(x,(1,1,n))))
        norm(ũ .- u.(x,x',reshape(x,(1,1,n))), 2)/(n+1)^(3/2)
    end for n in n]

    # Plot
    clf()
    loglog(n, error, label=L"\|u - u_n\|_{2,n}")
    loglog(n, 5e-1*n.^-2, "k--", label=L"O(n^{-2})")
    xlabel(L"n")
    legend(frameon=false)
    display(gcf())
end