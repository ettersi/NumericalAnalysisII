using PyPlot
using LinearAlgebra
using SparseArrays

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

function solve_poisson(f, n)
    x = LinRange(0,1,n+2)[2:end-1]
    Δ = laplacian_2d(n)
    b = vec(f.(x,x'))
    return x,reshape(-Δ\b, (n,n))
end

function example()
    f = (x1,x2)->x1*x2
    x,u = solve_poisson(f,300)

    clf()
    imshow(u, extent=(0,1,0,1), origin="bottom left")
    colorbar()
end

function convergence()
    f = (x1,x2)-> 2*π^2 * sin(π*x1) * sin(π*x2)
    u = (x1,x2)-> sin(π*x1) * sin(π*x2)

    n = 2 .^ (1:8)
    error = [begin
        x,ũ = solve_poisson(f,n)
        norm(ũ .- u.(x,x'), 2)/(n+1)
    end for n in n]

    clf()
    loglog(n, error)
    loglog(n, n.^-2, "k--")
    xlabel(L"n")
    ylabel("error")
end

function runtime_data()
    n = round.(Int, 10.0.^LinRange(1,3,10))
    runtime = [begin
        @show n
        Δ = laplacian_2d(n)
        b = zeros(n^2)
        @elapsed Δ\b
    end for n in n]
    return (n,runtime)
end

function runtime_plot((n,runtime))
    clf()
    loglog(n, runtime, "o-", ms=5)
    nn = [1e2,1e3]
    loglog(nn, 1e-7.*nn.^3, "k--")
    xlabel(L"n")
    ylabel("runtime [sec]")
end
