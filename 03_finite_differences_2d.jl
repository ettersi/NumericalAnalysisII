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
    Id = sparse(I, (n,n))
    return kron(Δ,Id) + kron(Id,Δ)
end

function solve_poisson(f, n)
    Δ = laplacian_2d(n)
    x = LinRange(0,1,n+2)[2:end-1]
    b = vec(f.(x,x'))
    u = -Δ\b
    return reshape(u, (n,n))
end

function example()
    f = (x1,x2)->x1*x2
    u = solve_poisson(f,300)

    clf()
    imshow(u, extent=(0,1,0,1), origin="bottom left")
    colorbar()
end
