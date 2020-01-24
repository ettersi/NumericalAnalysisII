using LinearAlgebra
using Random
using PyPlot

function generate_data(N)
    a = rand(N)
    b = 0.2 .+ 0.6.*a .+ 0.1*randn(N)
    return a,b
end

function solve_least_squares(a,b)
    A = [ones(length(a)) a]
    Q,R = qr(A)
    return R\Matrix(Q)'*b
    # Simpler and faster: return A\b
end

function example()
    Random.seed!(42)
    a,b = generate_data(100)
    c = solve_least_squares(a,b)

    clf()
    x = LinRange(0,1,1000)
    plot(a,b,"ko")
    plot(x, @.(c[1]+c[2]*x))
    ylim([0,1])
end


function modified_gram_schmidt(A)
    N,k = size(A)
    Q = zeros(N,k)
    R = zeros(k,k)

    for l = 1:k
        q̃ = A[:,l]
        for m = 1:k-1
            R[m,l] = Q[:,m]' * q̃
            q̃ = q̃ - Q[:,m] * R[m,l]
        end
        R[l,l] = norm(q̃)
        Q[:,l] = q̃/R[l,l]
    end

    return Q,R
end


using Test

function test(A)
    Q,R = modified_gram_schmidt(A)
    @test Q'*Q ≈ I
    @test all(tril(R,-1) .== 0)
    @test Q*R ≈ A
end

function test()
    Random.seed!(42)
    test(rand(4,1))
    test(rand(4,2))
    test(rand(4,3))
    test(rand(4,4))
    test(rand(10,4))
end
