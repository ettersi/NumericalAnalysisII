function conjugate_gradients(A,b,k)
    N = length(b)
    @assert size(A) == (N,N)

    x = zeros(N); r = b; p = r
    for l = 1:k
        α = (r'*r) / (p'*A*p)
        x = x + α*p
        denom = r'*r
        r = r - α*A*p
        β = (r'*r) / denom
        p = r + β*p
    end
    return x
end

function spd_rand(n)
    A = rand(n,n)
    return A'*A
end

using Test
using Random

function test()
    @testset "conjugate gradients" begin
        @testset for N = 2:6
            Random.seed!(42)
            A = spd_rand(N)
            b = rand(N)
            x = conjugate_gradients(A,b,N)
            @test A*x ≈ b
        end
        @testset "repeated eigenvalues" begin
            A = Diagonal([1,1,2,2])
            b = ones(4)
            x = conjugate_gradients(A,b,2)
            @test A*x ≈ b
        end
    end
end
