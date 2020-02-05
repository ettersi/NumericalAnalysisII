using LinearAlgebra
using Random
using Test

function lanczos(A,b,k)
    N = length(b)
    @assert size(A) == (N,N)

    Q = zeros(N,k+1)
    H = zeros(k+1,k)

    Q[:,1] = b / norm(b)
    for l = 1:k
        q̃ = A*Q[:,l]
        H[l,l] = dot(Q[:,l],q̃)
        if l == 1
            q̃ = q̃ - H[l,l]*Q[:,l]
        else
            q̃ = q̃ - H[l,l]*Q[:,l] - H[l-1,l]*Q[:,l-1]
        end
        H[l+1,l] = norm(q̃)
        if l < k
            H[l,l+1] = H[l+1,l]
        end
        Q[:,l+1] = q̃ / H[l+1,l]
    end
    return Q,H
end

function minres(A,b,k)
    Q,H = lanczos(A,b,k)
    b̂ = zeros(k+1); b̂[1] = norm(b)
    y = H \ b̂
    return Q[:,1:end-1]*y
end

function test()
    Random.seed!(42)
    @testset "minres" begin
        @testset for N = 2:9
            A = rand(N,N)
            A += A' # Make sure A is symmetric
            x = rand(N)
            x̃ = minres(A,A*x,N)
            @test x̃ ≈ x
        end
    end
end


function loss_of_orthogonality()
    N = 10
    A = rand(N,N); A = A'*A
    b = rand(N)
    Q,H = lanczos(A,b,N-1)
    display( round.(Q'*Q, sigdigits=2) )
end
