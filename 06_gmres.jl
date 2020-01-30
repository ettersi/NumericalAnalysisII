using LinearAlgebra

function krylov_vectors(A,b,k)
    N = length(b)
    @assert size(A) == (N,N)

    V = zeros(N,k)
    for l = 1:k
        V[:,l] = b
        b = A*b
    end
    return V
end

function gmres_unstable(A,b,k)
    V = krylov_vectors(A,b,k)
    y = ((A*V) \ b)
    return V * y
end

function normalised_krylov_vectors(A,b,k)
    V = krylov_vectors(A,b,k)
    return V ./ [norm(V[:,l]) for l = 1:k]'
end

function arnoldi(A,b,k)
    N = length(b)
    @assert size(A) == (N,N)

    Q = zeros(N,k+1)
    H = zeros(k+1,k)

    Q[:,1] = b / norm(b)
    for l = 1:k
        q̃ = A*Q[:,l]
        for m = 1:l
            H[m,l] = dot(Q[:,m],q̃)
            q̃ = q̃ - Q[:,m]*H[m,l]
        end
        H[l+1,l] = norm(q̃)
        Q[:,l+1] = q̃ / H[l+1,l]
    end
    return Q,H
end

function gmres_slow(A,b,k)
    Q,H = arnoldi(A,b,k)
    y = (A*Q[:,1:end-1]) \ b
    return Q[:,1:end-1]*y
end

function gmres(A,b,k)
    Q,H = arnoldi(A,b,k)
    b̂ = zeros(k+1); b̂[1] = norm(b)
    y = H \ b̂
    return Q[:,1:end-1]*y
end



using Test, Random

function test(gmres_implementation)
    @testset for N = 2:20
        A = rand(N,N)
        x = rand(N)
        x̃ = gmres_implementation(A,A*x,N)
        @test x̃ ≈ x
    end
end

function test()
    Random.seed!(42)
    @testset "gmres_implementations" begin
        @testset "gmres_unstable" begin test(gmres_unstable); end
        @testset "gmres_slow" begin test(gmres_slow); end
        @testset "gmres" begin test(gmres); end
    end
end
