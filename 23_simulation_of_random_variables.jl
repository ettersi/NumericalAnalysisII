using Random
using Statistics
using BenchmarkTools
using PyPlot

function rng_benchmarks()
    prng = MersenneTwister()
    trng = RandomDevice()
    println("Pseudo-random number generator:")
    @btime rand($prng)
    println()
    println("True random number generator:")
    @btime rand($trng)
end

function transformation_sampling()
    u = rand(1_000_000)
    x = sqrt.(u)

    clf()
    hist(x; bins= 100, density = true)
    xlabel(L"x")
    ylabel(L"\tilde f(x)")
    display(gcf())
end

function rejection_sampling()
    N = 1_000_000
    N_tries = 0
    f = zeros(N)
    for i = 1:length(f)
        g = rand(); N_tries += 1
        while rand() > g
            g = rand(); N_tries += 1
        end
        f[i] = g
    end

    println("Number of tries until accepted: ", N_tries/N)
    clf()
    hist(f; bins= 100, density = true)
    xlabel(L"x")
    ylabel(L"\tilde f(x)")
    display(gcf())
end

function importance_sampling()
    g = rand(1_000_000)
    E_F = mean(2.0.*g.^2)

    println("Exact expectation: ", 2/3)
    println("Estimated expectation: ", E_F)
end
