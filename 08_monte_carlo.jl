using PyPlot
using SpecialFunctions
using StaticArrays
using Printf
using Distributions


###############################################################################
# mnk-game

function mnk_probabilities()
    m,n,k = 3,3,3
    N = 100

    println("Runtimes:")
    time_sum = @elapsed p_sum = mnk_probabilities_via_summation(m,n,k)
    time_mc  = @elapsed p_mc  = mnk_probabilities_via_monte_carlo(m,n,k,N)

    println("    Summation: ", @sprintf("%.6f", time_sum), " sec")
    println("  Monte Carlo: ", @sprintf("%.6f", time_mc), " sec")

    s = Matrix{String}(undef,3,3)
    s[:,1] = (p->@sprintf("%.3f",p)).(p_sum)
    s[:,2] = (p->@sprintf("%.3f",p)).(p_mc)
    s[:,3] = (p->@sprintf("%.3f",p)).(abs.(p_mc .- p_sum))

    println()
    println("Results:")
    println("                   |  sum  |   MC  | error ")
    println("-------------------+-------+-------+-------")
    println("Win for player 1:  | ", s[1,1], " | ", s[1,2], " | ", s[1,3])
    println("Win for player 2:  | ", s[2,1], " | ", s[2,2], " | ", s[2,3])
    println("            Draw:  | ", s[3,1], " | ", s[3,2], " | ", s[3,3])
    println()
end

function mnk_probabilities_via_summation(m,n,k)
    board = zeros(Int,m,n)  # Current board state
    moves = vec([(i,j) for i = 1:m, j = 1:n])  # Remaining unoccupied squares
    player = 1  # Current player

    function recurse(board,moves,player)
        # If we run out of moves, then the game is a draw.
        if isempty(moves)
            return [0.0,0.0,1.0]
        end

        p = zeros(3)
        for idx = 1:length(moves)
            # Play out the current move
            i,j = pop_fast!(moves,idx)
            board[i,j] = player

            # Claim victory or recurse
            if winning_move(board,k,i,j)
                p[player] += 1/(length(moves)+1)
            else
                p .+= recurse(board,moves,mod1(player+1,2))./(length(moves)+1)
            end

            # Reset the board and list of moves
            board[i,j] = 0
            push!(moves,(i,j))
            moves[idx],moves[end] = moves[end],moves[idx]
        end
        return p
    end

    return recurse(board,moves,player)
end

function mnk_probabilities_via_monte_carlo(m,n,k,N)
    score = zeros(3)
    for i = 1:N
        board = zeros(Int,m,n)
        moves = vec([(i,j) for i = 1:m, j = 1:n])
        player = 1
        winner = 3

        # Make random moves until none left
        while !isempty(moves)
            i,j = pop_fast!(moves, rand(1:length(moves)))
            board[i,j] = player
            if winning_move(board,k,i,j)
                winner = player
                break
            end
            player = mod1(player+1,2)
        end
        score[winner] += 1
    end
    return score./N
end

function pop_fast!(vec,idx)
    vec[idx],vec[end] = vec[end],vec[idx]
    return pop!(vec)
end

"""
    count_equals(board,i,j,di,dj)

Find the largest integer `c` such that
```
    board[i+k*di,j+k*di] == board[i,j]  for k = 1:c
```
"""
function count_equals(board,i,j,di,dj)
    c = 1
    while i+c*di in 1:size(board,1) &&
          j+c*dj in 1:size(board,2) &&
          board[i+c*di,j+c*dj] == board[i,j]
        c += 1
    end
    return c - 1
end

"""
    winning_move(board,k,i,j)

Check whether the player who owns `board[i,j]` has won the game by taking
that square.
"""
winning_move(board,k,i,j) =
    count_equals(board,i,j, 0,-1) + count_equals(board,i,j, 0, 1) + 1 >= k ||
    count_equals(board,i,j,-1, 0) + count_equals(board,i,j, 1, 0) + 1 >= k ||
    count_equals(board,i,j,-1,-1) + count_equals(board,i,j, 1, 1) + 1 >= k ||
    count_equals(board,i,j,-1, 1) + count_equals(board,i,j, 1,-1) + 1 >= k

###############################################################################


###############################################################################
# High-dimensional integrals

function integral_via_quadrature()
    d = (2,4,8,16)
    N = round.(Int, 2.0.^LinRange(0,16,17))
    f = x-> exp(-sum(x.^2))
    I = sqrt(π)*erf(1)/2
    ylims = [1e-8,4e0]

    clf()
    title(L"Midpoint quadrature ($p = 2$)")
    for (i,d) in enumerate(d)
        n = round.(Int, N.^(1/d))
        loglog(n.^d, [abs(I^d - midpoint(f,d,n))/I^d for n in n], "C$(i-1)", label="d = $d")

        NN = (1e2,N[end])
        offset = (5e-2,1e-1,2e-1,5e-1)
        loglog(NN, offset[i].*inv.(NN).^(2/d), "C$(i-1)--", label=latexstring("O(N^{-$(2/d)})"));
    end
    xlabel(L"# function evaluations $N$")
    ylabel("Error")
    legend()
    ylim(ylims)
    display(gcf())
end

"""
    midpoint(f,d,n)

Compute the integral of `f` over `[0,1]^d` using the midpoint rule with `n`
quadrature points in each dimension.
"""
midpoint(f,d,n) = midpoint_nested(f,n,ntuple(k->n,d))
function midpoint_nested(f,n,nn)
    q = 0.0
    x = LinRange(0,1,2n+1)[2:2:end-1]
    for i in CartesianIndices(nn)
        q += f((ik->x[ik]).(i.I))
    end
    return q/n^length(nn)
end

function integral_via_monte_carlo()
    d = (2,4,8,16)
    N = round.(Int, 2.0.^LinRange(0,16,101))
    f = x-> exp(-sum(x.^2))
    I = sqrt(π)*erf(1)/2
    ylims = [1e-8,4e0]

    clf()
    title("Monte Carlo sampling")
    for d in d
        loglog(N, [abs(I^d - monte_carlo_integral(f,d,N))/I^d for N in N], label="d = $d")
    end
    NN = (1e2,N[end])
    loglog(NN, 6e0.*sqrt.(inv.(NN)), "k--", label=L"O(N^{-0.5})")
    xlabel(L"# function evaluations $N$")
    ylabel("Error")
    legend()
    ylim(ylims)
    display(gcf())
end

"""
    monte_carlo_integral(f,d,N)

Compute the integral of `f` over `[0,1]^d` using `N` uniformly distributed
samples.
"""
monte_carlo_integral(f,d,N) = monte_carlo_integral(f,Val(d),N)
function monte_carlo_integral(f,::Val{d},N) where {d}
    # The above `Val{d}` and the below `@SVector` are Julia-specific
    # performance optimisation. I recommend you ignore them.
    q = 0.0
    for i = 1:N
        q += f(@SVector rand(d))
    end
    return q/N
end



###############################################################################
# Random number generation

using BenchmarkTools
using Random

function rng_benchmark()
    println("Pseudo-random number generator:")
    rng = MersenneTwister()
    @btime rand($rng)
    println()
    println("True random number generator:")
    rng = RandomDevice()
    @btime rand($rng)
end

function inverse_transform()
    U = rand(10_000)
    X = sqrt.(U)

    clf()
    hist(X; bins = 100, density = true, label="empirical PDF")
    xx = LinRange(0,1,1000)
    plot(xx, 2xx, "k", label="theoretical PDF")
    xlabel(L"x")
    legend()
    display(gcf())
end

function bernoulli()
    p = 0.3
    U = rand(10_000)
    X = U .>= (1-p)

    println("Theoretical P(X = 1): ", p)
    println("  Empirical P(X = 1): ", mean(X))
end

function rejection_sampling()
    # Target distribution
    p = x -> ifelse(0 <= x <= 1, 2x, 0.0)
    # p = x -> ifelse(0 <= x <= 1, 6*x*(1-x), 0.0)
    # p = x -> pdf(Normal(0.5,0.1),x)

    # Proposal distribution
    DQ = Uniform(0,1)
    # DQ = Normal(0.5,0.25)
    q = x -> pdf(DQ,x)

    # Compute `M` such that `p(x) ≤ M*q(x)`
    x = rand(DQ,1_000_000)
    M = maximum(p.(x)./q.(x))

    # Log of proposals. For demonstration only
    Qlog = Float64[]

    # Do the rejection sampling
    function sample()
        while true
            Q = rand(DQ)
            push!(Qlog,Q)
            if rand() <= p(Q)/(M*q(Q))
                return Q
            end
        end
    end
    t = @elapsed X = [sample() for k = 1:10_000]

    println("Runtime: ", round(t, digits=4), " seconds")
    println()
    println("# proposals per until acceptance:")
    println("  Theoretical: ", M)
    println("    Empirical: ", length(Qlog)/length(X))

    clf()

    # Plot empirical PDF
    q̃,x = hist(Qlog; bins = 100, density = true, color="white");
    bar(x[1:end-1],M*q̃, diff(x), align="edge", color="C4", alpha=0.5, label="proposal PDF")
    hist(X; bins = 100, density = true, label="empirical PDF")

    # Plot theoretical PDF
    xx = LinRange(0,1,1000)
    plot(xx, p.(xx), "k", label="theoretical PDF")

    xlabel(L"x")
    legend()
    display(gcf())
end

function importance_sampling()
    Y = rand(10_000)
    E_X = mean(2.0.*Y.^2)

    println("    Exact expectation: ", 2/3)
    println("Estimated expectation: ", E_X)
end
