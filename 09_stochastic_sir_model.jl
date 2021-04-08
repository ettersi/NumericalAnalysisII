using PyPlot
using DifferentialEquations
using Random
using Statistics
using BenchmarkTools
using StatsBase


################################################################################
# Deterministic SIR

function solve_dSIR(N,I0,a,b,T)
    problem = ODEProblem(
        (y,p,t) -> begin
            S,I = y
            return [
                -a*S*I/N,
                +a*S*I/N - b*I
            ]
       end,
       Float64[N-I0,I0],
       (0.0,T)
    )
    solution = solve(problem, Heun())
    t = solution.t
    S = solution[1,:]
    I = solution[2,:]
    return t,S,I
end

function plot_dSIR()
    N = 1000
    I0 = round(Int, 0.01*N)
    a = 2
    b = 1
    T = 15.0

    clf()

    t,S,I = solve_dSIR(N,I0,a,b,T)
    R = N.-S.-I
    plot(t,S, "C0-", label="susceptible")
    plot(t,I, "C1-", label="infected")
    plot(t,R, "C2-", label="recovered")

    legend()
    ylim([0,N])
    xlabel("time")
    display(gcf())
end



################################################################################
# Stochastic SIR

function solve_sSIR(N,I0,a,b,T)
    # Current state. Only these variables are needed to run the algorithm.
    t = 0.0
    S = N-I0
    I = I0

    # History of past states for postprocessing
    t_history = [0.0]
    S_history = [S]
    I_history = [I0]

    # Play out infection and recovery events until we reach the final time
    while t < T
        # Sample times until the next infection and recovery events
        dt_infect = randexp() / (a/N*S*I)
        dt_recover = randexp() / (b*I)

        # Play out whichever event happens first
        if dt_infect < dt_recover
            t += dt_infect
            S -= 1
            I += 1
        else
            t += dt_recover
            I -= 1
        end

        # Push the current state onto the history
        push!(t_history,t)
        push!(S_history,S)
        push!(I_history,I)
    end

    # Final history entry corresponds to time > T
    # Throw out that entry and instead duplicate the penultimate entry
    t_history[end] = T
    S_history[end] = S_history[end-1]
    I_history[end] = I_history[end-1]

    return t_history,S_history,I_history
end

function plot_sSIR()
    if (regime = false)
        N = 1000
        I0 = round(Int, 0.01*N)
        a = 2
        b = 1
        T = 15.0
        ylims = [0,N]
    else
        Random.seed!(1)
        N = 1_000_000
        I0 = 10
        a = 1
        b = 1.1
        T = 50.0
        ylims = [0,60]
    end

    clf()

    t,S,I = solve_dSIR(N,I0,a,b,T)
    R = N.-S.-I
    plot(t,S, "C0--")
    plot(t,I, "C1--")
    plot(t,R, "C2--")

    t,S,I = solve_sSIR(N,I0,a,b,T)
    R = N.-S.-I
    step(t,S, color="C0", where="post", label="susceptible")
    step(t,I, color="C1", where="post", label="infected")
    step(t,R, color="C2", where="post", label="recovered")

    plot([],[], "k--", label="deterministic solutions")

    legend(frameon=false, bbox_to_anchor=(1.0,0.5), loc="center left")
    xlim([0,T])
    ylim(ylims)
    xlabel("time")

    display(gcf())
end



################################################################################
# Max I distribution

function max_I_distribution()
    # Model parameters
    N = 1_000_000
    I0 = 10
    a = 1
    b = 1.1
    T = Inf

    # Numerical parameters
    n = 10_000
    I_cut = 80
    ylims = [1e-5,1]

    # Play out many epidemics and keep counts of the outcomes
    p_dir = zeros(I_cut)
    for i = 1:n
        t,S,I = solve_sSIR(N,I0,a,b,T)
        max_I = maximum(I)
        if max_I <= I_cut
            p_dir[max_I] += 1/n
        end
    end

    clf()
    plot(1:length(p_dir),p_dir, "C0")
    if (errors = true)
        e = @. 3 * sqrt( p_dir * (1-p_dir) / n )
        fill_between(1:I_cut, p_dir.-e, p_dir.+e, color="C0", alpha=0.5)
    end
    yscale("log")
    xlabel("[max I]")
    ylabel("P([max I])")
    xlim([I0,I_cut])
    ylim(ylims)
    display(gcf())
end


################################################################################
# Importance sampling

function solve_isSIR(N,I0,a,b,p_floor,T)
    # Current state. Only these variables are needed to run the algorithm.
    t = 0.0
    S = N-I0
    I = I0
    p = 1.0

    # History of past states for postprocessing
    t_history = [0.0]
    S_history = [S]
    I_history = [I0]
    p_history = [p]

    # Play out infection and recovery events until we reach the final time
    while t < T
        # Precompute infection and recovery propensities
        A = a/N*S*I
        B = b*I

        # Sample the time until the next event (of either kind)
        dt = randexp() / (A+B)
        t += dt

        # Decide whether the next event should be an infection or recovery.
        # We want to simulate trajectories of varying but not exceedingly small
        # probability. I therefore introduce a probability floor `p_floor`
        # and proceed as follows.

        # Check whether infection is at all possible
        if S > 0
            # If the current probability budget allows it...
            q = 0.5
            if p*min(A/q,B/(1-q))/(A+B) > p_floor
                # ... pick an event type at random and update the
                # probability score using the importance sampling formula
                infect = rand() < q
                p *= ifelse(infect,A,B)/(A+B)
            else
                # Otherwise, follow the natural dynamics so we do not lose any
                # further probability mass
                infect = rand() < A/(A+B)
            end
        else
            infect = false
        end

        # Play out the event
        if infect
            S -= 1
            I += 1
        else
            I -= 1
        end

        # Push the current state onto the history
        push!(t_history,t)
        push!(S_history,S)
        push!(I_history,I)
        push!(p_history,p)
    end

    # Final history entry corresponds to time > T
    # Throw out that entry and instead duplicate the penultimate entry
    t_history[end] = T
    S_history[end] = S_history[end-1]
    I_history[end] = I_history[end-1]
    p_history[end] = p_history[end-1]

    return t_history,S_history,I_history,p_history
end

function plot_isSIR()
    # Model parameters
    N = 1_000_000
    I0 = 10
    a = 1
    b = 1.1
    T = 15.0

    # Numerical parameters
    p_floor = 1e-3

    clf()

    subplot(2,1,1)
    t,S,I = solve_dSIR(N,I0,a,b,T)
    R = N.-S.-I
    plot(t,I, "k--", label=L"deterministic $I$")

    t,S,I = solve_sSIR(N,I0,a,b,T)
    step(t,I, color="C1", where="post", label=L"unbiased $I$")

    t,S,I,p = solve_isSIR(N,I0,a,b,p_floor,T)
    step(t,I, color="C3", where="post", label=L"biased $I$")

    legend(frameon=false, bbox_to_anchor=(1.0,0.5), loc="center left")
    xlim([0,T])
    gca().get_xaxis().set_visible(false)

    subplot(2,1,2, sharex=gca())
    step(t,p, label="probability")
    legend(frameon=false, bbox_to_anchor=(1.0,0.5), loc="center left")
    yscale("log")
    xlabel("time")

    display(gcf())
end

function max_I_distribution_with_bias()
    # Model parameters
    N = 1_000_000
    I0 = 10
    a = 1
    b = 1.1
    T = Inf

    # Numerical parameters
    n = 10_000
    I_cut = 80
    p_floor = 0.1

    # Play out many epidemics and keep counts of the outcomes
    p_dir = zeros(I_cut)
    for i = 1:n
        t,S,I = solve_sSIR(N,I0,a,b,T)
        max_I = maximum(I)
        if max_I <= I_cut
            p_dir[max_I] += 1/n
        end
    end

    # Same as above, but now with importance sampling bias
    p_imp = zeros(I_cut)
    p2_imp = zeros(I_cut)
    for i = 1:n
        t,S,I,p = solve_isSIR(N,I0,a,b,p_floor,T)
        max_I = maximum(I)
        if max_I <= I_cut
            p_imp[max_I] += p[end]/n
            p2_imp[max_I] += p[end]^2/n
        end
    end
    @show sum(psums)/n

    clf()

    plot(1:I_cut),p_dir, "C0")
    e = @. 3 * sqrt( p_dir * (1-p_dir) / n )
    fill_between(1:I_cut, p_dir.-e, p_dir.+e, color="C0", alpha=0.5)

    plot(1:I_cut),p_imp, "C1")
    e = @. 3 * sqrt( (p2_imp - p_imp^2) / n )
    fill_between(1:I_cut), p.-e, p.+e, color="C1", alpha=0.5)
    yscale("log")
    xlabel("[max I]")
    ylabel("P([max I])")
    xlim([I0,I_cut])
    ylim(ylims)
    display(gcf())
end


################################################################################
# Sampling Exponential(p)

function exp_sampling()
    p = 2
    U = rand(10_000)
    X = @. -log(1-U)/p

    clf()
    xx = LinRange(0,maximum(X),1000)
    plot(xx, @.(p*exp(-p*xx)), "k-", label="theoretical PDF")
    hist(X; bins = 100, density = true, label="empirical PDF")
    legend()
    display(gcf())
end

function exp_benchmark()
    println("Runtime log(1-rand())")
    @btime log(1-rand())
    println()
    println("Runtime randexp()")
    @btime randexp()
end