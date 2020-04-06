using PyPlot
using DifferentialEquations
using Random
using Statistics
using Printf

function sir_ode(N,I0,k,T)
    problem = ODEProblem(
        (y,p,t) -> begin
            S,I = y
            return [
                -k*S*I/N,
                +k*S*I/N - I
            ]
       end,
       Float64[N-I0,I0],
       (0.0,T)
    )
    solution = solve(problem, Midpoint())
    t = solution.t
    S = solution[1,:]
    I = solution[2,:]
    return t,S,I
end

function sir_gillespie(N,I0,k,T)
    t = 0.0
    S = N-I0
    I = I0

    t_history = [0.0]
    S_history = [S]
    I_history = [I0]

    while t < T
        if I == 0
            t = Inf
        else
            t += randexp()/(k*S*I/N + I)
            if (k*S*I/N+I) * rand() < k*S*I/N
                S -= 1
                I += 1
            else
                I -= 1
            end
        end

        push!(t_history,t)
        push!(S_history,S)
        push!(I_history,I)
    end
    t_history[end] = T

    return t_history,S_history,I_history
end

function simulate()
    N = 1000
    I0 = round(Int, 0.01*N)
    k = 2
    T = 15.0

    clf()

    t,S,I = sir_ode(N,I0,k,T)
    R = N.-S.-I
    plot(t,S, "C0--")
    plot(t,I, "C1--")
    plot(t,R, "C2--")

    t,S,I = sir_gillespie(N,I0,k,T)
    R = N.-S.-I
    plot(t,S, "C0", label="susceptible")
    plot(t,I, "C1", label="infected")
    plot(t,R, "C2", label="recovered")

    legend(loc="best")
    ylim([0,N])
    xlabel("time")
    display(gcf())
end


function exp_sampling()
    λ = 2
    u = rand(1_000_000)
    x = @. -log(1-u)/λ

    clf()
    xx = LinRange(0,maximum(x),1000)
    hist(x; bins = 100, density = true)
    plot(xx, @.(λ*exp(-λ*xx)), "k-")
    xlabel(L"x")
    ylabel(L"\tilde f(x)")
    display(gcf())
end

function bernoulli_sampling()
    p = 0.3
    u = rand(1_000_000)
    x = Int.(u .< p)

    p̃ = mean(x)
    println("Empirical p: ", p̃)
end


function Rinf_distribution()
    N = 100
    I0 = round(Int, 0.01*N)
    k = 2
    n_samples = 100_000

    Rinf = zeros(Int,n_samples)
    for i = 1:n_samples
        t,S,I = sir_gillespie(N,I0,k,Inf)
        Rinf[i] = N-S[end]
    end

    clf()
    hist(Rinf, 0:N, density = true)
    xlabel(L"R(\infty)")
    ylabel(L"P(R(\infty))")
    display(gcf())
end

function p_epidemic()
    N = 100
    I0 = round(Int, 0.01*N)
    k = 2
    Rthres = 40
    n_samples = 100_000

    Rinf = zeros(n_samples)
    for i = 1:n_samples
        t,S,I = sir_gillespie(N,I0,k,Inf)
        Rinf[i] = N-S[end]
    end

    X = Rinf .>= Rthres
    E_X = mean(X)
    var_X = mean((X .- E_X).^2)
    error = sqrt(var_X / n_samples)

    println("        P(Rinf > Rthres): ", E_X)
    println("Estimated absolute error: ", error)
end
