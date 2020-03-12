using PyPlot

function euler_step(f,y0,t)
    return y0 + f(y0)*t
end

function midpoint_step(f,y0,t)
    ỹ = y0 + f(y0)*t/2
    return y0 + f(ỹ)*t
end

function integrate(f,y0,T, n,step)
    y = Vector{typeof(y0)}(undef,n)
    y[1] = y0
    for i = 2:n
        y[i] = step(f,y[i-1],T/(n-1))
    end
    return y
end

function integrate_adaptively(f,y0,T, tol,step)
    y1 = integrate(f,y0,T,2,step)
    y2 = integrate(f,y0,T,3,step)
    n = 3
    while maximum(abs.(y1.-y2[1:2:end])) > tol
        n = 2n-1
        y1 = y2
        y2 = integrate(f,y0,T,n,step)
    end
    return y2
end

function integrate_adaptively_2(f,y0,T, tol)
    n_trial = 20
    y1 = integrate(f,y0,T,n_trial,euler_step)
    y2 = integrate(f,y0,T,n_trial,midpoint_step)
    error = maximum(abs.(y1 - y2))
    n_ideal = round(Int, n_trial * (error/tol))  # p = 1 for Euler's method
    return integrate(f,y0,T,n_ideal,euler_step)
end

function plot_solution()
    # Model parameters
    λ = 1.0
    f = y->-λ*y
    y0 = 1.0
    T = 10.0

    # Numerical parameters
    step = euler_step
    n = 100_000

    # Solve
    t = LinRange(0,T,n)
    ỹ = integrate(f,y0,T,n, step)

    # Plot
    clf()
    plot(t, @.(exp(-λ*t)), "k-")
    plot(t,ỹ)
    xlabel(L"t")
    ylabel(L"y(t)")
    display(gcf())
end

function plot_adaptive_solution()
    # Model parameters
    λ = 1.0
    f = y->-λ*y
    y0 = 1.0
    T = 10.0

    # Numerical parameters
    step = euler_step
    tol = 1e-3

    # Solve
    ỹ = integrate_adaptively(f,y0,T, tol,step)
    # ỹ = integrate_adaptively_2(f,y0,T, tol)
    t = LinRange(0,T,length(ỹ))

    # Plot
    clf()
    plot(t, @.(exp(-λ*t)), "k-")
    plot(t,ỹ)
    xlabel(L"t")
    ylabel(L"y(t)")
    display(gcf())
end

function plot_error()
    # Model parameters
    λ = 1.0
    f = y->-λ*y
    y0 = 1.0
    T = 10.0

    # Numerical parameters
    step = euler_step
    tol = 1e-2

    # Solve
    ỹ = integrate_adaptively(f,y0,T, tol,step)
    n = length(ỹ)
    t = LinRange(0,T,length(ỹ))

    # Compute local errors
    local_error = zeros(length(ỹ))
    for i = 2:length(ỹ)
        local_error[i] = abs(ỹ[i] - ỹ[i-1]*exp(-λ*(t[i]-t[i-1])))
    end

    # Plot
    clf()
    semilogy(t,local_error, label="Local error")
    semilogy(t,@.(abs(ỹ - exp(-λ*t))), label="Total error")
    xlabel(L"t")
    legend(loc="best")
    display(gcf())
end


using DifferentialEquations

function plot_de_solution()
    # Model parameters
    λ = 1.0
    f = (y,p,t)->-λ*y  # `p` are model parameter, `t` is the time
    y0 = 1.0
    T = 10.0

    # Numerical parameters
    tol = 1e-3
    step = Midpoint()
    # See https://docs.juliadiffeq.org/stable/solvers/ode_solve/#Full-List-of-Methods-1
    # for a list of available stepping algorithms

    # Solve
    problem = ODEProblem(f, y0, (0.0,T))
    solution = solve(problem, step, abstol=tol)
    t = solution.t
    ỹ = solution.u

    # Plot
    clf()
    plot(t, @.(exp(-λ*t)), "k-")
    plot(t,ỹ)
    xlabel(L"t")
    ylabel(L"y(t)")
    display(gcf())
end

function plot_stepsize()
    # Model parameters
    λ = 1.0
    f = (y,p,t)->-λ*y  # `p` are model parameter, `t` is the time
    y0 = 1.0
    T = 100.0

    # Numerical parameters
    tol = 1e-3
    step = Midpoint()
    # See https://docs.juliadiffeq.org/stable/solvers/ode_solve/#Full-List-of-Methods-1
    # for a list of available stepping algorithms

    # Solve
    problem = ODEProblem(f, y0, (0.0,T))
    solution = solve(problem, step, abstol=tol, reltol=0.0)
    t = solution.t

    # Plot
    clf()
    semilogy(t[2:end], diff(t), "o-", ms=3)
    display(gcf())
end

function adaptive_example()
    T = 200.0
    problem = ODEProblem(
        (y,p,t) -> cos(y)^2,  # f(y)
        -1.56,                # y0
        (0.0,T)               # [0,T]
    )
    solution = solve(
        problem,
        Midpoint()
    )
    t = solution.t
    ỹ = solution.u

    Δt = minimum(diff(t))
    n_fixed = round(Int, T/Δt)
    println("Number of adaptive steps: ", length(t)-1)
    println("Number of fixed-size steps: ", n_fixed)

    # Plot
    clf()
    plot(t, zero.(t), "ko", ms=3)
    plot(t,ỹ, "-")
    xlabel(L"t")
    ylabel(L"y(t)")
    display(gcf())
end
