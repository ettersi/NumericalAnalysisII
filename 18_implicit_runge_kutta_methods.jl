using PyPlot
using Roots

function explicit_euler_step(f,y0,t)
    return y0 + f(y0)*t
end

function explicit_midpoint_step(f,y0,t)
    ỹ = y0 + f(y0)*t/2
    return y0 + f(ỹ)*t
end

function implicit_euler_step(f,y0,t)
    return find_zero(ỹ->y0 + f(ỹ)*t - ỹ, (0,1), Bisection())
    # `find_zero(f,(a,b), Bisection())` determines the root of `f` in the
    # bracketing interval `(a,b)` using the bisection method.
    # The choice `(a,b) = (0,1)` is problem-specific and will not work in general.
end

function implicit_midpoint_step(f,y0,t)
    ỹ = find_zero(y->y0 + f(y)*t/2 - y, (-1,1), Bisection())
    # `find_zero(f,(a,b), Bisection())` determines the root of `f` in the
    # bracketing interval `(a,b)` using the bisection method.
    # The choice `(a,b) = (0,1)` is problem-specific and will not work in general.
    return y0 + f(ỹ)*t
end

function integrate(f,y0,T,m,step)
    y = Vector{typeof(y0)}(undef,m)
    y[1] = y0
    for i = 2:m
        y[i] = step(f,y[i-1],T/(m-1))
    end
    return y
end

function harmonic_oscillator_solution()
    f = y -> [y[2],-y[1]]
    y0 = [1.0,0.0]
    T = 4π

    n = 4
    step = explicit_euler_step

    t = LinRange(0,T,n)
    ỹ = integrate(f,y0,t[end],n, step)
    x̃ = [ỹ[i][1] for i = 1:n]

    @show size(t)
    @show size(x̃)

    clf()
    tt = LinRange(0,T,1000)
    plot(tt, cos.(tt), "k-")
    plot(t, x̃)
    xlabel(L"t")
    ylabel(L"\tilde y(t)")
    display(gcf())
end

function harmonic_oscillator_divergence()
    f = y -> [y[2],-y[1]]
    y0 = [1.0,0.0]
    n = 30
    T = 10π
    t = LinRange(0,T,n)

    clf()
    for (name,step) in (
        ("Euler", explicit_euler_step),
        ("midpoint", explicit_midpoint_step),
    )
        ỹ = integrate(f,y0,t[end],n, step)
        semilogy(t, [abs(ỹ[i][1]) for i = 1:n], label=name)
    end
    z = im*T/(n-1)
    semilogy(t, abs(1+z).^(0:n-1), "k--")
    semilogy(t, abs(1+z+z^2/2).^(0:n-1), "k--")
    xlabel(L"t")
    ylabel(L"|\tilde y(t)|")
    legend(loc="best")
    display(gcf())
end

using DifferentialEquations

function plot_stepsize()
    # Model parameters
    λ = 1.0
    f = (y,p,t)->-λ*y  # `p` are model parameter, `t` is the time
    y0 = 1.0
    T = 100.0

    # Numerical parameters
    tol = 1e-3
    step = ImplicitEuler()
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

function convergence()
    λ = -1.0
    f = y->λ*y
    y0 = 1.0
    T = 10
    y = t->exp(λ*t)

    clf()
    m = round.(Int, 10.0.^LinRange(0,4,100))
    for (name,step) in (
        ("explicit Euler", explicit_euler_step),
        ("explicit midpoint", explicit_midpoint_step),
        ("implicit Euler", implicit_euler_step),
        ("implicit midpoint", implicit_midpoint_step),
    )
        error = [begin
            ỹ = integrate(f,y0,T,m, step)
            abs(y(T) - ỹ[end])
        end for m in m]
        loglog(m, error, label=name)
    end
    nn = (1e2,1e4)
    loglog(nn, 1e-2.*inv.(nn), "k--")
    loglog(nn, 4e-2.*inv.(nn).^2, "k-.")
    legend(loc="best")
    xlabel(L"n")
    ylabel(L"|\tilde y(T) - y(T)|")
    display(gcf())
end
