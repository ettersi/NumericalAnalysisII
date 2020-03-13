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
    return find_zero(y->y0 + f(y)*t - y, (0,1))
    # `find_zero(f,(a,b))` determines the root of `f` in the interval `(a,b)`.
    # The choice `(a,b) = (0,1)` is problem-specific and will not work in general.
end

function implicit_midpoint_step(f,y0,t)
    ỹ = find_zero(y->y0 + f(y)*t/2 - y, (-1,1))
    # `find_zero(f,(a,b))` determines the root of `f` in the interval `(a,b)`.
    # The choice `(a,b) = (-1,1)` is problem-specific and will not work in general.
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
