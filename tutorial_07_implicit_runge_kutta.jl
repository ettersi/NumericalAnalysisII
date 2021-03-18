using PyPlot
using LinearAlgebra

function euler_step(f,y0,t)
    return y0 + f(y0)*t
end

function trapezoidal_step(f,y0,t)
    f1 = t*f(y0)
    f2 = t*f(y0+f1)
    return y0 + (f1 + f2)/2
end

function semi_implicit_euler_step(f,y0,t)
    D,g = f.D,f.g
    x1,x2,v1,v2 = y0
    v = sqrt(v1^2 + v2^2)
    new_x1 = # TODO: Your code here
    new_x2 = # TODO: Your code here
    new_v1 = # TODO: Your code here
    new_v2 = # TODO: Your code here
    return [new_x1,new_x2,new_v1,new_v2]
end

function propagate(f,y0,T,n,step)
    y = Vector{typeof(y0)}(undef,n)
    y[1] = y0
    for i = 2:n
        y[i] = step(f,y[i-1],T/(n-1))
    end
    return y
end

function cannonball_f(y,D,g)
    x1,x2,v1,v2 = y
    v = sqrt(v1^2 + v2^2)
    return [ v1, v2, -D*v*v1, -D*v*v2-g ]
end

function cannonball_trajectory()
    # Physical parameters
    D = 5.0
    g = 1.0
    y0 = Float64[0,0,1,2]
    T = 2.0

    # Numerical parameters
    n = 100
    step = euler_step

    # Solve the ODE
    y = propagate(
        y->cannonball_f(y,D,g),
        y0, T, n, step
    )
    x1 = [y[i][1] for i = 1:length(y)]
    x2 = [y[i][2] for i = 1:length(y)]
    v1 = [y[i][3] for i = 1:length(y)]
    v2 = [y[i][4] for i = 1:length(y)]

    # Plot the solution
    clf()
    plot(x1,x2)
    axis("equal")
    xlabel(L"x_1")
    ylabel(L"x_2")
    display(gcf())
end

function stability()
    D = 5.0
    g = 1.0
    y0 = Float64[0,0,0, -0.99*sqrt(g/D)]
    n = 100

    clf()
    for (name, step, dt) in (
        ("Euler", euler_step, NaN), # TODO: replace `NaN` with your value for `dt`
        ("trapezoidal", trapezoidal_step, NaN), # TODO: replace `NaN` with your value for `dt`
        # ("semi-implicit Euler", semi_implicit_euler_step, 1e3),
    )
        y = propagate(
            y->cannonball_f(y,D,g),
            y0,dt*(n-1),n,step
        )
        v2 = [y[i][4] for i = 1:length(y)]
        semilogy(abs.(v2 .+ sqrt(g/D)), label=name)
    end
    legend(loc="best")
    xlabel("step number")
    ylabel(L"|v_2 - v_{F,2}|")
    display(gcf())
end

function convergence()
    # Model parameters
    D = 1.0
    g = 1.0
    y0 = Float64[0,0,1,1]
    T = 2.0

    # Compute reference solution
    y = propagate(
        y->cannonball_f(y,D,g),
        y0,T,10000,trapezoidal_step
    )

    clf()
    n = round.(Int, 10.0.^LinRange(0,3,30))
    for (name,step) in (
        ("Euler", euler_step),
        ("trapezoidal", trapezoidal_step),
        # ("semi-implicit Euler", semi_implicit_euler_step),
    )
        error = [begin
            ỹ = propagate(
                y->cannonball_f(y,D,g),
                y0,T,n, step
            )
            norm(y[end] - ỹ[end])
        end for n in n]
        loglog(n, error, label=name)
    end
    loglog(n, inv.(n), "k--")
    loglog(n, inv.(n).^2, "k-.")
    legend(loc="best")
    xlabel(L"n")
    ylabel(L"\|\tilde y(T) - y(T)\|")
    display(gcf())
end
