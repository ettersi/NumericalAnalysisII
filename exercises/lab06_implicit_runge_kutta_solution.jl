using PyPlot
using LinearAlgebra

struct CannonballF
    D::Float64
    g::Float64
end

function (f::CannonballF)(y)
    D,g = f.D,f.g
    x1,x2,v1,v2 = y
    v = sqrt(v1^2 + v2^2)
    return [v1,v2, -D*v*v1, -D*v*v2-g]
end

function euler_step(f,y0,t)
    return y0 + f(y0)*t
end

function midpoint_step(f,y0,t)
    f0 = f(y0)
    f1 = f(y0+f0*t/2)
    return y0 + f1*t
end

function semi_implicit_euler_step(f,y0,t)
    D,g = f.D,f.g
    x1,x2,v1,v2 = y0
    v = sqrt(v1^2 + v2^2)
    new_x1 = x1 + v1*t
    new_x2 = x2 + v2*t
    new_v1 = v1/(1+D*v*t)
    new_v2 = (v2 - g*t)/(1+D*v*t)
    return [new_x1,new_x2,new_v1,new_v2]
end

function integrate(f,y0,T,n,step)
    y = Vector{typeof(y0)}(undef,n)
    y[1] = y0
    for i = 2:n
        y[i] = step(f,y[i-1],T/(n-1))
    end
    return y
end

function split(y)
    x1 = [y[i][1] for i = 1:length(y)]
    x2 = [y[i][2] for i = 1:length(y)]
    v1 = [y[i][3] for i = 1:length(y)]
    v2 = [y[i][4] for i = 1:length(y)]
    return x1,x2,v1,v2
end

function plot_trajectory()
    D = 5.0
    g = 1.0
    f = CannonballF(D,g)
    y0 = Float64[0,0,1,2]
    T = 2.0

    n = 100
    step = euler_step
    y = integrate( f, y0, T, n, step )
    x1,x2,_ = split(y)

    clf()
    plot(x1,x2)
    axis("equal")
    xlabel(L"x_1")
    ylabel(L"x_2")
    display(gcf())
end

function convergence()
    # Model parameters
    D = 1.0
    g = 1.0
    f = CannonballF(D,g)
    y0 = Float64[0,0,1,1]
    T = 2.0

    # Compute reference solution
    y = integrate(f,y0,T,10000,midpoint_step)

    clf()
    n = round.(Int, 10.0.^LinRange(0,3,30))
    for (name,step) in (
        ("Euler", euler_step),
        ("midpoint", midpoint_step),
        ("semi-implicit Euler", semi_implicit_euler_step),
    )
        error = [begin
            ỹ = integrate(f,y0,T,n, step)
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

function stability()
    D = 1.0
    g = 0.5

    f = CannonballF(D,g)
    y0 = Float64[0,0,0,-0.99*sqrt(g/D)]
    n = 100

    clf()
    for (name, step, dt) in (
        ("Euler", euler_step, 1/sqrt(g*D)),
        ("midpoint", midpoint_step, 1/sqrt(g*D)),
        ("semi-implicit Euler", semi_implicit_euler_step, 1e3),
    )
        y = integrate(f,y0,dt*(n-1),n,step)
        x1,x2,v1,v2 = split(y)
        semilogy(abs.(v2 .+ sqrt(g/D)), label=name)
    end
    legend(loc="best")
    xlabel("step number")
    ylabel(L"|v_2 + g/D|")
    display(gcf())
end
