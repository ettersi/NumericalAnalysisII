using PyPlot

function cannonball_f(y,D,g)
    # TODO
    # Hint: it may be useful to split y into its components using
    #    x1,x2,dx1,dx2 = y
end

function euler_step(f,y0,t)
    return y0 + f(y0)*t
end

function midpoint_step(f,y0,t)
    # TODO
end

function ssprk3_step(f,y0,t)
    # TODO
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
    dx1 = [y[i][3] for i = 1:length(y)]
    dx2 = [y[i][4] for i = 1:length(y)]
    return x1,x2,dx1,dx2
end

function plot_trajectory()
    D = 5.0
    g = 1.0
    y0 = Float64[0,0,1,2]
    T = 2.0

    n = 100
    step = euler_step
    y = integrate(
        y->cannonball_f(y,D,g),
        y0, T, n, step
    )
    x1,x2,_ = split(y)
    clf()
    plot(x1,x2)
    axis("equal")
    xlabel(L"x_1")
    ylabel(L"x_2")
    display(gcf())
end

function convergence()
    f = y->y^2
    y0 = 1.0
    T = 0.5
    y = t-> y0/(1-y0*t)

    clf()
    n = round.(Int, 10.0.^LinRange(0,3,30))
    for (name,step) in (
        ("Euler", euler_step),
        ("midpoint", midpoint_step),
        ("SSPRK3", ssprk3_step),
    )
        error = [begin
            ỹ = integrate(f,y0,T,n, step)
            abs(y(T) - ỹ[end])
        end for n in n]
        loglog(n, error, label=name)
    end
    loglog(n, inv.(n), "k--")
    loglog(n, inv.(n).^2, "k-.")
    loglog(n, inv.(n).^3, "k:")
    legend(loc="best")
    xlabel(L"n")
    ylabel(L"|\tilde y(T) - y(T)|")
    display(gcf())
end
