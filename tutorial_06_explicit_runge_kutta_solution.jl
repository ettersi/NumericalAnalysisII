using PyPlot

#############
# Exercise 1

function euler_step(f,y0,t)
    return y0 + f(y0)*t
end

function integrate(f,y0,T,n,step)
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
    y = integrate(
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



#############
# Exercise 2

function midpoint_step(f,y0,t)
    f0 = f(y0)
    f1 = f(y0+f0*t/2)
    return y0 + f1*t
end

function ssprk3_step(f,y0,t)
    f0 = f(y0)
    f1 = f(y0+f0*t)
    f2 = f(y0+f0*t/4+f1*t/4)
    return y0 + f0*t/6 + f1*t/6 + f2*2t/3
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
    xlabel("Number of time steps")
    ylabel(L"Error at final time")
    display(gcf())
end
