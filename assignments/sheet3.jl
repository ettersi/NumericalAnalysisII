using PyPlot
using Roots # Install by typing `] add Roots`

function trapezoidal_step(f,y0,t)
    # TODO
end

function implicit_trapezoidal_step(f,y0,t)
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

function convergence()
    f = y->y^2
    y0 = 1.0
    T = 0.5
    y = t-> y0/(1-y0*t)

    clf()
    n = round.(Int, 10.0.^LinRange(1,3,30))
    for (name,step) in (
        ("explicit", trapezoidal_step),
        # ("implict", implicit_trapezoidal_step),
    )
        error = [begin
            ỹ = integrate(f,y0,T,n, step)
            abs(y(T) - ỹ[end])
        end for n in n]
        loglog(n, error, label=name)
    end
    loglog(n, inv.(n).^2, "k--")
    legend(loc="best")
    xlabel(L"n")
    ylabel(L"|\tilde y(T) - y(T)|")
    display(gcf())
end

R_trapezoidal(z) = xxx  # TODO
R_implicit_trapezoidal(z) = xxx  # TODO

function stability()
    f = y -> -y
    y0 = 1.0
    n = 100
    dt = 2.1

    clf()
    for (i,(name, step, R)) in enumerate((
        ("explicit", trapezoidal_step, R_trapezoidal),
        # ("implicit", implicit_trapezoidal_step, R_implicit_trapezoidal),
    ))
        # Compute the solution
        y = integrate(f,y0,dt*(n-1),n,step)

        # Plot the absolute value of the solution
        semilogy(abs.(y), "C$(i-1)-", label=name)

        # Plot `R(-dt)^k` for reference
        nn = (5,n)
        semilogy(nn, 1e1.*abs(R(-dt)).^nn, "C$(i-1)--")
    end
    legend(loc="best")
    ylims = ylim()
    ylim(clamp.(ylims,1e-18,Inf))
    xlabel("step number")
    ylabel(L"|y(k \, \Delta t)|")
    display(gcf())
end
