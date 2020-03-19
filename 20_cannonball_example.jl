# Copy-paste the following into the REPL to ensure that you have all the
# required packages:
#    ] add DifferentialEquations Roots ForwardDiff
using PyPlot
using DifferentialEquations
using Roots
using ForwardDiff

function trajectory(v0,D,g)
    problem = ODEProblem(
        (y,p,t) -> begin        # Definition of ODE
            x1,x2,v1,v2 = y
            v = sqrt(v1^2 + v2^2)
            return [v1,v2, -D*v*v1, -D*v*v2-g]
        end,
        [0.0,0.0,v0[1],v0[2]],  # Initial conditions
        (0.0,Inf)               # Time interval
    )

    callback = ContinuousCallback(
        # Condition: when to intervene in the timestepping
        (y,t,integrator) -> y[2],
        # Affect: what to do if the condition is met
        terminate!
    )

    return solve(
        problem,
        Midpoint(),
        callback=callback
    )
end

function plot_trajectory(x)
    plot(x[1,[1,end]], [0,0], "k-")
    plot(x[1,:],x[2,:])
    axis("equal")
    xlabel(L"x_1")
    ylabel(L"x_2")
end

function example()
    v0 = [1.0,2.0]
    D = 5.0
    g = 1.0

    x = trajectory(v0,D,g)
    clf()
    plot_trajectory(x)
    display(gcf())
end

function optimal_angle(D,g)
    F = θ->begin
        θε = ForwardDiff.Dual(θ,1.0)
        s = trajectory([cos(θε),sin(θε)], D,g)

        x1 = ForwardDiff.value(s[1,end])
        x2 = ForwardDiff.value(s[2,end])
        v1 = ForwardDiff.value(s[3,end])
        v2 = ForwardDiff.value(s[4,end])
        dx1 = ForwardDiff.partials(s[1,end],1)  # dx1/dθ(T)
        dx2 = ForwardDiff.partials(s[2,end],1)  # dx2/dθ(T)

        return dx1 - v1/v2*dx2
    end
    return find_zero(F, π/4)
end

function example2()
    D = 5.0
    g = 1.0

    θopt = optimal_angle(D,g)

    clf()
    for θ in θopt .+ (0.0,-0.2,0.2)
        v0 = [cos(θ),sin(θ)]

        x = trajectory(v0,D,g)
        plot_trajectory(x)
    end
    display(gcf())
end
