using PyPlot
using DifferentialEquations
using Random
using LinearAlgebra
using Statistics
using Printf

using PyCall
FuncAnimation = pyimport("matplotlib.animation").FuncAnimation



################################################################################
# Equations of motion

function kinetic_energy(v)
    # `v` is an `2×n` matrix such that `v[:,i] = (velocity of atom i)`
    n = size(v,2)
    E = 0.0
    for i = 1:n
        E += (v[1,i]^2 + v[2,i]^2)/2
    end
    return E
end

function potential_energy(x,L)
    # `x` is an `2×n` matrix such that `x[:,i] = (position of atom i)`
    # `[0,√3*L] × [0,L]` is the simulation box
    n = size(x,2)
    E = 0.0
    for i = 1:n
        for j = 1:i-1
            # Distance between atoms `i` and `j` under periodic boundary conditions.
            # Ignore at first and simply assume `r2 = norm(x[:,i] - x[:,j])^2`.
            r2 = minimum(
                (x[1,i] - x[1,j] + sqrt(3)*L*s1)^2 +
                (x[2,i] - x[2,j] +         L*s2)^2
                for s1 = (-1,0,1), s2 = (-1,0,1)
            )
            ir2 = inv(r2)
            ir6 = ir2^3
            E += ir6^2 - 2*ir6
        end
    end
    return E
end

energy(x,v,L) = kinetic_energy(v) + potential_energy(x,L)

################################################################################


################################################################################
# Initial conditions

function initial_positions(n,L)
    x = Matrix{Float64}(undef, 2,2n^2)
    o = (L/2-n/2) .* (sqrt(3),1)
    for i = 1:n
        for j = 1:n
            x[:,2n*(i-1) + 2j-1] .= (sqrt(3)*(i-1  ), j-1  ) .+ o
            x[:,2n*(i-1) + 2j  ] .= (sqrt(3)*(i-0.5), j-0.5) .+ o
        end
    end
    return x
end

function initial_velocities(n,E)
    v = randn(2,2n^2)
    v .-= mean(v,dims=2)
    EE = kinetic_energy(v,)
    v .*= sqrt(E/EE)
    return v
end

################################################################################



function map_to_box!(x,L) # Ignore at first
    # Enfore periodic boundary conditions: when an atom leaves the simulation
    # box on one side, make it reappear on the other side.
    n = size(x,2)
    for i = 1:n
        if x[1,i] < 0
            x[1,i] += sqrt(3)*L
        elseif x[1,i] > sqrt(3)*L
            x[1,i] -= sqrt(3)*L
        end
        if x[2,i] < 0
            x[2,i] += L
        elseif x[2,i] > L
            x[2,i] -= L
        end
    end
    return x
end

function simulate()
    # Need to escape from VS Code to make the animation work
    pygui(true)

    # Physical parameters
    n = 3  # Number of atoms per dimension
    L = 5*n  # Sidelength of simulation box
    E = 2n^2 * 1.0  # Initial kinetic energy
             # 0.0 -> solid
             # 1.0 -> liquid
             # 2.0 -> gas

    # Numerical parameters
    step = Heun()  # Runge-Kutta method
                   # `Heun()` is same as `trapezoidal_step()` from Lecture 6
    dt = 0.01  # Step size

    # Assemble ODE problem and solver
    Random.seed!(42)
    x = initial_positions(n,L)
    v = initial_velocities(n,E)

    problem = HamiltonianProblem(
        (x,v,_)->energy(x,v,L),
        x,v, (0,Inf)
    )

    # Ignore at first
    enforce_pbc = DiscreteCallback(
        (u,t,integrator) -> true,
        integrator -> begin
            x = integrator.u.x[1]
            map_to_box!(x,L)
        end
    )

    integrator = init(
        problem,
        step,
        adaptive=false, dt = dt,
        alias_u0 = true,
        save_on = false,
        callback=enforce_pbc
    )

    # Plot initial configuration
    clf()
    x = integrator.u.x[1]
    p, = plot(x[1,:],x[2,:], "o", ms = 10)
    xlim([0, sqrt(3)*L])
    ylim([0, L])
    gca().set_aspect("equal","box")

    tlabel = text(0.05*sqrt(3)*L,0.95*L, "", ha="left", va="top")

    # Start animation
    dt_per_frame = max(0.1,dt)
    frames_per_second = 20
    FuncAnimation(
        gcf(),
        i->begin
            step!(integrator, dt_per_frame)
            p.set_data(integrator.u.x[1])
            tlabel.set_text(@sprintf("t = %.1f", i*dt_per_frame))
            return (p,tlabel)
        end,
        interval = 1000/frames_per_second,
        blit=true,
        init_func=()->(p,tlabel)
    )
end
