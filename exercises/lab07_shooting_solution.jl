using PyPlot
using DifferentialEquations
using Roots

function main()
    # Define exact solution and right-hand side
    u = x->sin(π*x)
    ddu = x->-π^2*sin(π*x)
    f = x -> u(x)^2 - ddu(x)

    # Define the mapping d0 -> ODE solution
    solve_ode = d0 -> begin
        problem = ODEProblem(
            (y,p,x) -> begin
                local u,du,up,dup = y
                # du = du/dt
                # up = du/d(d0)
                # dup = d^2u/(dt*d(d0))

                return [du, -f(x) + u^2, dup, 2u*up]
            end,
            [0.0,d0,0.0,1.0], # Initial conditions
            (0.0,1.0)           # Time span
        )
        return solve(problem, Midpoint())
    end

    # Determine d0
    d0 = find_zero(
        d0 -> begin
            sol = solve_ode(d0)
            return (
                sol[1,end],             # f(x)
                sol[1,end]/sol[3,end]   # f(x)/f'(x)
            )
        end,
        float(π),
        Roots.Newton()
    )

    # Compute solution and plot
    sol = solve_ode(d0)
    x = sol.t
    ũ = sol[1,:]

    clf()
    plot(x,u.(x), "k-");
    plot(x,ũ)
    display(gcf())
end
