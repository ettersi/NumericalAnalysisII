using PyPlot
using LinearAlgebra

laplacian(n) = (n+1)^2*Tridiagonal(
    fill( 1.0,n-1), # subdiagonal
    fill(-2.0,n),   # diagonal
    fill( 1.0,n-1)  # superdiagonal
)

function solve_poisson(f, n)
    x = LinRange(0,1,n+2)[2:end-1]
    Δ = laplacian(n)

    # The . on the following line means to apply f to each element.
    # Example: f.([x1,x2,x3]) == [f(x1),f(x[2]),f(x[3])]
    b = f.(x)

    return x, -Δ\b
end

function example()
    f = x->2
    x,u = solve_poisson(f,1000)

    # "clf" stands for "clear figure".
    # If we removed `clf()`, running `example()` again would
    # add another line rather than replacing the previous one.
    clf()

    # The actual plotting command
    plot(x,u)

    # Remember: in Juno you need to run `gcf()` (get current figure)
    # without `;` at the end to update the plot.
end

function convergence()
    # Define problem and solution
    f = x -> π^2 * sin(π*x)
    u = x -> sin(π*x)

    # Compute errors
    n = 2 .^ (1:15)
    error = [begin
        x,ũ = solve_poisson(f,n)
        norm(ũ .- u.(x), 2)/sqrt(n+1)
    end for n in n]

    # Plot
    clf()
    loglog(n, error)
    loglog(n, @.( π^6/24 * (n+1)^-2 ), "k--")
    # `@.(expression)` adds a `.` to every function call in `expression`.
    # Example: `@.(f(g(x)) -> f.(g.(x))`
end
