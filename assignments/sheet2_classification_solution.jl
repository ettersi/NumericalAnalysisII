using Random
using PyPlot
using Test

function generate_dataset(n)
    t = [8,4,-6]
    xs = [rand(2) for i = 1:n]
    ls = [t[1]*x[1] + t[2]*x[2] + t[3] > rand() for x in xs]
    return xs,ls
end

h(x,t) = t[1]*x[1] + t[2]*x[2] + t[3]
dh(x,t) = [x[1], x[2], 1.0]

g(xs,ls,t) = sum(
    (1-ls[i])*log(1+exp( h(xs[i],t)))
      +ls[i] *log(1+exp(-h(xs[i],t)))
    for i = 1:length(xs)
)
dg(xs,ls,t) = sum(
    (1-ls[i])*inv(1+exp(-h(xs[i],t)))*dh(xs[i],t)
      -ls[i] *inv(1+exp( h(xs[i],t)))*dh(xs[i],t)
    for i = 1:length(xs)
)

function gradient_descent(g,dg,x,nsteps)
    for i = 1:nsteps
        α = 1.0; dx = dg(x)
        while g(x) < g(x - α*dx)
            α /= 2
        end
        x -= α*dx
    end
    return x
end

function solve_and_plot()
    Random.seed!(42)
    n = 1000
    xs,ls = generate_dataset(n)
    t = gradient_descent(
        t-> g(xs,ls,t),  # Objective function
        t->dg(xs,ls,t),  # Derivative
        zeros(3),        # Initial guess
        100              # Number of steps
    )
    @show t

    clf()
    plot([xs[i][1] for i = 1:n if  ls[i]], [xs[i][2] for i = 1:n if  ls[i]], "C2o", ms=4)
    plot([xs[i][1] for i = 1:n if !ls[i]], [xs[i][2] for i = 1:n if !ls[i]], "C3o", ms=4)
    plot([0,1], [-t[3]/t[2],-(t[1]+t[3])/t[2]], "k", lw=2)
    gca().set_aspect("equal", "box")
    xlim([0,1])
    ylim([0,1])
    display(gcf())
end

function test_gradient_descent()
    g = x->(x<=0) + 2*(x>0.25)
    dg = x->-1.0
    @test gradient_descent(g,dg,0.0,1) == 0.25
end

function test_dg()
    @test dg([[1,2],[3,4]],[0,1], [5,6,7]) ≈ [1,2,1]
    @test dg([[4,3],[2,1]],[1,0], [7,6,5]) ≈ [2,1,1]
end
