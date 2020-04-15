using PyPlot
using Statistics

sin_pdf(x) = π/2*sin(π*x)
quad_pdf(x) = 6*x*(1-x)

rand_sin(n) = [rand_sin() for i = 1:n]
rand_sin() = NaN # TODO

rand_quad(n) = [rand_quad() for i = 1:n]
function rand_quad()
    # TODO
end

function histogram(rand_fun)
    n = 1_000_000
    if rand_fun == rand_sin
        pdf = sin_pdf
    elseif rand_fun == rand_quad
        pdf = quad_pdf
    else
        error("Invalid argument rand_fun = $rand_fun")
    end

    clf()
    hist(rand_fun(n); bins = 100, density = true)
    x = LinRange(0,1,1000)
    plot(x, pdf.(x), "k-")
    xlabel(L"x")
    ylabel(L"f(x)")
    display(gcf())
end

function importance_sampling()
    N = 1_000_000
    F = rand_quad(N)
    G = rand_sin(N)

    println("          Var[F] = ", round(var(F), sigdigits=3))
    println("        Var[F]*M = ", round(var(F)*12/π^2, sigdigits=3))
    println("Var[G*f(G)/g(G)] = ", round(var(G.*quad_pdf.(G)./sin_pdf.(G)), sigdigits=3))
end
