using PyPlot
using Statistics

sin_pdf(x) = π/2*sin(π*x)
quad_pdf(x) = 6*x*(1-x)

rand_sin(n) = [rand_sin() for i = 1:n]
rand_sin() = acos(1-2rand())/π

rand_quad(n) = [rand_quad() for i = 1:n]
function rand_quad()
    M = 12/π^2
    while true
        G = rand_sin()
        if rand() < quad_pdf(G)/(M*sin_pdf(G))
            return G
        end
    end
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

    # Example output:
    #             Var[F] = 0.05
    #           Var[F]*M = 0.0608
    #   Var[G*f(G)/g(G)] = 0.0537
    #
    # Recall that the expected squared error of a Monte Carlo estimator for E[X]
    # is Var[X] / N, which shows that the error is smaller the smaller Var[X].
    # Since Var[F] < Var[G*f(G)/g(G)], it may thus seem like direct Monte Carlo
    # is more efficient than importance sampling Monte Carlo, but this
    # comparison ignores that every sample of F on average requires M samples of
    # G due to rejection sampling.
    #
    # To make a fair comparison, we should compare the error of importance
    # sampling Monte Carlo with N samples to that of direct Monte Carlo with N/M
    # samples. Inserting these numbers into the Var[X] / N error estimates leads
    # us to comparing
    #          Var[F] * M / N for direct Monte Carlo, and
    #    Var[G*f(G)/g(G)] / N for importance sampling Monte Carlo.
    # Since Var[F] * M > Var[G*f(G)/g(G)] in our case, we conclude that
    # importance sampling Monte Carlo is actually more efficient than direct
    # sampling.
end
