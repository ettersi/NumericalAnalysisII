using PyPlot
using LinearAlgebra

function newton(f,df, x, kmax)
    for k = 1:kmax
        x -= df(x)\f(x)
    end
    return x
end

function broyden(f, x1, kmax)
    n = length(x1)
    df_inv = Matrix(I, (n,n))
    for k = 1:kmax
        x1,x0 = x1 - df_inv*f(x1),x1

        dx = x1 - x0
        df = f(x1) - f(x0)
        df_inv += (dx - df_inv*df)*dx'*df_inv / (dx'*df_inv*df)
    end
    return x1
end

function gradient_descent(f, df, x, α, kmax)
    for k = 1:kmax
        x -= α*df(x)'*f(x)
    end
    return x
end

function convergence()
    # Use newton() and broyden() to compute complex square roots.
    w = exp((0.5)*π*im)
    kmax = 20
    α = 0.3
    x0 = [real(w),imag(w)]
    f = x -> [
        x[1]^2 - x[2]^2 - real(w),
            2*x[1]*x[2] - imag(w)
    ]
    df = x-> [
        2x[1]   -2x[2]
        2x[2]    2x[1]
    ]


    f_tracked = x->(push!(xvals,x[1]+x[2]*im); f(x))

    xvals = ComplexF64[]
    newton(f_tracked, df, x0, kmax)
    x_newton = unique(xvals)

    xvals = ComplexF64[]
    broyden(f_tracked, x0, kmax)
    x_broyden = unique(xvals)

    xvals = ComplexF64[]
    gradient_descent(f_tracked, df, x0, α, kmax)
    x_gd = unique(xvals)

    clf()
    semilogy(abs.(x_newton.-sqrt(w)), label="Newton")
    semilogy(abs.(x_broyden.-sqrt(w)), label="Broyden")
    semilogy(abs.(x_gd.-sqrt(w)), label="Gradient descent")
    xlabel(L"k")
    ylabel(L"\|x_k - x^\star\|")
    legend(loc="best")
    display(gcf())
end
