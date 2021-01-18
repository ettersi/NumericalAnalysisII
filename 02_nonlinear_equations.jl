function bisection(f, a::Float64,b::Float64) # Ignore the `::Float64` here
    # Check that `[a,b]` is a bracketing interval
    fa = f(a)
    fb = f(b)
    @assert sign(fa) != sign(fb)

    # Ensure `a < b`
    if b < a; a,b = b,a; end

    # Do the bisection
    while b > nextfloat(a)
        # Bisect `[a,b]` such that `[a,m]` and `[m,b]`
        # contain the same number of `Float64`
        m = bisect(a,b)
        fm = f(m)

        # Decide which interval to pursue
        if sign(fa) != sign(fm)
            a,b = a,m
            fa,fb = fa,fm
        else
            a,b = m,b
            fa,fb = fm,fb
        end
    end

    # Return the result. Either `a` or `b` could be returned here.
    return b
end

function bisection_demo()
    # Compute the root of `exp(x) - 2` (i.e. `log(2)`) using the bisection
    # method. In addition, count the number of function evaluations to
    # demonstrate that this number is upper-bounded by roughly 64.
    # (It is not exactly 64 because the `bisect()` function is not perfect, and
    # because we do check that `[a,b]` is a bracketing interval.)

    count = 0  # counts the number of calls to `f(x)`
    f = x -> begin
        count +=1
        return exp(x)-2
    end
    x = bisection(f, -Inf,Inf)
    println("        Error: ", x - log(2))
    println("# evaluations: ", count)
end


#########################################################################
# WARNING: The following code uses Julia and floating-point features not
#          discussed in this module. You are not expected to understand
#          how it works (but you are welcome to come talk to me if you
#          would like to learn more).

# Convert `a` and `b` to `Float64`
# This method is needd to ensure that e.g. `bisection(sin,3,4)` works since `3`
# and `4` are ints, not floats.
bisection(f,a,b) = bisection(f, Float64(a),Float64(b))

"""
    bisect(a::Float64, b::Float64)

Compute `m::Float64` such that `[a,m]` and `[m,b]` contain the same number of
`Float64` (± some small constant).
"""
function bisect(a::Float64,b::Float64)
    a,b = to_int.((a,b))
    m = (a&b) + xor(a,b) >> 1 # Compute `(a+b)÷2` without overflow
    return to_float(m)
end

"""
    to_int(float) -> int

Bijective, monotonous map from `[-Inf,Inf]` to `{-a,...,a}` for some `a > 0`.
Effectively an enumeration of all non-`NaN` `Float64`.
"""
function to_int(float::Float64)
    int = reinterpret(Int64,float)
    if signbit(float)
        int = xor(int, 2^63-1)
    end
    return int
end

"""
    to_float(int)

Inverse of `to_int()`.
"""
function to_float(int::Int64)
    if signbit(int)
        int = xor(int, 2^63-1)
    end
    return reinterpret(Float64,int)
end

#########################################################################



using Printf

function bisection_convergence()
    # Problem parameters
    f = sin
    a,b = 3,4

    for k = 0:12
        # Print the width of the current bracketing interval
        @printf("error(k = %2.d) = %.10f", k, b-a)
        println()
        if mod(k,3) == 0; println(); end

        # Bisect
        m = (b + a)/2
        if sign(f(a)) != sign(f(m))
            a,b = a,m
        else
            a,b = m,b
        end
    end
end

function newton_convergence()
    # Problem parameters
    f = sin
    df = cos
    x = big(1.0)

    for k = 0:6
        # Print the current error
        @printf("error(k = %d) = %.100f", k, abs(x))
        println()

        # Do a Newton step
        x -= f(x) / df(x)
    end
end



using PyPlot

function newton_linear_convergence()
    clf()
    for (i,m) = enumerate(2:4)
        # Problem parameters
        f = x -> x^m
        df = x -> m*x^(m-1)
        x = 1.0
        n = 100

        # Run Newton's method and keep a history of the iterates
        x_hist = zeros(n)
        for k = 1:n
            x_hist[k] = x
            x -= f(x) / df(x)
        end

        # Plot the convergence history
        nn = [50,n]
        r = 1 - inv(m)
        s = 2 * abs(x_hist[nn[1]-1]) *  r^-nn[1]
        semilogy(0:n-1, abs.(x_hist), "C$(i-1)", label=latexstring("f(x) = x^$m"))
        semilogy(nn, s.*r.^nn, "C$(i-1)--", label=latexstring("O(($(m-1)/$m)^{k})"))
    end
    xlabel(L"k")
    ylabel(L"|x_k - x^\star|")
    legend(frameon=false)
    display(gcf())
end



using Printf

function newton_termination()
    # Compute `sqrt(2) = root(x->x^2-2)` using Newton's method and print the
    # error in each iteration.
    # We observe that after six iterations, Newton's method simply jumps back
    # and forth between two values.

    f = x -> x^2 - 2
    df = x -> 2x
    x = 2.0
    for k = 1:15
        @printf("error(k = %2d) = % .2e\n", k, (x - sqrt(big(2))))
        x -= f(x) / df(x)
    end
end



using Roots

function roots_examples()
    @show find_zero(sin, (3.0,4.0), Roots.Bisection())
    @show find_zero((sin,cos), 3.0, Roots.Newton())
    println()
    # @show find_zero(sin, (3.0,4.0), Roots.Bisection(), xatol=0.1)
    # @show find_zero((sin,cos), 4.0, Roots.Newton(), xatol=0.1, atol=1e-3)
    # ^ `atol = 1e-3` is required due to some quirks in how `find_zero()`
    # determines convergence.
end
