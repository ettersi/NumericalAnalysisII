bisection_point(f,a,b) = (b+a)/2
false_position_point(f,a,b) = (a*f(b) - b*f(a)) / (f(b) - f(a))

"""
    bracket(f,a,b,splitting_point)

Run the bracketing method defined by `splitting_point()` until a tightest
possible bracketing interval is found.
"""
function bracket(f,a,b, splitting_point)
    if f(a) == 0; return a,a; end         # B
    if f(b) == 0; return b,b; end         # C
    for i = 1:10_000
        if nextfloat(a) == b; return a,b; end
        m = splitting_point(f,a,b)
        if f(m) == 0; return m,m; end     # A
        if m == a; m = nextfloat(m); end  # E
        if m == b; m = prevfloat(m); end  # F
        if sign(f(a)) == sign(f(m))
            a,b = m,b
        else
            a,b = a,m
        end
    end
    error(
        "Performed an excessive number of iterations.\n"*
        "Final bracketing interval is [$a,$b]."
    )
end

bisection(f,a,b) = bracket(f,a,b, bisection_point)
false_position(f,a,b) = bracket(f,a,b, false_position_point)

macro test_bracketing_method(method,fstring,a,b)
    # The following code uses advanced programming tools so we can get nice
    # error messages. You are not expected to understand how this code works.
    quote
        method = $method
        fstring = $fstring
        f = $(Meta.parse(fstring))
        a,b = $a,$b
        try
            a,b = method(f,a,b)
            if sign(f(a)) != -sign(f(b))
                error(
                    "Result is not a bracketing interval.\n"*
                    "a = $a, b = $b, f(a) = $(f(a)), f(b) = $(f(b))"
                )
            end
            if a != b && nextfloat(a) != b
                error(
                    "Interval is not tight.\n"*
                    "a = $a, b = $b"
                )
            end
        catch e
            if isa(e, ErrorException)
                println(
                    "$method failed for f = $fstring:\n"*
                    e.msg*"\n"
                )
            else
                bt = catch_backtrace()
                showerror(stderr, e, bt)
            end
        end
    end
end

function test_bracketing_methods()
    @test_bracketing_method(bisection, "x->sin(x)", 3.0, 4.0)
    @test_bracketing_method(bisection, "x->x", -1.0, 1.0)   # Fixed by line A
    @test_bracketing_method(bisection, "x->x",  0.0, 1.0)   # Fixed by line B
    @test_bracketing_method(bisection, "x->x", -1.0, 0.0)   # Fixed by line C

    @test_bracketing_method(false_position, "x->sin(x)", 3.0, 4.0)   # Fixed by lines E and F
    # The following test is in principle fixed by lines E and F, but the fix is
    # not practical.
    @test_bracketing_method(false_position, "x->ifelse(x<=0,-1.0,floatmin())", -1.0, 1.0)
end



function newton(f,df,x, kmax)
    for i = 1:kmax
        x -= f(x)/df(x)
    end
    return x
end

function secant(f,x0,x1, kmax)
    for i = 1:kmax
        if f(x1) == f(x0); return x1; end    # Avoid division by 0
        x1,x0 = (x0*f(x1) - x1*f(x0)) / (f(x1) - f(x0)), x1
    end
    return x1
end

using Test
function test_newton_like_methods()
    @test newton(sin,cos, 3.1, 10) ≈ π
    @test secant(sin, 3.1,3.2, 10) ≈ π
end



function divergence()
    xvals = Float64[]
    newton(
        x->(
            push!(xvals,x);
            atan(x)
        ),
        x->1/(x^2+1),
        3.0, 10
    )
    println("Newton: ", round.(unique(xvals), sigdigits=2))

    xvals = Float64[]
    secant(
        x->(
            push!(xvals,x);
            atan(x)
        ),
        3.0, 3.1, 10
    )
    println("Secant: ", round.(unique(xvals), sigdigits=2))
end


using Printf

function newton_convergence()
    xvals = BigFloat[]
    newton(
        x->(
            push!(xvals,x);
            x^2 - 1
        ),
        x->2x,
        big(2.0), 10
    )
    for x in unique(xvals)
        @printf("%.100f", x)
        println()
    end
end

function secant_convergence()
    xvals = BigFloat[]
    secant(
        x->(
            push!(xvals,x);
            x^2 - 1
        ),
        big(2.0), big(2.1), 10
    )
    for x in unique(xvals)
        @printf("%.100f", x)
        println()
    end
end


using PyPlot

function false_position_convergence()
    f = x->x^2 - 2
    df = x->2x
    d2f = x->2
    xx = sqrt(2)
    a,b = 1.4,1.6
    xvals = Float64[]
    false_position(
        x->(
            push!(xvals,x);
            f(x)
        ),
        a,b
    )
    error = abs.(unique(xvals) .- xx)

    clf()
    semilogy(error)
    idx = 2:13
    semilogy(idx, 1e1*((b-xx)*d2f(xx)/df(xx)/2).^(idx), "k--")
    ylabel(L"|a_k - x^\star|")
    xlabel(L"k")
    display(gcf())
end

function newton_convergence_slow()
    xvals = Float64[]
    newton(
        x->(
            push!(xvals,x);
            x^2
        ),
        x->2x,
        1.0, 100
    )
    error = abs.(xvals)

    clf()
    semilogy(1e4*2.0.^.-(1:length(error)), "k--")
    semilogy(error)
    ylabel(L"|x_k - x^\star|")
    xlabel(L"k")
    display(gcf())
end
