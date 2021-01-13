function matrix_product(A,B)
    @assert size(A,2) == size(B,1)
    # `@assert condition` does nothing if `condition` is `true`, and throws an
    # error otherwise.

    C = zeros(size(A,1),size(B,2))
    for i = 1:size(A,1)
        for j = 1:size(B,2)
            for k = 1:size(A,2)
                C[i,j] += A[i,k]*B[k,j]
            end
        end
    end
    return C
end

function matrix_product()
    n = 1_000
    A = zeros(n,n)
    B = zeros(n,n)
    println(" Naive matrix product: ", @elapsed(matrix_product(A,B)), " seconds")
    println("Export matrix product: ", @elapsed(A*B), " seconds")
end



function sum_ij(A)
    s = 0.0
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            s += A[i,j]
        end
    end
    return s
end

function sum_ji(A)
    s = 0.0
    for j = 1:size(A,2)
        for i = 1:size(A,1)
            s += A[i,j]
        end
    end
    return s
end

function matrix_sum()
    n = 10_000
    A = zeros(n,n)
    println("Summing over columns first: ", @elapsed(sum_ji(A)), " seconds")
    println("   Summing over rows first: ", @elapsed(sum_ij(A)), " seconds")
end



function sum_if(x,y)
    s = 0.0
    for i = 1:length(x)
        if x[i]
            s += y[i]
        end
    end
    return s
end

function branch_prediction()
    # Based on https://stackoverflow.com/q/11227809

    n = 10_000_000
    x_rand = rand(Bool,n)
    x_sort = sort(x_rand)
    y = rand(n)

    println("Conditional sum with random input: ", @elapsed(sum_if(x_rand,y)), " seconds")
    println("Conditional sum with sorted input: ", @elapsed(sum_if(x_sort,y)), " seconds")
    # `sum_if(x,y)` is much faster if `x` is sorted because branch prediction
    # is easier if `x` is of the form `x = [0,0, ..., 0,0,1,1, ... 1,1]`
end



using BenchmarkTools
# Defines the `@belapsed()` macro which estimates the runtime of very short
# operations more accurately by running the operation several times.

function preasymptotic()
    x = zeros(32)
    println("Summing ", length(x)," numbers: ", 1e9*@belapsed(sum($x)), " nanoseconds")
    x = zeros(64)
    println("Summing ", length(x)," numbers: ", 1e9*@belapsed(sum($x)), " nanoseconds")
end

#=
Example output:

Summing 32 numbers: 13.683049147442327 nanoseconds
Summing 64 numbers: 16.608040201005025 nanoseconds

Twice the amount of computations in same amount of time!
=#



function asymptotic()
    x = zeros(1024)
    println("Summing ", length(x)," numbers: ", 1e9*@belapsed(sum($x)), " nanoseconds")
    x = zeros(2048)
    println("Summing ", length(x)," numbers: ", 1e9*@belapsed(sum($x)), " nanoseconds")
    x = zeros(4096)
    println("Summing ", length(x)," numbers: ", 1e9*@belapsed(sum($x)), " nanoseconds")
end

#=
Example output:

Summing 1024 numbers: 77.34751037344398 nanoseconds
Summing 2048 numbers: 162.4756258234519 nanoseconds
Summing 4096 numbers: 325.6359649122807 nanoseconds

No
    twice the number of operations => twice the runtime.
=#



using PyPlot

function algebraic_scaling()
    f = x -> (2 + sin(x))*sqrt(x)
    # `x -> ...` defines an anonymous function with a single argument.
    # In Matlab, the equivalent syntax is `@(x) ...`.
    # In Python, the equivalent syntax is `lambda x: ...`.

    clf()
    if true
        # Bad: Plot `f(x)` using linear axes.
        # Linear axes do not allow us to reliably distinguish between, say,
        # `O(x^(1/2))`, `O(x^(1/3))` or even `O(log(x))`.
        x = LinRange(1,100,1000)
        plot(x, f.(x), label=L"f(x)")  # `f.(x)` applies `f` to every element of `x`.
    else
        # Good: Plot `f(x)` using doubly logarithmic axes
        # `f(x)` is then upper-bounded by a straight line whose slope indicates
        # the power `p` in `f(x) = O(x^p)`.
        x = 10.0.^LinRange(0,3,10000)
        loglog(x, f.(x), label=L"f(x)")
        loglog(x, 4e0.*sqrt.(x), "k--", label=L"O(\sqrt{x})")
    end
    xlabel(L"x")  # `L"[maths]"` is like writing `$[maths]` in LaTeX
    legend()
    display(gcf())
end

function exponential_scaling()
    f = x -> (2 + sin(x))*exp(x)

    clf()
    if true
        # Bad: Plot `f(x)` using linear axes.
        # Linear axes do not allow us to reliably distinguish between, say,
        # `O(x^2)` or `O(exp(x))`.
        x = LinRange(1,5,1000)
        plot(x, f.(x), label=L"f(x)")
    else
        # Good: Plot `f(x)` using a linear x-axis and a logarithmic y-axis.
        # `f(x)` is then upper-bounded by a straight line whose slope indicates
        # the base `a` in `f(x) = O(a^x)`.
        x = LinRange(1,100,1000)
        semilogy(x, f.(x), label=L"f(x")
        semilogy(x, 1e2.*exp.(x), "k--", label=L"O(\exp(x))")
    end
    xlabel(L"x")
    legend()
    display(gcf())
end



function machine_precision()
    # Compute Euler's number using the series expansion of `exp(x)`.
    exp_sum = n -> sum(1/factorial(big(k)) for k = 0:n)

    # Plot the convergence
    n = 1:30
    clf()
    semilogy(n, abs.(exp(1) .- exp_sum.(n)), "o-", label=L"\mathrm{error}(n)")
    semilogy(n[[1,end]], eps().*[1,1], "k--", label=L"\mathrm{eps}()}")
    xlabel(L"n")
    legend()
    display(gcf())
    # Note how the error stagnates at roughly eps() == 2e-16.
end