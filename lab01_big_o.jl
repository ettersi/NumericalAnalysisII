using PyPlot

function fibonacci_sequence()
    n = 30

    # Evaluate
    f = ones(n)
    for i = 3:n
        f[i] = f[i-1] + f[i-2]
    end

    # Plot
    clf()
    semilogy(1:n, f, "-o", label=L"f(n)")
    xlabel(L"n")
    ylabel(L"f(n)")
    legend()
    display(gcf())
    # Conclusion: `f(n)` looks like a straight line on a `semilogy` plot;
    # hence `f(n)` scales exponentially.
end

function triangular_loop()
    n = round.(Int, 10.0.^LinRange(0,4,20))

    # Evaluate
    f = [sum(sum(1 for j = 1:i) for i = 1:n) for n in n]

    # Plot
    clf()
    loglog(n, f, "-o", label=L"f(n)")
    loglog(n, n.^2, "k--", label=L"O(n^2)")
    xlabel(L"n")
    ylabel(L"f(n)")
    legend()
    display(gcf())
    # Conclusion: `f(n)` looks like a straight line on a `loglog` plot with
    # slope parallel to `n^2`. hence `f(n)` scales algebraically with order 2.
end

function geometric_series()
    n = 1:20

    # Evaluate
    f = [sum(2.0.^(1:n)) for n in n]

    # Plot
    clf()
    semilogy(n, f, "-o", label=L"f(n)")
    xlabel(L"n")
    ylabel(L"f(n)")
    legend()
    display(gcf())
    # Conclusion: `f(n)` looks like a straight line on a `semilogy` plot;
    # hence `f(n)` scales exponentially.
end

function recursive()
    n = 2 .^ (1:20)

    # Define function
    f = n -> (if n == 1; return 1; else return 2*f(n÷2)+1; end)

    # Plot
    clf()
    loglog(n, f.(n), "-o", label=L"f(n)")
    loglog(n, n, "k--", label=L"O(n)")
    xlabel(L"n")
    ylabel(L"f(n)")
    legend()
    display(gcf())
    # Conclusion: `f(n)` looks like a straight line on a `loglog` plot with
    # slope parallel to `n`. hence `f(n)` scales algebraically with order 1.
end

function exp_sqrt()
    n = 1:100

    # Evaluate
    f = @. exp(sqrt(n))

    # Plot
    clf()

    subplot(1,2,1)
    loglog(n, f, "-o", label=L"f(n)")
    xlabel(L"n")
    ylabel(L"f(n)")

    subplot(1,2,2)
    semilogy(n, f, "-o", label=L"f(n)")
    xlabel(L"n")
    ylabel(L"f(n)")
    legend()
    display(gcf())
    # Conclusion: `f(n)` is curved upwards in the `loglog` plot but curved
    # downwards in the `semilogy` plot; hence `f(n)` scales super-algebraically
    # but sub-exponentially.
end

function recursive_2()
    n = 2 .^ (1:20)

    # Define function
    f = n -> (if n == 1; return 1; else return f(n÷2)+1; end)

    # Plot
    clf()
    semilogx(n, f.(n), "-o", label=L"f(n)")
    xlabel(L"n")
    ylabel(L"f(n)")
    legend()
    display(gcf())
    # Conclusion: `f(n)` looks like a straight line on a `semilogx` plot;
    # hence it scales logarithmically.
end