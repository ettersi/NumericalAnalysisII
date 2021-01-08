# Type `] activate .; instantiate` in the REPL to install these packages
using PyPlot
using Distributions

function run_scenario(p1,p2)
    n = [1]
    while 0 < n[end] < 100
        push!(n,
            n[end]
            + rand(Binomial(n[end], p1)) # newly infected
            - rand(Binomial(n[end], p2)) # recovered
        )
    end
    return n
end

function pandemic()
    p1 = 0.5
    p2 = 0.45
    n_samples = 10
    # n_samples = 10_000

    clf() # "clf" == "clear figure"
    count = 0
    for i = 1:n_samples
        n = run_scenario(p1,p2)
        if i < 20 # Don't plot an excessive number of scenarios
            plot(n, "o-", ms=4)
        end
        count += (n[end] > 0)
    end
    xlabel("Day")
    ylabel("Active cases")
    display(gcf()) # Make figure appear in VSCode. "gcf" == "get current figure"
    println("Estimated probability for pandemic: ", count/n_samples)
end