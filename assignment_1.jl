function newton_inv(y)
    x = 24/17 - (8/17) * y
    for k = 1:TODO   # Your code here!
        x = x * (2 - y*x)
    end
    return x
end

function bisection_inv(y)
    a,b = TODO   # Your code here!
    for k = 1:TODO # Your code here!
        m = (a+b)/2
        if y*m > 1
            a,b = a,m
        else
            a,b = m,b
        end
    end
    return (a+b)/2
end

using Printf
using BenchmarkTools

# Increase the number of benchmarking samples to get more accurate runtimes
BenchmarkTools.DEFAULT_PARAMETERS.samples = 1_000_000

function benchmark()
    @printf("      inv runtime: %6.2f nanoseconds\n", @belapsed(          inv($(Ref(1.0))[]))*1e9)
    @printf("   Newton runtime: %6.2f nanoseconds\n", @belapsed(   newton_inv($(Ref(1.0))[]))*1e9)
    @printf("Bisection runtime: %6.2f nanoseconds\n", @belapsed(bisection_inv($(Ref(1.0))[]))*1e9)
end