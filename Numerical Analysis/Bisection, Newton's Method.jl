# Find all of the roots of the function f(x) = x^5 − 3x^2 + 1 inside the interval [−2, 2].
# Given that the roots of this function are separated by at least 0.25, implement the bisection
# algorithm to obtain estimates of the roots that are accurate to within 0.1.
# Implement Newton’s method to refine the values of the roots so that the absolute difference
# between successive iterates is at most 10^-12


function f(x)
    x^5-3x^2+1
end

function derf(x)
    5x^4-6x
end

function bisection(f,lower, upper)
    # Check that upper and lower are opposite signs
    if f(upper)*f(lower) > 0
        return;
    end

    mid = (lower + upper)/2
    if (upper-lower) < 0.2
        return mid;
    elseif f(upper) * f(mid) > 0
        return bisection(f, lower, mid);
    elseif f(lower)*f(mid)>0
        return bisection(f, mid, upper);
    end
end

function newton(f, derf, Xk)
    print("From ", Xk,  ", ")
    not_done = true
    while not_done
        Xkplus = Xk - (f(Xk)/derf(Xk))
        if abs(Xkplus-Xk) < 10^-12
            not_done = false
        end
        Xk = Xkplus
    end
    println("the refined root is ", Xk)
    return Xk;

end

roots = []

lower = -2
upper = -2

results = []
using Printf
println("Lower   ", "Upper   ", "Bisection Result")
while upper < 2
    if lower >= 0
        print(" ")
    end
    global upper = lower + .25
    rootest = bisection(f, lower, upper)
    @printf("%.2f", lower)
    print("   ")
    if upper >= 0
        print(" ")
    end
    @printf("%.2f", upper)
    print("   ")
    println(rootest)
    if rootest !== nothing
        push!(roots, rootest)
    end
    global lower = upper


end

rootfinal = []

println()
println("Running Newton’s Method on Root Estimates")
for rootest in roots
    fin = newton(f, derf, rootest)
    push!(rootfinal, fin)
end


println()







    


