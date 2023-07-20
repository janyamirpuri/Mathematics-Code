# Newton's Method for the function f(x) = (x − 0.5)^2
# Starting with the initial guess x0 = 0.25

# Newton's Method for the function f(x) = x^3 − 1. 
# Use the starting value x0 = 0.5.

using Printf

function newton(f, derf, Xk)
    println("# ", "Previous Guess (X_k)       ", "Current Guess (X_k+1)     ", "Difference in Guesses")
    not_done = true
    i = 1
    while not_done
        @printf("%2d", i)
        i = i+1
        change = (f(Xk)/derf(Xk))
        Xkplus = Xk - change
        if Xk >= 0
            print(" ")
        end
        @printf("%.20f", Xk)
        print("   ")
        if Xkplus >= 0
            print(" ")
        end
        @printf("%.20f", Xkplus)
        print("   ")
        
        print(" ")
        
        @printf("%.20f", abs(change))
        print("   ")
        println()


        if abs(Xkplus-Xk) < 10^-12
            not_done = false
        end

        Xk = Xkplus
    end
    println()
    println("The refined root is ", Xk)
    println()

    return Xk;

end


function f(x)
    (x - 0.5)^2
end
println()
function derf(x)
    2*(x - 0.5)
end


newton(f, derf, 0.25)


function f2(x)
    (x^3)-1
end

function derf2(x)
    3*x^2
end

newton(f2, derf2, 0.5)
