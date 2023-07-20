
using Plots: plot
using Plots: plot!
using Plots: display
using Plots: savefig
using PDFmerger: append_pdf!
using Plots

# Quadratic Approximation of sin(x) based on L2 norm

function quad_approx(x)
    y = ((45*pi^2-540)/(3*pi^3))*((4*x^2/pi^2)-(4*x/pi) +(2/3)) + (2/pi)
    return y
end

function solve_lagrange(x)
    y = 0
    for k in 1:n
        pro = 1
        for j in 1:n
            if j!=k
                pro = pro*(x - xi[j])/(xi[k]-xi[j])
            end
        end
        y = y + yi[k]*pro       
    end
    return y
end


xi = [0, (pi/2), (pi)]
yi = sin.(xi)
n = 3
x = range(0, pi, length=1500)
yvals0 = sin.(x)
yvals1 = quad_approx.(x)
yvals2 = solve_lagrange.(x)

plt = plot(x, yvals0, label = "True Function", legend=:bottom, fmt =:pdf)
plt = plot!(x, yvals1, label = "Quadratic Approximation in L2", fmt=:pdf)
plt = plot!(x, yvals2, label = "Standard Lagrange with n = 3", fmt =:pdf)

plt = Plots.xaxis!("x")
plt = Plots.yaxis!("f(x)")
plt = Plots.title!("Comparison of Different Estimations of sin(x) on [0, pi]")


savefig(plt, "Question 6 Results.pdf")