using Plots: plot
using Plots: plot!
using Plots: display
using Plots: savefig
using PDFmerger: append_pdf!
using Plots

#Lagrange Interpolating Polynomials of Different Forms with Xi as 2n+1 Nodes

function generate_xi(n)
    nodes = zeros(Float64, 1, 2*n + 1)
    for i in 0:2*n
        x_i = -1 + (i/n)
        nodes[1, i+1] = x_i
    end
    return nodes 
end

function f(x)
    y = x^9 -5*x^4 +3*x-1
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

function solve_bary1(x)
    for i in 1:n
        if x == xi[i]
            return yi[i]
            break
        end
    end
    y = 1
    for i in 1:n
        y = y*(x-xi[i])
    end
    sum = 0 
    for k in 1:n
        pro = 1/(x-xi[k])
        for j in 1:n
            if j!=k
                pro = pro*(1)/(xi[k]-xi[j])
            end
        end
        sum = sum + yi[k]*pro       
    end
    y = y * sum
    return y
end
function solve_bary2(x)
    for i in 1:n
        if x == xi[i]
            return yi[i]
            break
        end
    end
    sum1 = 0
    sum2 = 0
    for k in 1:n
        pro = 1/(x-xi[k])
        for j in 1:n
            if j!=k
                pro = pro*(1)/(xi[k]-xi[j])
            end
        end
        sum1 = sum1 + yi[k]*pro
        sum2 = sum2 + pro   
    end
    y = sum1/sum2
    return y
end

xi = generate_xi(4)
yi = f.(xi)
n = 9


xtest = collect(Float64, -1:0.002:1)
yvals0 = f.(xtest)
yvals1 = solve_lagrange.(xtest)
yvals2 = solve_bary1.(xtest)
yvals3 = solve_bary2.(xtest)

plt = plot(xtest, yvals0, label = "True Function", fmt =:pdf)
plt = plot!(xtest, yvals1, label = "Standard Lagrange", fmt=:pdf)
plt = plot!(xtest, yvals2, label = "First Barycentric", fmt =:pdf)
plt = plot!(xtest, yvals3, label = "Second Barycentric", fmt =:pdf)
plt = Plots.xaxis!("x")
plt = Plots.yaxis!("f(x)")
plt = Plots.title!("b) Lagrange Interpolations (n = 4/ 9 nodes)")

plt2 = plot(xtest, log10.(abs.(yvals1-yvals0)), label = "Standard Lagrange Error", fmt =:pdf)
plt2 = plot!(xtest, log10.(abs.(yvals2-yvals0)), label = "First Barycentric Error", fmt =:pdf)
plt2 = plot!(xtest, log10.(abs.(yvals3-yvals0)), label = "Second Barycentric Error", fmt =:pdf)
plt2 = Plots.xaxis!("x")
plt2 = Plots.yaxis!("Log10 Error of Interpolated Value")
plt2 = Plots.title!("b) Errors of Lagrange Interpolations (n = 4/ 9 nodes)")

xi = generate_xi(5)
yi = f.(xi)
n = 11

yvals01 = f.(xtest)
yvals11 = solve_lagrange.(xtest)
yvals21 = solve_bary1.(xtest)
yvals31 = solve_bary2.(xtest)

plt3 = plot(xtest, yvals01, label = "True Function", fmt =:pdf)
plt3 = plot!(xtest, yvals11, label = "Standard Lagrange", fmt=:pdf)
plt3 = plot!(xtest, yvals21, label = "First Barycentric", fmt =:pdf)
plt3 = plot!(xtest, yvals31, label = "Second Barycentric", fmt =:pdf)
plt3 = Plots.xaxis!("x")
plt3 = Plots.yaxis!("f(x)")
plt3 = Plots.title!("c) Lagrange Interpolations (n = 5/ 11 nodes)")

plt4 = plot(xtest, log10.(abs.(yvals11-yvals01)), label = "Standard Lagrange Error", fmt =:pdf)
plt4 = plot!(xtest, log10.(abs.(yvals21-yvals01)), label = "First Barycentric Error", fmt =:pdf)
plt4 = plot!(xtest, log10.(abs.(yvals31-yvals01)), label = "Second Barycentric Error", fmt =:pdf)
plt4 = Plots.xaxis!("x")
plt4 = Plots.yaxis!("Log10 Error of Interpolated Value")
plt4 = Plots.title!("c) Errors of Lagrange Interpolations (n = 5/ 11 nodes)")

for p in [plt, plt2, plt3, plt4]
    savefig(p, "temp.pdf")
    append_pdf!("bary.pdf", "temp.pdf", cleanup=true)
end






# ###NOT IMPORTANT

# function foil(x1, y1) #when multiplying any polynomial by (a + bx) ==> [a , b]
#     x2 = zeros(Float64, 1, length(x1)+1)
#     for i in range(1, length(x1))
#         x2[i] = y1[1]*x1[i]
#     end
#     for i in range(1, length(x1))
#         x2[i+1] = x2[i+1] + y1[2]*x1[i]
#     end
#     if x1 == [1,0]
#         x2 = y1
#     end
#     return x2
# end

# function generate_lagrange_coefs(xi, yi, n) #creates the coefficients of the polynomial
#     totals = zeros(Float64,1,n)
#     for k in 1:n
#         cons = 1
#         x1 = [1, 0]
#         for j in 1:n
#             if j!=k
#                 x1 = foil(x1, [-xi[j], 1])
#                 cons = cons*(xi[k]-xi[j])
#             end
#         end
#         x1 = (1/cons)*x1
#         totals = totals + yi[k]*x1
#     end
#     return totals
# end
# lag = generate_lagrange_coefs(xi, yi, 9)
# function gen_func_lag(x)
#     y = 0
#     for i in 0:8
#         y = y + lag[i+1]*(x^i) 
#     end
#     return y
# end

# println("Lagrange Coefficients (x^0 to x^10): ", generate_lagrange_coefs(xi, yi, 11))
# yvalsa = gen_func_lag.(xtest)
# plt5 = plot(xtest, yvals0, label = "og", fmt =:pdf)
# plt5 = plot!(xtest, yvalsa, label = "lagrange", fmt=:pdf)

# plt6 = plot(xtest, log10.(abs.(yvalsa-yvals01)), label = "lagrange error", fmt =:pdf)

# for p in [plt5, plt6]
#     savefig(p, "temp.pdf")
#     append_pdf!("testing2.pdf", "temp.pdf", cleanup=true)
# end

