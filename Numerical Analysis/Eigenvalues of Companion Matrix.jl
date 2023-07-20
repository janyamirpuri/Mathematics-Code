using Printf
using LinearAlgebra

function foil(x1, y1) #when multiplying any polynomial by (a + bx) ==> [a , b]
    x2 = zeros(Float64, 1, length(x1)+1)
    for i in range(1, length(x1))
        x2[i] = y1[1]*x1[i]
    end
    for i in range(1, length(x1))
        x2[i+1] = x2[i+1] + y1[2]*x1[i]
    end
    return x2
end

function pwcoefs(n)
    x1 = [-1, 1]
    for i in range(2,n)
        x1 = foil(x1, [-i, 1])
    end
    return x1
end

coefs15 = pwcoefs(15)

for i in range(1, length(coefs15))
    print("Coefficient of the x^", i-1, " term: ", coefs15[i], "\n")
end
print("\n")

A = zeros(Float64, 15, 15)

A[:, 15] = -1 * coefs15[1:15]

for i in range(2, 15)
    A[i,i-1] = 1
end



eigendata = LinearAlgebra.eigen(A)

for eig in eigendata.values
    println("Relative Error of Root ", eig, " is ", abs(eig - round(eig))/round(eig))
end



