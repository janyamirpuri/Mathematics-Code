using LinearAlgebra: norm
using LinearAlgebra: dot
using LinearAlgebra: I
using LinearAlgebra: cond

function classical_gs(n, A)
    nums = 1:n
    Q = zeros(Float64, n, n)
    for j in nums
        v_j = A[:,j]
        nums2 = 1:j-1
        for i in nums2
            r = dot(Q[:,i], A[:,j])
            v_j = v_j - r*Q[:,i]
        end
        Q[:, j] = v_j/(norm(v_j))
    end
    return Q
end

println("Classical GS Test for Random Matrix, n = 10:")
n = 10
A = rand(n,n)
Q = classical_gs(n, A)
err = norm(transpose(Q) * Q - I)

println("Error: " , err)

function modified_gs(n, A)
    nums = 1:n
    V = zeros(Float64, n, n)
    Q = zeros(Float64, n, n)
    for i in nums
        V[:,i] = A[:,i]
    end
    for i in nums
        Q[:,i] = V[:,i]/(norm(V[:,i]))
        nums2 = i+1:n
        for j in nums2
            r = dot(Q[:,i], V[:,j])
            V[:,j] = V[:,j] - r*Q[:,i]
        end
    end
    return Q
end
println()
println("Modified GS Test for the same Random Matrix, n = 10:")

Q = modified_gs(n, A)
err = norm(transpose(Q) * Q - I)

println("Error: " , err)
println()

function hilbert(n)
    H = zeros(Float64, n, n)
    nums = 1:n
    for i in nums
        for j in nums
            H[i,j] =1/(i+j-1)
        end
    end
    return H
end



prac = [25,50,100,200,400]

for z in prac
    H = hilbert(z)
    Q1= classical_gs(z, H)
    Q2= modified_gs(z, H)
    err1 = norm(transpose(Q1) * Q1 - I)
    err2 = norm(transpose(Q2) * Q2 - I)

    println("n = " ,z,"\nCondition Number of Matrix (2-Norm): ", cond(H,2), "\nClassical GS err: " ,err1,"\nModified GS err: " ,err2, "\nDifference in errors: ", err1-err2)
    println()
end

