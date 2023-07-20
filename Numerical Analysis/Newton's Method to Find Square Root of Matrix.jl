using LinearAlgebra: norm

function mat_sqrt(A, n, tol)
    X_n = zeros(Float64, n, n)
    for i in 1:n
        X_n[i, i] = 1
    end
    not_done = true
    while not_done
        X_n = 0.5 * (X_n + A*inv(X_n))
        if norm((X_n*X_n - A), 2) < tol
            not_done = false
        end
    end
    return X_n 
end

A = zeros(Float64, 6, 6)
for i in 1:6
    for j in 1:6
        if i == j
            A[i,j] = 5 + sqrt((i^2 + j^2))
        else
            A[i,j] = sqrt((i^2 + j^2))
        end
    end
end

sqrrt = mat_sqrt(A, 6, 1e-8)
for i in 1:6
    println(sqrrt[i,:])
end
println()
print("The error of this is: ", norm((sqrrt*sqrrt - A), 2))