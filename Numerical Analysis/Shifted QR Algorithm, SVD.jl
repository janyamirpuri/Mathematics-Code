using LinearAlgebra: norm
using LinearAlgebra: qr
using LinearAlgebra: diag
using LinearAlgebra: eigen
using LinearAlgebra: LowerTriangular
using LinearAlgebra: Diagonal
using LinearAlgebra: I
using Plots: plot
using Plots: display
using Plots: savefig
using Plots


# Finding Singular Values of Matrix with Shifted QR Algorithm

function bidiagmatrix(n)
    A = zeros(Float64, n, n)
    for i in range(1, n-1)
        A[i,i] = sqrt(i)
        A[i,i+1] = (1+i)^(0.03)
    end
    A[n,n] = sqrt(n)
    return A
end

function permutation_mat(n)
    P = zeros(Float64, 2*n, 2*n)
    for i in 1:n
        P[2*i-1, i] = 1
        P[2*i, n+i] = 1
    end
    return P
end

function shifted_QR_alg(A,n)
    A_prev = copy(A)
    A_update = A_prev
    eig = zeros(n)
    if n == 1
        eig[1] = A_prev[1,1]
    else
        not_done = true
        while not_done
            A_prev = A_update
            mu = A_prev[n,n]
            (Q, R) = qr(A_prev-mu*I)
            A_update = R*Q + mu*I
            if norm(A_prev[n,1:n-1]) <= 1e-10 && abs(A_prev[n,n] - A_update[n,n]) <= 1e-6
                not_done = false
            end
        end
        eig = [A_prev[n,n]; shifted_QR_alg(A_prev[1:n-1, 1:n-1], n-1)]
    end
    return eig
end

function svd_eval(A, n) #assuming A is already bidiagonal
    println("Computing Singular Values for n = ", n ," Bidiagonal Matrix...")
    z = zeros(Float64, n, n)
    H = [z A'; A z]
    P = permutation_mat(n)
    T = P*H*P'
    Eig = shifted_QR_alg(T,2*n)
    z = sort(abs.(Eig), rev=false)
    eigenvals = [z[2*i] for i in 1:n]
    println(eigenvals)
    println()
    return eigenvals
end



A = bidiagmatrix(25)
si = svd_eval(A, 25)
println()
tests = [10,20,40,80,160,320,640]
results = [@elapsed svd_eval(bidiagmatrix(i), i) for i in tests]
println(results)


plt = plot(tests, results, label = "Run Time", fmt=:pdf)
plt = Plots.xaxis!("n")
plt = Plots.yaxis!("Run Time (seconds)")
plt = Plots.title!("SVD of nxn Bidiagonal Matrix Run Time")
savefig(plt, "bidiag.pdf")




