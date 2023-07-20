using LinearAlgebra


function pow_method(n, A)
    w = rand(n,1)
    w = w/norm(w,2)
    x = A * w
    not_done = true

    lam = 0
    lam_prev = lam
    lam_prev_prev = 0
    while not_done
        lam_prev_prev = lam_prev
        lam_prev = lam
        w = x/norm(x,2)
        lam = dot(x, w)
        x = A * w
        if abs(lam - lam_prev) <= 1e-12
            not_done = false
        end
    end
    rate = abs(lam -lam_prev)/abs(lam_prev-lam_prev_prev)
    return lam, w, rate
end

function inverse_with_shift_pow(n, A, s)
    #returns the largest eigenvalue and corresponding eigenvector of B = (A-sI)^-1
    B = inv((A - s*I))
    eg_valB, eg_vec, rate = pow_method(n, B)
    eg_valA = (1/eg_valB)+s
    return eg_valA, eg_vec, rate
end

epsilon = 0.25
A = zeros(Float64, 3,3)
A[1,1] = 1
A[2,2] = 3
A[3,3] = 6
A[1,3] = 2*epsilon
A[3,1] = 2*epsilon
A[2,3] = epsilon
A[3,2] = epsilon


#println(pow_method(3, A))
println()
egval, egvec, rate = inverse_with_shift_pow(3, A, 3)

println("The estimated eigenvalue is ", egval)
println("The estimated eigenvector is ", egvec)
println("The estimated rate of convergence is ", rate)

