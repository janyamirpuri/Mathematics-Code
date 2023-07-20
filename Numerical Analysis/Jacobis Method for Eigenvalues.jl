#we are assuming that the input matrix A is symmetric

using LinearAlgebra

#function to find the location of the largest off diagonal value
function maxx(A, n)
    maxx = 0
    location = [0,0]
    z = diag(A)
    C = Diagonal(z)
    J = A-C
    for i in range(1, n)
        for j in range(1, n)
            if abs(J[i,j]) >= maxx
                maxx = abs(J[i,j])
                location = [i,j]
            end
        end
    end
    return location
end

#computing phi
function computephi(A, loc)
    apq = A[loc[1], loc[2]]
    app = A[loc[1], loc[1]]
    aqq = A[loc[2], loc[2]]
    phii = 0.5*atan(2apq,(aqq-app))
    return phii
end

#generating R_pq

function genR_pq(phi, loc, n)
    R_pq = Matrix{Float64}(LinearAlgebra.I, n, n)
    R_pq[loc[1], loc[1]] = cos(phi)
    R_pq[loc[1], loc[2]] = sin(phi)
    R_pq[loc[2], loc[1]] = -sin(phi)
    R_pq[loc[2], loc[2]] = cos(phi)
    return R_pq
end
#function to return the sqrt of sum of squares of off_diagonals
function zeroerror(A, n)
    z = diag(A)
    C = Diagonal(z)
    return norm(A - C)
end

function jacobi_eigen(A, n, tolerance)
    A_update = copy(A)
    not_done = true

    while not_done
        location = maxx(A_update, n)
        phii = computephi(A_update, location)
        R_pq = genR_pq(phii, location, n)
        A_update = R_pq'*A_update*R_pq
        if (zeroerror(A_update, n) <= tolerance)
            not_done = false
            break    
        end

    end
    return A_update
end

function sqrtmatrix(n)
    A = zeros(Float64, n, n)
    for i in range(1, n)
        for j in range(1, n)
            A[i,j] = sqrt((i^2)+(j^2))
        end
    end
    return A
end


for n in [10,20,40]
    A = sqrtmatrix(n)
    mat = jacobi_eigen(A, n, 1e-8)
    println("Top Eigenvalues with Largest Magnitude for ", n, " size matrix ", sort(diag(mat), by = abs, rev = true)[1:5])
    println()
    println("Top Most Positive Eigenvalues for ", n, " size matrix ", sort(diag(mat), rev = true)[1:5])
    println()

end

