using LinearAlgebra: opnorm
using LinearAlgebra: inv
using LinearAlgebra: cond

# Condition Number of Vandermonde Matrix
function generate_nodes(n)
    nodes = zeros(Float64, 1, n+1)
    for i in 0:n
        x_i = -1 + (2*i/n)
        nodes[1, i+1] = x_i
    end
    return nodes 
end

function generate_V(nodes, n)
    V = zeros(Float64, n+1, n+1)
    for i in 1:n+1
        node = nodes[i]
        for j in 1:n+1
            V[i,j] = node^(j-1)
        end
    end
    return V
end

function condition_num(A)
    magA = abs(opnorm(A, 2))
    magAIn = abs(opnorm(inv(A), 2))
    K = magA * magAIn
    return K
end




for i in [5, 15, 25]
    nodes = generate_nodes(i)

    x = generate_V(nodes, i)

    println("Condition Number for Transformation matrix, n = ", i, ": ", condition_num(x))

end

println()


