#Calculating Machine Precisions for Different bit floating point numbers

#Machine Precision for 16 bit floating point numbers (half precision)

macheps = Float16(1.0)

for i = 1:100
    global macheps = Float16(2.0)^(-i)
    if (Float16(1.0) == Float16(1.0)+macheps)
        macheps = Float16(2.0)^(-i+1)
        println("Machine Precision for 16 bit floating point numbers = ", macheps)
        break
    end
end

#Machine Precision for 32 bit floating point numbers (half precision)

macheps = Float32(1.0)

for i = 1:100
    global macheps = Float32(2.0)^(-i)
    if (Float32(1.0) == Float32(1.0)+macheps)
        macheps = Float32(2.0)^(-i+1)
        println("Machine Precision for 32 bit floating point numbers = ", macheps)
        break
    end
end

#Machine Precision for 64 bit floating point numbers (half precision)

macheps = Float64(1.0)

for i = 1:100
    global macheps = Float64(2.0)^(-i)
    if (Float64(1.0) == Float64(1.0)+macheps)
        macheps = Float64(2.0)^(-i+1)
        println("Machine Precision for 64 bit floating point numbers = ", macheps)
        break
    end
end

