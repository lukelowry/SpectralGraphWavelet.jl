using SparseArrays
using LinearAlgebra
using SuiteSparse

"""
    VFKern

Vector Fitting Kernel structure containing poles (Q) and residues (R).
"""
struct VFKern
    Q::Vector{Float64}
    R::Vector{Float64}
end

function impulse(L::AbstractSparseMatrix, n::Int=1, ntime::Int=1)
    b = zeros(Float64, size(L, 1), ntime)
    b[n, :] .= 1.0
    return b
end

# Regarding the inplce solve (ldiv!)
#       Seems to significantly improve memory usage in some cases 
# NOTE Could improve further by optimizing workspace Y and E but
#       Julia hides those parameters it passes to CHOLMOD, 
#       so unless we want to rewrite the Julia port, we are out of of luck there.

"""
    DyConvolve

Struct to hold the dynamic convolution state, including the Laplacian, poles, and factorizations.
"""
mutable struct DyConvolve
    nBus::Int
    L::SparseMatrixCSC{Float64, Int}
    poles::Vector{Float64}
    factors::Vector{SuiteSparse.CHOLMOD.Factor{Float64}}

    function DyConvolve(L::SparseMatrixCSC{Float64, Int}, poles::Vector)
        n = size(L, 1)
        factors = SuiteSparse.CHOLMOD.Factor{Float64}[]
        
        for q in poles
            F = ldlt(Symmetric(L); shift=q)
            push!(factors, F)
        end
        
        new(n, L, poles, factors)
    end
end

DyConvolve(L::SparseMatrixCSC, K::VFKern) = DyConvolve(L, K.Q)

(conv::DyConvolve)(B, K::VFKern) = convolve(conv, B, K)

function convolve(conv::DyConvolve, B, K::VFKern)
    nDim = size(B, 2)
    rows = size(B, 1)
    W = zeros(Float64, rows, nDim)
    
    # Pre-allocate scratch buffer once
    X1 = similar(B)

    for (i, (q, r)) in enumerate(zip(K.Q, K.R))
        F = conv.factors[i]
        
        # In-place solve: X1 = F \ B
        ldiv!(X1, F, B)
        
        W .+= X1 .* r
    end

    return W
end

function lowpass(conv::DyConvolve, B)
    W = typeof(B)[]
    
    # Pre-allocate buffer
    X1 = similar(B)
    
    for (i, q) in enumerate(conv.poles)
        F = conv.factors[i]
        
        # In-place solve
        ldiv!(X1, F, B)
        
        # Copy to output and scale
        # We must copy() here because X1 is reused in the next iteration
        res = copy(X1)
        res .*= q 
        push!(W, res)
    end

    return W
end

function bandpass(conv::DyConvolve, B)
    W = typeof(B)[]
    
    # Pre-allocate buffers
    X1 = similar(B)
    X2 = similar(B)
    
    for (i, q) in enumerate(conv.poles)
        F = conv.factors[i]
        
        # Step 1: X2 = F \ B
        ldiv!(X2, F, B)
        
        # Step 2: X1 = F \ X2
        ldiv!(X1, F, X2)
        
        # Step 3: X2 = L * X1
        # Use in-place matrix multiplication
        mul!(X2, conv.L, X1)
        
        # Scale
        X2 .*= (4.0 * q)
        push!(W, copy(X2))
    end

    return W
end

function highpass(conv::DyConvolve, B)
    W = typeof(B)[]

    # Pre-allocate buffers
    X1 = similar(B)
    X2 = similar(B)

    for (i, q) in enumerate(conv.poles)
        F = conv.factors[i]
        
        # Step 1: X1 = F \ B
        ldiv!(X1, F, B)
        
        # Step 2: X2 = L * X1
        mul!(X2, conv.L, X1)
        
        push!(W, copy(X2))
    end

    return W
end

function addbranch!(conv::DyConvolve, i::Int, j::Int, w::Float64)
    conv.L[i, j] -= w
    conv.L[j, i] -= w
    conv.L[i, i] += w
    conv.L[j, j] += w

    ws = sqrt(w)
    C = sparsevec([i, j], [ws, -ws], conv.nBus)

    for F in conv.factors
        SuiteSparse.CHOLMOD.updown!(F, C, true)
    end

    return true
end