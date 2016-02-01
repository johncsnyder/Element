

# Hamiltonian operator
h(φ::Float64, Δφ::Float64, v::Float64) = -Δφ/2 + v*φ

function h!(dst::AbstractArray, φ::AbstractArray, Δφ::AbstractArray, v::AbstractArray)
    @assert length(dst) == length(φ) == length(Δφ) == length(v)
    for i in 1:length(dst)
        dst[i] = h(φ[i],Δφ[i],v[i])
    end
end

function h!(dst::AbstractMatrix, φ::AbstractMatrix, Δφ::AbstractMatrix, v::AbstractArray)
    @assert size(dst) == size(φ) == size(Δφ)
    @assert size(dst,1) == length(v)
    for i in 1:size(dst,2)
        for j in 1:size(dst,1)
            dst[j,i] = h(φ[j,i],Δφ[j,i],v[j])
        end
    end
end

function h(φ::AbstractArray, Δφ::AbstractArray, v::AbstractArray)
    dst = similar(φ)
    h!(dst,φ,Δφ,v)
    dst
end

function h!(dst::SparseMatrixCSC, φ::SparseMatrixCSC, Δφ::SparseMatrixCSC, v::AbstractArray)
    @assert size(dst) == size(φ) == size(Δφ)
    @assert size(dst,1) == length(v)
    for i in 1:length(φ.nzval)
        j = φ.rowval[i]
        dst.nzval[i] = h(φ.nzval[i],Δφ.nzval[i],v[j])
    end
end

function h(φ::SparseMatrixCSC, Δφ::SparseMatrixCSC, v::AbstractArray)
    dst = SparseMatrixCSC(φ.m,φ.n,φ.colptr,φ.rowval,similar(φ.nzval))
    h!(dst,φ,Δφ,v)
    dst
end



