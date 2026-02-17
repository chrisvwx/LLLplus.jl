"""
    B,Z = bkz(H,K)

Block Korkine-Zolotarev (BKZ) reduction of the lattice (matrix) H with
blocksize K.

BKZ was originaly described in "A hierarchy of polynomial time basis
reduction algorithms" by Claus Peter Schnorr, Theoretical Computer Science,
1987. This function is not well tested and could be wrong or otherwise
different from Schnorr's BKZ.

Unfortunately I have forgotten from what paper or source I go the BKZ
algorithm implemented below.

This function appears to give reduced bases than are significantly worse
than those of fplll's BKZ.
# Examples
```
julia> H = latticegen_qary_b(24,12,10,1);

julia> B=LLLplus.bkz(H,2);

```
"""
function bkz(H::AbstractArray{Td,2},K) where {Td<:Number}
    m,n = size(H)
    Ti= LLLplus.getIntType(Td)

    B,Z = lll(H)
    Ti = eltype(Z)
    k = 1
    clean = zeros(n)
    norms = zeros(n) .+ Inf
    prevclean = false
    while !all(clean .==2)
        kx = k:min(k+K-1,n)
        Bk = view(B,:,kx)
        Zk = view(Z,:,kx)
        if kx[end]==n
            Bh,Zh = hkz(Bk)
            bhnorms = diag(Bh'*Bh)
            ih = sortperm(diag(Bh'*Bh))
            bhnorms = bhnorms[ih]
            Bh = Bh[:,ih]
            Zh = Zh[:,ih]
            if all(bhnorms .>= norms[kx])
                clean[kx] .= prevclean ? 2 : 1
            else
                clean[kx] .= 0
                Bk[:,:] = Bh
                Zh[:,:] = Zh
                norms[kx] .=bhnorms
            end
            k=1
            prevclean = all(clean .>=0)
        else
            u = LLLplus.svpu(Bk)
            b = Bk*u
            bnrm = b'*b
            z = Zk*u
            kk =findfirst(u .!=zero(Ti))
            if kk==1
                if (abs(u[1])==one(Ti) && all(u[2:end] .==zero(Ti))) || bnrm>=norms[k]
                    clean[k]= prevclean ? 2 : 1
                else
                    clean[k]=0
                    B[:,k] = b
                    Z[:,k] = z
                    norms[k] = bnrm
                end
            else
                if bnrm>=norms[k]
                    clean[k]= prevclean ? 2 : 1
                else
                    clean[k]=0
                    B[:,[k,k+kk-1]] = [b B[:,k]]
                    Z[:,[k,k+kk-1]] = [z Z[:,k]]
                    norms[[k,k+kk-1]] = [bnrm norms[k]]
                end
            end
            k=k+1
        end
    end
    return B,Z
end
