# comparisons between LLLplus and the LLL from Nemo

# Older versions of this file used Cxx.jl to run functions from fplll.
# Unfortunately Cxx.jl is no longer maintained and even when it was it was
# awkward to use.  The current version of the file focuses on Nemo and
# LLLplus.

# everything in this file is a bit fragile and not especially well explained

Pkg.add("Nemo")
Pkg.add("Plots")


########################################
# Int64 bases: plot of time vs dimenstion
########################################

import Nemo
using LLLplus
import Statistics
using Plots
gr();

macro mytwice(ex)
    quote
        local t0 = time_ns()
        local val0 = $(esc(ex))
        local t1 = time_ns()
        local val1 = $(esc(ex))
        local t2 = time_ns()
        min(t1-t0,t2-t1)/1e9
    end
end

# I typically copy-and-pasted the code from here to the plot cmd multiple
# times when changing things like the bit depth or largest dimension.


δ=.69;
η=.51
bits = 30;
ctx = Nemo.lll_ctx(δ,η)
Nsamp = 5
Nvec = Int64.(round.(2 .^(2:7)))
Nnum = length(Nvec);
times = zeros(2,Nnum,Nsamp);
small = zeros(2,Nnum,Nsamp);
for nx=1:Nnum
    N = Nvec[nx]
    k = Int64(floor(N/2))

    #b = zeros(Int128,N,N);
    b = zeros(Int64,N,N);

    S = Nemo.matrix_space(Nemo.ZZ, N,N);
    # Gtype=dataTypeForGram(bits,N);
    # @show Gtype

    for sx = 1:Nsamp
        b = LLLplus.gen_qary_b(Int64,N,k,bits) # no good for bits>31
        bn = Nemo.ZZMatrix(collect(b'));

        times[1,nx,sx] = @mytwice B,_ = LLLplus.lll(b,δ)
        small[1,nx,sx] = B[:,1]'*B[:,1]
        times[2,nx,sx] = @mytwice Bn, Tn = Nemo.lll_with_transform(bn,ctx);
        small[2,nx,sx] = 0.0;
        for ix = 1:N
            small[2,nx,sx] += Float64(Bn[1,ix]*conj(Bn[1,ix]))
        end

    end
tave  = Statistics.mean(times,dims=3)[:,:,1]
display(tave)
end

tave  = Statistics.mean(times,dims=3)[:,:,1];
save  = Statistics.mean(small,dims=3)[:,:,1];
[Nvec'; tave]
[Nvec'; save]

plot(Nvec,tave',label=permutedims(["lll (Julia) Int64",
                       "FLINT (C) GMP"]),
                       legend=:topleft,
                       legendfontsize=14,
                       linewidth=5,size=(700,525))
xaxis!("Basis dimension", :log10,guidefontsize=14)
yaxis!("time (sec)", :log10)
plot!(xticks = (Nvec, [string(ix) for ix ∈ Nvec]),
      xtickfontsize=14,ytickfontsize=14)
#savefig("timeVdim_$(bits)bitsInt64.png")




########################################
# Int64 bases: plot of time vs length of smallest basis as delta is varied
########################################

import Nemo
using LLLplus
import Statistics
using Plots
gr();

macro mytwice(ex)
    quote
        local t0 = time_ns()
        local val0 = $(esc(ex))
        local t1 = time_ns()
        local val1 = $(esc(ex))
        local t2 = time_ns()
        min(t1-t0,t2-t1)/1e9
    end
end



# I typically copy-and-pasted the code from here to the plot cmd multiple
# times when changing things like the bit depth or largest dimension.

δv=[collect(.29:.20:.89); .999]
#δv=[collect(.49:.15:.94); .999]
η=.51
Nδ = length(δv);
bits = 25
Nsamp = 20
N = 64
k = Int64(floor(N/2))
times = zeros(2,Nδ,Nsamp);
small = zeros(2,Nδ,Nsamp);
b = zeros(Int64,N,N);
S = Nemo.matrix_space(Nemo.ZZ, N,N);
# Gtype=dataTypeForGram(bits,N);

for nx=1:Nδ
    δ = δv[nx]
    ctx = Nemo.lll_ctx(δ, η)

    for sx = 1:Nsamp
        b = LLLplus.gen_qary_b(Int64,N,k,bits) # no good for bits>31

        bn = Nemo.ZZMatrix(collect(b'));

        times[1,nx,sx] = @mytwice Bl,T = LLLplus.lll(b,δ)
        small[1,nx,sx] = Bl[:,1]'*Bl[:,1]
        times[2,nx,sx] = @mytwice Bn, Tn = Nemo.lll_with_transform(bn,ctx);
        small[2,nx,sx] = 0.0;
        for ix = 1:N
            small[2,nx,sx] += Float64(Bn[1,ix]*conj(Bn[1,ix]))
        end
    end

tave  = Statistics.mean(times,dims=3)[:,:,1];
display([δv'; tave])
end
tave  = Statistics.mean(times,dims=3)[:,:,1];
save  = Statistics.mean(small,dims=3)[:,:,1];

plot(save',tave',label=permutedims(["lll (Julia) Int64",
                                    "FLINT (C) GMP"]),
                 legend=:topright,
                 legendfontsize=14,
                 xtickfontsize=14,ytickfontsize=14,
                 linewidth=5,size=(700,525),titlefontsize=14)
if N==16
    xaxis!("sqared 2-norm of smallest vector",
         :log10,guidefontsize=14,xlims=(10^7.3,10^9.8))
elseif N==32
    xaxis!("sqared 2-norm of smallest vector",
         :log10,guidefontsize=14,xlims=(10^7.5,10^10.2))
elseif N==64
    xaxis!("sqared 2-norm of smallest vector",
           :log10,guidefontsize=14,xlims=(10^8.0,10^11.1))
elseif N==128
    xaxis!("sqared 2-norm of smallest vector",
           :log10,guidefontsize=14,xlims=(10^9.6,10^14.4))
end
yaxis!("time (sec)", :log10)
