# Cxx.jl may require Julia v1.2.x; it hasn't worked on all versions of Julia
# Pkg.add("Cxx")
# Pkg.add("Nemo")
# # Pkg.add("LLLplus")
# Pkg.add("Plots")



########################################
# Int64 bases: plot of time vs dimenstion
########################################

using Cxx
using LoopVectorization
using Libdl
import Nemo
using LLLplus
import Statistics
using Plots
gr();
include("src/l2.jl")

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

# following assumes fplll is in a sibling dir. It could also be in a system directory.
const path_to_lib = "../lll/fplll/"
addHeaderDir(path_to_lib, kind=C_System)
addHeaderDir("/usr/local/include", kind=C_System)
Libdl.dlopen(path_to_lib * "/lib/libfplll.dylib", Libdl.RTLD_GLOBAL)

cxxinclude("fplll.h")
cxxinclude("fplll/main.h")
cxx"""#include <iostream>"""


cxx"""ZZ_mat<long> b;"""
cxx"""ZZ_mat<long> u;"""
cxx"""ZZ_mat<long> u_inv;"""
cxx"""Options o;"""
cxx"""int status, flags;"""

icxx"""flags = LLL_DEFAULT;""";
icxx"""o.action = ACTION_LLL;""";
icxx"""o.delta = LLL_DEF_DELTA;"""; # 0.99
icxx"""o.eta = LLL_DEF_ETA;""";     # 0.51; see defs.h
icxx"""o.method = LM_FAST;""";    # see defs.h
icxx"""o.int_type = ZT_LONG;""";      # use double  (Int64?). See defs.h line 197
icxx"""o.float_type = FT_DOUBLE;""";  # use double  (Float64?). See defs.h line 206
# ?does generating or passing out unimodular slow the lll down?
function fplll()
    icxx"""status=lll_reduction(b,u,u_inv,o.delta,o.eta,o.method,o.float_type,o.precision,flags);"""
end

# I typically copy-and-pasted the code from here to the plot multiple times
# when changing things like the bit depth or largest dimension.


δ=.99;
η=.51
bits = 25;
ctx = Nemo.lll_ctx(δ,η)
Nsamp = 5
Nvec = Int64.(round.(2 .^(2:7)))
Nnum = length(Nvec);
times = zeros(4,Nnum,Nsamp);
small = zeros(4,Nnum,Nsamp);
for nx=1:Nnum
    N = Nvec[nx]
    k = Int64(floor(N/2))

    icxx"""b.resize($N,$N);""";
    icxx"""u.resize($N,$N);""";

    #b = zeros(Int128,N,N);
    b = zeros(Int64,N,N);
    bfp = zeros(Int64,N,N);

    S = Nemo.MatrixSpace(Nemo.ZZ, N,N);
    Gtype=dataTypeForGram(bits,N);
    @show Gtype

    for sx = 1:Nsamp
        # b = LLLplus.gen_qary_b(Int64,N,k,bits) # no good for bits>31
        # for ix = 1:N, jx = 1:N
        #      icxx"""b[$(ix-1)][$(jx-1)]= $(b[jx,ix]);"""
        # end
        icxx"""b.gen_qary($k,$bits);""";
        for ix = 1:N, jx = 1:N
            b[jx,ix] = icxx"""b[$(ix-1)][$(jx-1)].get_si();"""
        end
        bn = S(Nemo.fmpz.(b')); # fmpz is GMP

        times[1,nx,sx] = @mytwice B,T = LLLplus.lll(b,δ)
        small[1,nx,sx] = B[:,1]'*B[:,1]
#        times[2,nx,sx] = @mytwice B = LLLplus.l2(b,Gtype,δ,η)
        times[2,nx,sx] = @mytwice B = l2avx(b,Gtype,δ,η)
        small[2,nx,sx] = B[:,1]'*B[:,1]
        times[3,nx,sx] = @mytwice Bn, Tn = Nemo.lll_with_transform(bn,ctx);
        small[3,nx,sx] = 0.0;
        for ix = 1:N
            small[3,nx,sx] += Float64(Bn[1,ix]*Bn[1,ix]')
        end
        times[4,nx,sx] = @elapsed fplll()
        for ix = 1:N, jx = 1:N
            bfp[jx,ix] = icxx"""b[$(ix-1)][$(jx-1)].get_si();"""
        end
        small[4,nx,sx] = bfp[:,1]'*bfp[:,1]

    end
tave  = Statistics.mean(times,dims=3)[:,:,1]
display(tave)
end

tave  = Statistics.mean(times,dims=3)[:,:,1];
save  = Statistics.mean(small,dims=3)[:,:,1];
[Nvec'; tave]
[Nvec'; save]

plot(Nvec,tave',label=permutedims(["lll (Julia) Int64",
                       "l2avx (Julia) Int64",
                       "FLINT (C) GMP",
                       "fplll (C++) Int64"]),
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

using Cxx
using LoopVectorization
using Libdl
import Nemo
using LLLplus
import Statistics
using Plots
gr();
include("src/l2.jl")

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



# following assumes fplll is in a sibling dir. It could also be in a system directory.
const path_to_lib = "../lll/fplll/"
addHeaderDir(path_to_lib, kind=C_System)
addHeaderDir("/usr/local/include", kind=C_System)
Libdl.dlopen(path_to_lib * "/lib/libfplll.dylib", Libdl.RTLD_GLOBAL)

cxxinclude("fplll.h")
cxxinclude("fplll/main.h")
cxx"""#include <iostream>"""


cxx"""ZZ_mat<long> b;"""
cxx"""ZZ_mat<long> u;"""
cxx"""ZZ_mat<long> u_inv;"""
cxx"""Options o;"""
cxx"""int status, flags;"""

icxx"""flags = LLL_DEFAULT;""";
icxx"""o.action = ACTION_LLL;""";
icxx"""o.eta = LLL_DEF_ETA;""";     # 0.51; see defs.h
icxx"""o.method = LM_FAST;""";    # see defs.h
icxx"""o.int_type = ZT_LONG;""";      # use double  (Int64?). See defs.h line 197
icxx"""o.float_type = FT_DOUBLE;""";  # use double  (Float64?). See defs.h line 206
# ?does generating or passing out unimodular slow the lll down?
function fplll()
    icxx"""status=lll_reduction(b,u,u_inv,o.delta,o.eta,o.method,o.float_type,o.precision,flags);"""
end

# I typically copy-and-pasted the code from here to the plot multiple times
# when changing things like the bit depth or largest dimension.

δv=[collect(.29:.20:.89); .999]
#δv=[collect(.49:.15:.94); .999]
η=.51
Nδ = length(δv);
bits = 25
Nsamp = 10
N = 32
k = Int64(floor(N/2))
times = zeros(4,Nδ,Nsamp);
small = zeros(4,Nδ,Nsamp);
icxx"""b.resize($N,$N);""";
icxx"""u.resize($N,$N);""";
b = zeros(Int64,N,N);
bfp = zeros(Int64,N,N);
S = Nemo.MatrixSpace(Nemo.ZZ, N,N);
Gtype=dataTypeForGram(bits,N);

for nx=1:Nδ
    δ = δv[nx]
    ctx = Nemo.lll_ctx(δ, η)

    icxx"""o.delta = $(δ);""";

    for sx = 1:Nsamp

        icxx"""b.gen_qary($k,$bits);""";
        for ix = 1:N, jx = 1:N
            b[jx,ix] = icxx"""b[$(ix-1)][$(jx-1)].get_si();"""
        end
        bn = S(Nemo.fmpz.(b')); # fmpz is GMP

        times[1,nx,sx] = @mytwice Bl,T = LLLplus.lll(b,δ)
        small[1,nx,sx] = Bl[:,1]'*Bl[:,1]
#        times[2,nx,sx] = @mytwice B2 = LLLplus.l2(b,Gtype,δ,η)
        times[2,nx,sx] = @mytwice B2 = l2avx(b,Gtype,δ,η)
        small[2,nx,sx] = B2[:,1]'*B2[:,1]
        times[3,nx,sx] = @mytwice Bn, Tn = Nemo.lll_with_transform(bn,ctx);
        small[3,nx,sx] = 0.0;
        for ix = 1:N
            small[3,nx,sx] += Float64(Bn[1,ix]*Bn[1,ix]')
        end
        times[4,nx,sx] = @elapsed fplll()
        for ix = 1:N, jx = 1:N
            bfp[jx,ix] = icxx"""b[$(ix-1)][$(jx-1)].get_si();"""
        end
        small[4,nx,sx] = bfp[:,1]'*bfp[:,1]
    end

tave  = Statistics.mean(times,dims=3)[:,:,1];
display([δv'; tave])
end
tave  = Statistics.mean(times,dims=3)[:,:,1];
save  = Statistics.mean(small,dims=3)[:,:,1];

plot(save',tave',label=permutedims(["lll (Julia) Int64",
                                    "l2avx (Julia) Int64",
                                    "FLINT (C) GMP",
                                    "fplll (C++) Int64"]),
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
# plot(save',tave',label=permutedims(["LLLplus.lll Int64",
#                                     "LLLplus.l2+avx Int64",
#                                     "FLINT (C++) GMP",
#                                     "fplll (C++) Int64"]),
#                           linewidth=3,legend=:topright,
#                           legendfontsize=13,
#                           xtickfontsize=13,ytickfontsize=13)
# xaxis!("sqared 2-norm of smallest vector", :log10,guidefontsize=13)
#yaxis!("time (sec)", :log10,ylims=(10^-2.7,10^-1.3))
yaxis!("time (sec)", :log10)
#title!("fplll.gen_qary($k,$bits) used for random basis")
#savefig("timeVsmallest_$(bits)bitsInt64.png")

