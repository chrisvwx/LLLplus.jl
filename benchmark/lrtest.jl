"""
    lrtest(Ns::Int, N::Array{Int,1}, L::Array{Int,1}, dataType::Array{DataType,1}, distType)

Monte-Carlo evaluation of lattice reduction techniques. NxN matrices of
dataType are decomposed with each LR technique, with execution time,
memory, orthogonality defect, and Hermite factor shown.  Accepted
datatypes include Integer types such as Int8; FloatingPoint types such
as Float32; and Complex versions of these types. Accepted distTypes are
"randn" (Gaussian distributed with variance L), "rand" (uniformly
distributed over -L:L), "colUnif" (first column is uniform over -L:L).

# Examples
  # As function of Integer width
  lrtest(50,2.^[4;],[100;],[Int32,Int64,Int128,Float64,BigInt,BigFloat],"rand")
  # As function of N
  lrtest(50,2.^[0:6;],[100;],[Float64],"rand")
"""
function lrtest(Ns::Int,N::Array{Int,1},L::Array{Int,1},
                dataType::Array{DataType,1},distType)

#lrAlgs = [lll, lllrecursive,seysen]
lrAlgs = [lll, seysen]

@printf("      Ns      N      L   dataType")
for ax = 1:min(length(lrAlgs),6)
    @printf(" %13s",lrAlgs[ax])
end
@printf("\n")

xtickStrs = [];
out = []
if length(N)>1
    for s = 1:length(N)
        push!(out,lrsim(Ns,N[s],L[1],dataType[1], distType, lrAlgs));
    end
    xscale=:log10;
    yscale=:log10;
    xval = N;
    xlab = "Matrix Size";
    tstr = @sprintf("Ns=%d,Type=%s,dist=%s",
                    Ns,string(dataType[1]),distType);
elseif length(L)>1
    for s = 1:length(L)
        push!(out,lrsim(Ns,N[1],L[s],dataType[1], distType, lrAlgs));
    end
    xscale = :identity;
    yscale = :log10;
    xval = L;
    xlab = "L";
    tstr = @sprintf("Ns=%d,N=%d,Type=%s,dist=%s",
                    Ns,N[1],string(dataType[1]),distType);
elseif length(dataType)>1
    for s = 1:length(dataType)
        push!(out,lrsim(Ns,N[1],L[1],dataType[s], distType, lrAlgs));
    end
    xscale = :identity;
    yscale = :log10;
    xval = 1:length(dataType)
    xtickStrs = map(string,dataType)
    xlab = "dataType";
    tstr = @sprintf("Ns=%d,N=%d,L=%d,dist=%s",Ns,N[1],L[1],distType);
    for n=xval
        if dataType[n]== DoubleFloat{Float64}
            xtickStrs[n] = "Double64"
        end
    end
end

pColor = ["r-","b.-","k-","g-","c-","m-",
          "r--","b.-.","k--","g.-.","c--","m*-.",
          "rs--","gp-.","bv-","kh--","c+-.","m.-",];
pIdx   = 1;

Nout = size(out,1);


times = zeros(Nout);
orthf = zeros(Nout);
pltT = plot(legend=:topleft);
#pltO = plot(legend=:topleft);

for a=1:length(out[1][2])
    for k=1:Nout
        times[k] = out[k][2][a];
        orthf[k] = out[k][3][a];
    end
    plot!(pltT,xval,times,xscale=xscale,yscale=yscale,label=out[1][1][a],linewidth=3);
#    plot!(pltO,xval,orthf,xscale=xscale,yscale=yscale,label=out[1][1][a],linewidth=3);
    pIdx = pIdx==length(pColor) ? 1 : pIdx + 1;
end

plot!(pltT,xlabel=xlab,ylabel="execution time (sec)")
#plot!(pltO,xlabel=xlab,ylabel="orthogonalization factor")

if ~isempty(xtickStrs)
    plot!(pltT,xticks=(xval,xtickStrs))
#    plot!(pltO,xticks=(xval,xtickStrs))
end

display(pltT)
#display(pltO)

end


function lrsim(Ns,N,L,dataType,distType,lrAlgs)
# sortDict = Dict("Insertion" => d->sort(d,alg=InsertionSort),

distDict = Dict("rand"=> d -> rand!(d,-L:L),
                "randmax"=> d -> rand(dataType,N,N),
                "randn"=> d -> L*randn!(d::AbstractArray{dataType,2}),
                 );
randf! = distDict[distType];

Nalgs = length(lrAlgs);
times = ones(Nalgs,Ns)*Inf;
hermitef = ones(Nalgs,Ns)*Inf;
orthf = ones(Nalgs,Ns)*Inf;
#algNames = cell(Nalgs)
algNames = Array{Any,1}(undef,Nalgs)
data = Array{dataType,2}(undef,N,N);

for ix = 1:Ns
    cnum = Inf;
    while cnum > 1e5  # With Ints need to check for singular matrices
        data = randf!(data);
        if !(dataType<:BigInt || dataType<:BigFloat)
            cnum = cond(data)
        else
            cnum = 0.0
        end
    end

    #
    # Lattice-reduction algorithms
    #
    for ax = 1:length(lrAlgs)
        # CPUtic();
        # (B,T) = lrAlgs[ax](data);
        # times[ax,ix] = CPUtoq();
        times[ax,ix] = @elapsed (B,T) = lrAlgs[ax](data) 
        detB = abs(det(B))
        # Hermite factor
        hermitef[ax,ix] = norm(B[:,1])/detB^(1/N)
        # Orthogonality defect
        orthf[ax,ix] = prodBi(B,N)/detB

        # if abs(abs(det(T))-1.0)>1e-6
        #     println("For $(lrAlgs[ax]), det(T) is $(abs(det(T)))")
        #     times[ax,ix] = Inf
        # elseif abs(sum(B-data*T))>1e-6
        #     println("For $(lrAlgs[ax]), sum(B-H*T) "*
        #             "is $(abs(sum(B-data*T))), ")
        #     times[ax,ix] = Inf
        # end
    end
end

#algNames[Nalgs] = "randf";
for ax = 1:length(lrAlgs)
    algNames[ax] = string(lrAlgs[ax]);
end
@printf("%8d %6d %6d %10s", Ns,  N,  L, string(dataType)[1:min(end,10)])

#outtimes = mean(times,2);
outtimes = zeros(Nalgs,1);
outHF = zeros(Nalgs,1);
outOF = zeros(Nalgs,1);
mean(x) = sum(x)/length(x)
for ax = 1:Nalgs
#    stimes = sort(vec(times[ax,:]));
#    outtimes[ax] = mean(stimes[1:floor(Ns*2/3)]);
#    outtimes[ax] = minimum(times[ax,:]);
    outtimes[ax] = mean(times[ax,:]);
    outHF[ax] = mean(hermitef[ax,:]);
    outOF[ax] = mean(orthf[ax,:]);
end

for ax = 1:min(Nalgs,7)
    @printf("    %10.4f",outtimes[ax]*1000)
end
@printf("\n")

return algNames, outtimes, outOF, outHF
end

#######################################################################
function prodBi(B,N)
    prod = 1.0
    for n=1:N
        pn = 0.0
        for m=1:N
            pn+=conj(B[m,n])*B[m,n]
        end
        # println("pn=$(pn), typeof(pn)=$(typeof(pn))")
        prod*=sqrt(pn)
    end
    return prod
end
