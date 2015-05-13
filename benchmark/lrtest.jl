# @doc doc"""function lrtest(Ns::Int,N::Array{Int,1},L::Array{Int,1},
#                dataType::Array{DataType,1},distType)
#   Monte-Carlo evaluation of lattice reduction techniques. NxN matrices of
#   dataType are decomposed with each LR technique, with execution time,
#   memory, orthogonality defect, and Hermite factor shown.  Accepted
#   datatypes include Integer types such as Int8; FloatingPoint types such
#   as Float32; and Complex versions of these types. Accepted distTypes are
#   "randn" (Gaussian distributed with variance L), "rand" (uniformly
#   distributed over -L:L), "colUnif" (first column is uniform over -L:L).
#
# Examples
#   # As function of Integer width
#   lrtest(50,2.^[4;],[100;],[Int32,Int64,Int128,Float64,BigInt,BigFloat],"rand")
#   # As function of N
#   lrtest(50,2.^[0:6;],[100;],[Float64],"rand")
# """
function lrtest(Ns::Int,N::Array{Int,1},L::Array{Int,1},
                dataType::Array{DataType,1},distType)
#                dataType,distType)

# Packages that need to be loaded for lrtest to work include PyPlot, and
# CPUTime.
    
#lrAlgs = [lll, lllrecursive,seysen]
lrAlgs = [lll, seysen]

@printf("      Ns      N      L   dataType")
for ax = 1:min(length(lrAlgs),6)
    @printf(" %13s",lrAlgs[ax])
end
@printf("\n")

xtickStrs = [];
out = cell(0)
if length(N)>1
    for s = 1:length(N)
        push!(out,lrsim(Ns,N[s],L[1],dataType[1], distType, lrAlgs));
    end
    plotfun = PyPlot.loglog;
    xval = N;
    xlab = "Matrix Size";
    tstr = @sprintf("Ns=%d,L=%d,Type=%s,dist=%s",
                    Ns,L[1],string(dataType[1]),distType);
elseif length(L)>1
    for s = 1:length(L)
        push!(out,lrsim(Ns,N[1],L[s],dataType[1], distType, lrAlgs));
    end
    plotfun = semilogy;
    xval = L;
    xlab = "L";
    tstr = @sprintf("Ns=%d,N=%d,Type=%s,dist=%s",
                    Ns,N[1],string(dataType[1]),distType);
elseif length(dataType)>1
    for s = 1:length(dataType)
        push!(out,lrsim(Ns,N[1],L[1],dataType[s], distType, lrAlgs));
    end
    plotfun = semilogy;
    xval = 1:length(dataType)
    xtickStrs = map(string,dataType)
    xlab = "dataType";
    tstr = @sprintf("Ns=%d,N=%d,L=%d,dist=%s",Ns,N[1],L[1],distType);
end

pColor = {"r-","b.-","k-","g-","c-","m-",
          "r--","b.-.","k--","g.-.","c--","m*-.",
          "rs--","gp-.","bv-","kh--","c+-.","m.-",};
pIdx   = 1;

Nout = size(out,1);

f=figure(num=1,figsize=(6.5,4.5))
clf()
plt.style[:use]("ggplot")

times = zeros(Nout);
for a=1:length(out[1][2])
    for k=1:Nout
        times[k] = out[k][1][a];
    end
    plotfun(xval,times,pColor[pIdx],label=out[1][2][a]);
    hold(true);
    pIdx = pIdx==length(pColor)? 1: pIdx + 1;
end
hold(false);

xlabel(xlab);
ylabel("execution time (sec)");
legend(loc=2);
#title(tstr)
if ~isempty(xtickStrs)
    xticks(xval,xtickStrs)
end
grid(figure=f,which="both",axis="y")
grid(figure=f,which="major",axis="x")
plt.tight_layout()

return
end

#######################################################################
function lrsim(Ns,N,L,dataType,distType,lrAlgs)
# sortDict = Dict("Insertion" => d->sort(d,alg=InsertionSort),

distDict =     {"rand"=> d -> rand!(-L:L,d),
                "randmax"=> d -> rand(dataType,N,N),
                "randn"=> d -> L*randn!(d::AbstractArray{dataType,2}),
                 };
# distDict =     {"rand"=> d -> rand!(-L:L,d),
#                 "randmax"=> d -> rand!(dataType,d),
#                 "randn"=> d -> L*randn!(d::AbstractArray{dataType,2}),
#                  };
randf! = distDict[distType];

Nalgs = length(lrAlgs);
times = ones(Nalgs,Ns)*Inf;
hermitef = ones(Nalgs,Ns)*Inf;
orthf = ones(Nalgs,Ns)*Inf;
algNames = cell(Nalgs)

data = Array(dataType,N,N);

for ix = 1:Ns
    cnum = Inf;
    while cnum>1e5  # With Ints need to check for singular matrices
        data = randf!(data);
        if !(dataType<:BigInt || dataType<:BigFloat)
            cnum = cond(data)
        else
            cnum = 0
        end
    end

    #
    # Lattice-reduction algorithms
    #
    for ax = 1:length(lrAlgs)
        CPUtic();
        (B,T) = lrAlgs[ax](data);
        times[ax,ix] = CPUtoq();
        # detB = abs(det(B))
        # # Hermite factor
        # hermitef[ax,ix] = norm(B[:,1])/detB^(1/N)
        # # Orthogonality defect
        # orthf[ax,ix] = prodBi(B,N)/detB

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
@printf("%8d %6d %6d %10s", Ns,  N,  L, string(dataType))

#outtimes = mean(times,2);
outtimes = zeros(Nalgs,1);
outHF = zeros(Nalgs,1);
outOF = zeros(Nalgs,1);
for ax = 1:Nalgs
#    stimes = sort(vec(times[ax,:]));
#    outtimes[ax] = mean(stimes[1:floor(Ns*2/3)]);
    outtimes[ax] = mean(times[ax,:]);
    outHF[ax] = mean(hermitef[ax,:]);
    outOF[ax] = mean(orthf[ax,:]);
end

# BOLD= "\x1B[1;30m"
# RESET="\x1B[0;0m"
# (mn,mx)=findmin(outtimes)
for ax = 1:min(Nalgs,7)
    # if ax==mx
    #     @printf("    %s%10.4f%s",BOLD,outtimes[ax]*1000,RESET)
    # else
        @printf("    %10.4f",outtimes[ax]*1000)
#        @printf("    %10.4f",outHF[ax])
    # end 
end
@printf("\n")

return outtimes, algNames
end

#######################################################################
function prodBi(B,N)
    prod = 1.0
    for n=1:N
        pn = 0.0
        for m=1:N
            pn+=conj(B[m,n])*B[m,n]
        end
        println("pn=$(pn), typeof(pn)=$(typeof(pn))")
        prod*=sqrt(pn)
    end
    return prod
end
