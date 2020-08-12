var documenterSearchIndex = {"docs":
[{"location":"functions/#Functions","page":"Functions","title":"Functions","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"DocTestSetup = quote\n    using LLLplus\n    using LinearAlgebra\nend","category":"page"},{"location":"functions/","page":"Functions","title":"Functions","text":"    lll\n    cvp\n    svp\n    brun\n    gauss\n    seysen\n    vblast\n    subsetsum\n    lagariasodlyzko\n    mdsubsetsum\n    integerfeasibility\n    rationalapprox\n    spigotBBP\n    hard_sphere\n    issizereduced\n    islllreduced\n    orthogonalitydefect\n    hermitefactor\n    seysencond\n    gen_qary_b","category":"page"},{"location":"functions/#LLLplus.lll","page":"Functions","title":"LLLplus.lll","text":"B,T,Q,R = lll(H,δ=3/4)\n\nDo Lenstra–Lenstra–Lovász lattice reduction of matrix H using optional parameter δ.  The output is B, an LLL-reduced basis; T, a unimodular (meaning det(T)=+/-1) transformation matrix such that B= H*T; and finally Q and R which are a QR decomposition of B.  So H = B*inv(T) = Q*R*inv(T).\n\nFollows D. Wuebben, et al, \"Lattice Reduction - A Survey with Applications in Wireless Communications\", IEEE Signal Processing Magazine, 2011. When comparing with academic lattice reduction papers, it is likely closest to the floating-point algorithm of C. P. Schnorr. \"A more efficient algorithm for lattice basis reduction\". Journal of Algorithms, Vol 9, 1988.\n\nExamples\n\njulia> H= [1 2; 3 4];B,_ = lll(H); B\n2×2 Array{Int64,2}:\n 1  -1\n 1   1\n\njulia> H= BigFloat.([1.5 2; 3 4]) .+ 2im; B,_= lll(H); B\n2×2 Array{Complex{BigFloat},2}:\n 0.50+0.0im  0.0+1.0im\n  1.0+0.0im  0.0+0.0im\n\njulia> N=500;H = randn(N,N); B,T = lll(H);\n\n\n\n\n\n\n","category":"function"},{"location":"functions/#LLLplus.cvp","page":"Functions","title":"LLLplus.cvp","text":"x=cvp(y,R,infinite=Val{true},Umax=1)\n\nSolve the problem argmin_x ||z-Hx|| for integer x using the technique from the paper below, where H=QR and y=Q'*z. The input vector y is of length n, with H of dimension n by n, and the returned vector x of length n. If infinite==Val{true} then we search the (infinite) lattice, otherwise we search integers in [-Umax,Umax].  Note that cvp does not handle complex numbers.\n\nUses alg from \"Faster Recursions in Sphere Decoding\" Arash Ghasemmehdi, Erik Agrell, IEEE Transactions on Information Theory, vol 57, issue 6 , June 2011.\n\nExamples\n\njulia> H=[1 2; 3 4]; Q,R=qr(H); uhat = cvp(Q'*[0,2],R)\n2-element Array{Float64,1}:\n  2.0\n -1.0\n\njulia> n=100;H=randn(n,n);Q,R=qr(H);\n\njulia> u=Int.(rand(0:1e10,n));y=H*u+rand(n)/100;\n\njulia> uhat=cvp(Q'*y,R); sum(abs.(u-uhat))\n0.0\n\njulia> n=500;H=randn(n,n);Q,R=qr(H);\n\njulia> u=Int.(rand(-1:1,n));y=H*u+rand(n)/10;\n\njulia> uhat=cvp(Q'*y,R,Val{false},1); sum(abs.(u-uhat))\n0.0\n\n\n\n\n\n","category":"function"},{"location":"functions/#LLLplus.svp","page":"Functions","title":"LLLplus.svp","text":"b = svp(B)\n\nFind the shortest basis vector b for the lattice formed by the matrix B. This solves the 'shortest vector problem' (SVP). \n\nWe call the cvp function in the library n times for an n- dimensional lattice, so this is definitely not the fastest SVP solver :-) Roughly follows the CVP-to-SVP reduction in http://web.eecs.umich.edu/~cpeikert/lic15/lec06.pdf\n\nAlso, this function is not always correct. For example svp does not return a shortest vector for the following basis:\n\nB = [ 1    0    0    0\n      0    1    0    0\n    208  175  663    0\n    651  479    0  663];\n\nOne shortest vector for B is [-16; 19; -3; 11], the function currently returns [3; -11; 25; -1].\n\nExamples\n\njulia> H=[1 2; 3 4]; svp(H)\n2-element Array{Int64,1}:\n -1\n  1\n\njulia> H= BigFloat.([2.5 2; 3 4]); svp(H)\n2-element Array{BigFloat,1}:\n  0.50\n -1.0 \n\n\n\n\n\n\n","category":"function"},{"location":"functions/#LLLplus.brun","page":"Functions","title":"LLLplus.brun","text":"B, T = brun(H)\n\nBrun's integer-relations alg implemented as a matrix decomposition. Takes as input the matrix H and returns a reduced basis B and T, a unimodular transformation matrix such that B = H*T. Brun reduction is often done with pinv(H) as input to yield B = pinv(H)*T.\n\nSee V. Brun, \"En generalisation av kjedebrøken I,\" Skr. Vid ensk. Selsk. Kristiana, Mat. Nat. Klasse, 1919.  See https://archive.org/stream/skrifterutgitavv201chri#page/300/mode/2up\n\nFollows code from D. Wuebben, D. Seethaler, J. Jalden, and G. Matz, \"Lattice Reduction - A Survey with Applications in Wireless Communications\" IEEE Signal Processing Magazine, March 2011\n\nExamples\n\njulia> H=[1 2; 3 4]; B,T=brun(H); T\n2×2 Array{Int64,2}:\n  3  -1\n -2   1\n\n\n\n\n\n\n","category":"function"},{"location":"functions/#LLLplus.gauss","page":"Functions","title":"LLLplus.gauss","text":"B = gauss(H)\n\nDo Gauss/Lagrange reduction on the lattice defined by the two columns of H.\n\nFollows Fig 2.3 of \"Lattice Basis Reduction: An Introduction to the LLL Algorithm and Its Applications\" by Murray R. Bremner, CRC Press, 2012.\n\nExamples\n\njulia> H = [1 2; 3 3]; B = gauss(H)\n2×2 Array{Float64,2}:\n 1.0  0.0\n 0.0  3.0\n\n\n\n\n\n\n","category":"function"},{"location":"functions/#LLLplus.seysen","page":"Functions","title":"LLLplus.seysen","text":"B,T,B_dual,num_it = seysen(H::Array{Td,2}) where Td\n\nDo greedy  Seysen lattice reduction  on the  matrix H, returning  B, the reduced lattice basis;  T a unimodular matrix that reduces  H (i.e. B = H*T); B_dual, dual lattice basis (i.e., B_dual = pinv(B)); and num_it the number of iterations (basis updates). See also lll.\n\nFollows Seysen algorithm in \"Lattice Reduction - A Survey with Applications in Wireless Communications\" by D. Wuebben, et al, IEEE Signal Processing Magazine, 2011.\n\nExamples\n\njulia> H= [1 2; 3 4];B,T = seysen(H); B\n2×2 Array{Int64,2}:\n -1  1\n  1  1\n\njulia> H= BigFloat.([1.5 2; 3 4]) .+ 2im; B,_= seysen(H); B\n2×2 Array{Complex{BigFloat},2}:\n 0.0+1.0im  0.50+0.0im\n 0.0+0.0im   1.0+0.0im\n\n\n\n\n\n\n","category":"function"},{"location":"functions/#LLLplus.vblast","page":"Functions","title":"LLLplus.vblast","text":"W,P,B = vblast(H)\n\nFind a VBLAST decomposition of H such that H = pinv(W)*B*P' or B = W*H*P.  Here P is a permutation matrix, B is lower triangular with ones on the diagonal, and W has orthogonal rows.\n\nW,P,B = vblast(H,mu)\n\nIf an SNR argument mu is passed in, a regularized (\"MMSE\") decomposition is done, with the result that W will no longer have orthogonal rows and B is no longer lower triangular.\n\nExamples\n\njulia> H= [1. 2; 3 4];W,_ = vblast(H); W\n2×2 Array{Float64,2}:\n 1.5  -0.5\n 0.1   0.3\n\njulia> H= BigFloat.([1.5 2; 3 4]) .+ 2im; W,_= vblast(H); W\n2×2 Array{Complex{BigFloat},2}:\n      -2.0+3.0im            2.0-1.5im     \n 0.0779221-0.103896im  0.155844-0.103896im\n\n\n\n\n\n\n","category":"function"},{"location":"functions/#LLLplus.subsetsum","page":"Functions","title":"LLLplus.subsetsum","text":"x = subsetsum(a,s)\n\nFor a vector of integers a, and an integer s, try to find a binary vector x such that x'*a=s. We use the LLL algorithm to find the solution. This is not a robust tool, just a demo.\n\nThis function tries first the technique in the lagariasodlyzko function, and if it fails, a solution via mdsubsetsum is attempted.\n\nIt appears that this function can also solve some integer relations problems. See the first example.\n\nExamples\n\njulia> a=[1.5;.5;0;.1;.2]; s=2.2; x,_=subsetsum(a,s,true); s-x'*a\nA binary Lagarias-Odlyzko solution was found.\nA solution was found via lagariasodlyzko\n0.0\n\njulia> a=[32771,65543,131101,262187,524387,1048759, # from Bremner p 117\n          2097523,4195057,8390143,16780259,33560539,\n          67121039,134242091,268484171,536968403];\n\njulia> s=891221976; x,_=subsetsum(a,s,false); s-x'*a\n0.0\n\njulia> N=40;a=rand(1:2^BigInt(256),N);xtrue=rand(Bool,N); s=a'*xtrue; \n\njulia> setprecision(BigFloat,300); x,_=subsetsum(a,s,false); s-x'*a\n0.0\n\n\n\n\n\n\n","category":"function"},{"location":"functions/#LLLplus.lagariasodlyzko","page":"Functions","title":"LLLplus.lagariasodlyzko","text":"x = lagariasodlyzko(a,s)\n\nFor a vector of integers a, and an integer s, try to find a binary vector x such that x'*a=s. We use the LLL algorithm to find the solution. This is not a robust tool, just a demo.\n\nThis follows the technique described by Lagarias and Odlyzko  in  \"Solving Low-Density Subset Sum Problems\"  in Journal of ACM, Jan 1985. Code based on http://web.eecs.umich.edu/~cpeikert/lic15/lec05.pdf We can likely get better results using techniques described and referenced in https://www-almasty.lip6.fr/~joux/pages/papers/ToolBox.pdf\n\nIt's odd that permuting the a vector in the second example given below causes the alg to often not find a binary solution. Apparently this is a common oddity with lattice solvers.\n\nExamples\n\njulia> a=[1.5;.5;0;.1;.2]; s=2.2; x,_=lagariasodlyzko(a,s); s-x'*a\n0.0\n\njulia> a=[32771,65543,131101,262187,524387,1048759, # from Bremner p 117\n          2097523,4195057,8390143,16780259,33560539,\n          67121039,134242091,268484171,536968403];\n\njulia> s=891221976; x,_=lagariasodlyzko(a,s); s-x'*a\n0.0\n\njulia> N=40;a=rand(1:2^BigInt(256),N);xtrue=rand(Bool,N); s=a'*xtrue; \n\njulia> setprecision(BigFloat,300); x,_=lagariasodlyzko(a,s); s-x'*a\n0.0\n\n\n\n\n\n\n","category":"function"},{"location":"functions/#LLLplus.mdsubsetsum","page":"Functions","title":"LLLplus.mdsubsetsum","text":"x = mdsubsetsum(a,sM,ratio=.5,Kpm=3)\n\nFor a vector of integers a, and an integer sM, try to find a binary vector x such that x'*a=s using the technique from \"Multidimensional subset sum problem\" [1][2]. A major goal of the technique is to solve problems in which there are about 50% ones in x; other ratios of ones to zeros can be specified in ratio.  The thesis also suggests searching Kpm=3 values around the nominal k. This technique is related to that in subsetsum in that both use the LLL algorithm.  This is not a robust tool, just a demo.\n\n[1] https://scholarworks.rit.edu/theses/64/ [2] https://pdfs.semanticscholar.org/21a7/c2f9ff29507f1153aefcca04d1cd308e45c0.pdf\n\nExamples\n\njulia> a=[1.5;.5;0;.1;.2]; s=2.2; x=mdsubsetsum(a,s); s-x'*a\n0.0\n\njulia> a=[32771,65543,131101,262187,524387,1048759, # from Bremner p 117\n          2097523,4195057,8390143,16780259,33560539,\n          67121039,134242091,268484171,536968403];\n\njulia> sM=891221976; x=mdsubsetsum(a,sM); sM-x'*a\n0\n\njulia> setprecision(BigFloat,300); Random.seed!(0);\n\njulia> N=40;a=rand(1:2^BigInt(256),N);xtrue=rand(Bool,N); s=a'*xtrue;\n\njulia> x=mdsubsetsum(a,s); s-x'*a\n0\n\n\n\n\n\n","category":"function"},{"location":"functions/#LLLplus.integerfeasibility","page":"Functions","title":"LLLplus.integerfeasibility","text":"integerfeasibility(A,d,nullVecs=false)\n\nGiven a linear system Ax=d, return an integer vector x which satisfies the system.  This is the integer programming feasibility problem.\n\nIf nullVecs==true, then as well as returning a solution x, also return a matrix xNull of vectors in the null space of A which could be added to the x vector to find a solution which satisfies a constraint such as 0 .≤ x .≤ u; see the paper below.  \n\nThis is not a robust tool, just a demo.\n\n\"Solving A System Of Diophantine Equations With Bounds On The Variables\" by Karen Aardal, Cor Hurkens, and Arjen Lenstra in Integer Programming and Combinatorial Optimization, 6th International IPCO Conference, vol 1412, pp 229-242, 1998. See http://softlib.rice.edu/pub/CRPC-TRs/reports/CRPC-TR98782.pdf\n\nExamples\n\njulia> A=[10 1 -9; 1 8 8]; xtrue=[0; 2; 9]; d=A*xtrue;\n\njulia> integerfeasibility(A,d)\n3-element Array{Int64,1}:\n 0\n 2\n 9\n\njulia> A=[10 1.1 -9.1; 1 8 8]; d=A*xtrue;\n\njulia> integerfeasibility(A,d)\n3-element Array{Float64,1}:\n 0.0\n 2.0\n 9.0\n\njulia> n=20;m=30; A = rand(-10:10,n,m); xtrue = rand(0:10,m); d=A*xtrue;\n\njulia> sum(abs.(xtrue - integerfeasibility(A,d) ))\n0\n\n\n\n\n\n\n","category":"function"},{"location":"functions/#LLLplus.rationalapprox","page":"Functions","title":"LLLplus.rationalapprox","text":"rationalapprox(x::AbstractArray{<:Real,1},M,Ti=BigInt,verbose=false)\n\nFor a vector of Reals x, and an integer M, find an integer q such that maximum(abs.(x*q-round.(x*q))) is small; the vector x is approximated by round.(x*q)//q.  The integer q is less than or equal to M and the approximation satisfies max(abs.(x*q-round.(x*q)))≤sqrt(5)*2^(n/4 - 5)*M^(-1/n); this equation comes from the paper below.  The LLL algorithm reduction is used to find the solution. The approximation vector is returned. This is also known as \"simultaneous diophantine approximation\"; see for example the title of the Hanrot paper below.\n\nThis is not a robust tool, just a demo.\n\n\"LLL: A Tool for Effective Diophantine Approximation\" by Guillaume Hanrot in the book \"The LLL Algorithm: Survey and Applications\" edited by Phong Q. Nguyen and Brigitte Vallée, Springer, Heidelberg, 2010.\n\nSee also Chapter 9 of M. R. Bremner, \"Lattice Basis Reduction: An Introduction to the LLL Algorithm and Its Applications\" CRC Press, 2012.\n\nExamples\n\njulia> x = [0.3912641745333527; 0.5455179974014548; 0.1908698210882469];\n\njulia> rationalapprox(x,1e4,Int64)\n3-element Array{Rational{Int64},1}:\n 43//110\n  6//11\n 21//110\n\n\n\n\n\n\n","category":"function"},{"location":"functions/#LLLplus.spigotBBP","page":"Functions","title":"LLLplus.spigotBBP","text":"spigotBBP(α::Td,s,b,n,K,verbose=false) where {Td}\n\nCheck for a BBP-style [1] infinite series for the constant α.  These are \"spigot\" formulas that can be used to generate (for example) the millionth digit of the constant α without learning the previous digits. Specifically, given the constant α, and parameters b, n, and s, look for a vector of numbers a_1 through a_n that satisfies the following equation:\n\nalpha= sum_k=0^infty frac1b^k left( fraca_1(nk+1)^s + ldots + fraca_n(nk+n)^s right)\n\nBecause it's hard to sum to infinity, the sum is stopped at K. If a formula is found, it is printed to the screen in LaTeX and the coefficents a are returned as a vector.  An online LaTeX viewer such as https://www.latex4technics.com/ may be helpful.\n\nThis is not a robust tool, just a demo. For example, there may be a  problem with s≥2. See [2] for derivation of the technique used, and to  check whether a formula you find is new.\n\n[1] David Bailey, Peter Borwein, and Simon Plouffe. \"On the rapid computation of various polylogarithmic constants.\" Mathematics of Computation 66.218 (1997): 903-913. https://www.ams.org/journals/mcom/1997-66-218/S0025-5718-97-00856-9/\n\n[2] David Bailey, \"A Compendium of BBP-Type Formulas for Mathematical Constants\". https://www.davidhbailey.com//dhbpapers/bbp-formulas.pdf\n\nExample\n\njulia> spigotBBP(BigFloat(pi),1,16,8,45,true);\nA solution was found w error -4.728672e-60. In LaTeX form it is\n\\alpha= \\sum_{k=0}^\\infty \\frac{1}{16^k} \\left(\\frac{4}{8k+1}-\\frac{2}{8k+4}-\\frac{1}{8k+5}-\\frac{1}{8k+6}\\right)\n\nOther examples without output:\n\nspigotBBP(Float64(pi),1,-4,4,22,true);\nspigotBBP(log(2),1,2,2,30,true);\nspigotBBP(9*log(3),1,9,2,30,true);\nspigotBBP(atan(2)*8,1,16,8,30,true);\nspigotBBP(8*sqrt(2)*log(1+sqrt(2)),1,16,8,25,true);\n\nThere is a formula for pi^2 which the following command should find, but it does not find it. In fact the technique doesn't seem to work at all for  s>2; It's not obvious what the problem is\n\nspigotBBP(BigFloat(pi)*pi,2,64,6,25,true);\n\n\n\n\n\n","category":"function"},{"location":"functions/#LLLplus.hard_sphere","page":"Functions","title":"LLLplus.hard_sphere","text":"hard_sphere(...) = hardsphere(...)\n\nSee hardsphere; in a future version hard_sphere will be removed.\n\n\n\n\n\n","category":"function"},{"location":"functions/#LLLplus.issizereduced","page":"Functions","title":"LLLplus.issizereduced","text":"issizereduced(B)\n\nDetermine if the matrix B is size reduced or not.\n\nExamples\n\njulia> H= [1 2; 3 4];issizereduced(H)\nfalse\n\njulia> B,T = lll(H);issizereduced(B)\ntrue\n\n\n\n\n\n\n","category":"function"},{"location":"functions/#LLLplus.islllreduced","page":"Functions","title":"LLLplus.islllreduced","text":"islllreduced(B)\n\nDetermine if the matrix B is LLL reduced or not. See p 56 of Bremner for a definition. \n\nM. R. Bremner, \"Lattice Basis Reduction: An Introduction to the LLL  Algorithm and Its Applications\" CRC Press, 2012.\n\nExamples\n\njulia> H= [1 2; 3 4];islllreduced(H)\nfalse\n\njulia> B,T=lll(H);islllreduced(B)\ntrue\n\n\n\n\n\n\n","category":"function"},{"location":"functions/#LLLplus.orthogonalitydefect","page":"Functions","title":"LLLplus.orthogonalitydefect","text":"orthogonalitydefect(B)\n\nFind the orthogonality defect of the matrix B defined, for example, on page 2 of Bennet \n\nBennet\n\nExamples\n\njulia> H= [1 2; 3 4];B,T=lll(H);\n\njulia> [orthogonalitydefect(H) orthogonalitydefect(B)]\n1×2 Array{Float64,2}:\n 7.07107  1.0\n\n\n\n\n\n\n","category":"function"},{"location":"functions/#LLLplus.hermitefactor","page":"Functions","title":"LLLplus.hermitefactor","text":"hermitefactor(B)\n\nFind the Hermite factor of matrix B\n\nExamples\n\njulia> H= [1 2; 3 4];hermitefactor(H)\n1.5811388300841898\n\n\n\n\n\n\n","category":"function"},{"location":"functions/#LLLplus.seysencond","page":"Functions","title":"LLLplus.seysencond","text":"seysencond(B)\n\nSeysen condition number as on, for example, page 3 of Bennet \n\nBennet\n\nExamples\n\njulia> H= [1 2; 3 4];seysencond(H)\n2.8284271247461903\n\n\n\n\n\n\n","category":"function"},{"location":"functions/#LLLplus.gen_qary_b","page":"Functions","title":"LLLplus.gen_qary_b","text":"b= gen_qary_b(T, d::Int,k::Int,b::Int)\n\nGenerate a q-ary lattice given an element type T, dimension d, parameter k, and bit-depth b. Specifically, find a d by d matrix which has the block structure [I zeros(T,k,d-k); H q*I]], where the k by d-k matrix H is sampled from 0:q-1 and q is sampled uniformly from 1:big(2)^b-1\n\nThese bases correspond to the SIS/LWE q-ary lattices; see D. Micciancio and O. Regev. Post-Quantum Cryptography. Chapter of Lattice-based Cryptography, 147-191 (2009) and latticegen in https://github.com/fplll/fplll\n\nExamples\n\njulia> b=gen_qary_b(Int64,2,1,6)\n2×2 Array{Int64,2}:\n 1   0\n 7  32\n\n\n\n\n\n\n","category":"function"},{"location":"#LLLplus-README","page":"Home","title":"LLLplus README","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = LLLplus","category":"page"},{"location":"","page":"Home","title":"Home","text":"LLLplus includes Lenstra-Lenstra-Lovász (LLL), Brun, and Seysen lattice reduction; and shortest vector problem (SVP) and closest vector problem (CVP) solvers. These lattice reduction and related lattice tools are used in cryptography, digital communication, and integer programming. The historical and practical prominence of the LLL technique in lattice tools is the reason for its use in the name \"LLLplus\". This package is experimental; see fplll for a robust tool.","category":"page"},{"location":"","page":"Home","title":"Home","text":"LLL [1] lattice reduction is a powerful tool that is widely used in cryptanalysis, in cryptographic system design, in digital communications, and to solve other integer problems.  LLL reduction is often used as an approximate solution to the SVP. We also include Gauss/Lagrange, Brun [2] and Seysen [3] lattice reduction techniques. The LLL, Brun, and Seysen algorithms are based on [4]. The CVP solver is based on [5] and can handle lattices and bounded integer constellations. A slow SVP solver based on the CVP tool is included as well.","category":"page"},{"location":"","page":"Home","title":"Home","text":"We also include code to do a Vertical-Bell Laboratories Layered Space-Time (V-BLAST) [6] matrix decomposition which is used in digital communications. The LLL, Brun, Seysen, V-BLAST, and CVP functions are used to solve (exactly or approximately) CVP problems; the MUMIMO.jl package demostrates how these functions can be used in encoding and decoding multi-antenna signals.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Another important application is in cryptanalysis; as an example of a cryptanalytic attack, see the subsetsum function.  The LLL algorithm has been shown to solve the integer programming feasibility problem; see integerfeasibility. Lattice tools are often used to study and solve Diophantine problems; for example in  \"simultaneous diophantine approximation\" a vector of real numbers are approximated by rationals with a common deonminator. For a demo function, see rationalapprox. Finally, to see how the LLL can be used to find spigot formulas for irrationals, see spigotBBP.","category":"page"},{"location":"#Examples","page":"Home","title":"Examples","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Each function contains documentation and examples available via Julia's built-in documentation system, for example with ?lll. Documentation for all functions is available on pkg.julialang.org. A tutorial notebook is found in the docs directory or on nbviewer.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Here are a few examples of using the functions in the package on random lattices.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Pkg.add(\"LLLplus\")\nusing LLLplus\n# repeat the commands below to remove JIT compile time\n\n# Time the decomposition of a matrix with randn entries\nN = 100;\nH = randn(N,N);\n@time B,T = sizereduction(H);\n@time B,T = brun(H);\n@time B,T = lll(H);\n@time B,T = seysen(H);\n@time W,P,B = vblast(H);\n\n# check out the CVP solver\n@time Q,R=qr(H);\nu=Int.(rand(0:1e10,N));\ny=H*u+rand(N)/100;\n@time uhat=cvp(Q'*y,R);\nsum(abs.(u-uhat))","category":"page"},{"location":"#Execution-Time-results","page":"Home","title":"Execution Time results","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"In the first test we compare the lll function from LLLplus, the l2avx function in the src\\l2.jl file in LLLplus, the lll_with_transform function from Nemo (which uses FLINT), and the lll_reduction function from fplll. Nemo and fplll are written by number theorists and are good benchmarks against which to compare.  We first show how the execution time varies as the basis (matrix) size varies over [4 8 16 32 64]. For each matrix size, 20 random bases are generated using fplll's gen_qary function with depth of 25 bits, with the average execution time shown; the eltype is Int64 except for NEMO, which uses GMP (its own BigInt); in all cases the δ=.99. The vertical axis shows execution time on a logarithmic scale; the x-axis is also logarithmic. The generally linear nature of the LLL curves supports the polynomial-time nature of the algorithm. The LLLplus.lll function is slower, while l2avx is similar to fplll. Though not shown, using bases from gen_qary with bit depth of 45 gives fplll a larger advantage. This figure was generated using code in test/timeLLLs.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: Time vs basis size)","category":"page"},{"location":"","page":"Home","title":"Home","text":"One question that could arise when looking at the plot above is what the quality of the basis is. In the next plot we show execution time vs the norm of the first vector in the reduced basis, this first vector is typically the smallest; its norm is an rough indication of the quality of the reduced basis. We show results averaged over 20 random bases from gen_qary with depth 25 bits, this time with the dimension fixed at 32. The curve is created by varying the δ parameter from .29 to .99 in steps of .2; the larger times and smaller norms correspond to the largest δ values. Though the l2avx function is competitive with fplll in this case, in many other cases the fplll code is significantly faster.","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: Time vs reduction quality)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Finally, we show execution time for several built-in datatypes (Int32, Int64, Int128, Float32, Float64, BitInt, and BigFloat) as well as type from external packages (Float128 from Quadmath.jl and Double64 from DoubleFloat.jl) which are used to generate 40 128x128 matrices, over which execution time for the lattice reduction techniques is averaged.  The vertical axis is a logarithmic representation of execution time as in the previous figure. This figure was generated using code in test/perftest.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: Time vs data type)","category":"page"},{"location":"#Notes","page":"Home","title":"Notes","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"There are certainly many improvements and additions that could be made to LLLplus, such as adding Block-Korkin-Zolotarev (BKZ) lattice reduction with improvements as in [8]. Even so, it would be hard to compete with fplll on features. In fact, a Julia wrapper around fplll would be the most useful addition to lattice tools in Julia; it would provide funcionality not in LLLplus such as BKZ reduction.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The algorithm pseudocode in the monograph [7] and the survey paper [4] were very helpful in writing the lattice reduction tools in LLLplus and are a good resource for further study. If you are trying to break one of the Lattice Challenge records or are looking for robust, well-proven lattice tools, look at fplll. Also, for many number-theoretic problems the Nemo.jl package is appropriate; it uses the FLINT C library to do LLL reduction on Nemo-specific data types.  Finally, no number theorists have worked on LLLplus; please treat the package as experimental.","category":"page"},{"location":"#References","page":"Home","title":"References","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"[1] A. K. Lenstra; H. W. Lenstra Jr.; L. Lovász, \"Factoring polynomials with rational coefficients\". Mathematische Annalen 261, 1982.","category":"page"},{"location":"","page":"Home","title":"Home","text":"[2] V. Brun, \"En generalisation av kjedebrøken I,\" Skr. Vidensk. Selsk. Kristiana, Mat. Nat. Klasse, 1919.","category":"page"},{"location":"","page":"Home","title":"Home","text":"[3] M. Seysen, \"Simultaneous reduction of a lattice basis and its reciprocal basis\" Combinatorica, 1993.","category":"page"},{"location":"","page":"Home","title":"Home","text":"[4] D. Wuebben, D. Seethaler, J. Jalden, and G. Matz, \"Lattice Reduction - A Survey with Applications in Wireless Communications\". IEEE Signal Processing Magazine, 2011.","category":"page"},{"location":"","page":"Home","title":"Home","text":"[5] A. Ghasemmehdi, E. Agrell, \"Faster Recursions in Sphere Decoding\" IEEE Transactions on Information Theory, vol 57, issue 6 , June 2011.","category":"page"},{"location":"","page":"Home","title":"Home","text":"[6] P. W. Wolniansky, G. J. Foschini, G. D. Golden, R. A. Valenzuela, \"V-BLAST: An Architecture for Realizing Very High Data Rates Over the Rich-Scattering Wireless Channel\". Proc. URSI ISSSE: 295–300, 1998. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"[7] M. R. Bremner, \"Lattice Basis Reduction: An Introduction to the LLL  Algorithm and Its Applications\" CRC Press, 2012.","category":"page"},{"location":"","page":"Home","title":"Home","text":"[8] Y. Chen, P. Q. Nguyen, \"BKZ 2.0: Better Lattice Security Estimates\". Proc. ASIACRYPT 2011.","category":"page"},{"location":"#List-of-Functions","page":"Home","title":"List of Functions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"}]
}