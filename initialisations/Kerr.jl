# Add the path to your local package to the LOAD_PATH
push!(LOAD_PATH, "C:/Users/dwuuu/Documents/GitHub/ContourIntegrals.jl")
push!(LOAD_PATH, "C:/Users/dwuuu/Documents/GitHub/SchwarzschildQNMs")
push!(LOAD_PATH, "C:/Users/dwuuu/Documents/GitHub/KerrQuasinormalModes.jl")

using ContourIntegrals
using KerrQNMShifts


### Use expressions for the operators, noted in their respective csvs
file1 = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/SchwarzschildQNMShifts/OperatorShifts/Schwarzschild/dwOpluscoefficients.csv"
file2 = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/SchwarzschildQNMShifts/OperatorShifts/Schwarzschild/dwOminuscoefficients.csv"
file3 = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/SchwarzschildQNMShifts/OperatorShifts/Schwarzschild/Hpluscoefficients.csv"
file4 = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/SchwarzschildQNMShifts/OperatorShifts/Schwarzschild/Hminuscoefficients.csv"
file5 = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/SchwarzschildQNMShifts/OperatorShifts/Schwarzschild/Ipluscoefficients.csv"
file6 = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/SchwarzschildQNMShifts/OperatorShifts/Schwarzschild/Iminuscoefficients.csv"

∂ωOplus = OperatorShift(file1)
∂ωOminus= OperatorShift(file2)
Hplus = OperatorShift(file3)
Hminus = OperatorShift(file4)
Iplus = OperatorShift(file5)
Iminus = OperatorShift(file6)

# ψ = qnmfunctionnew(-2,2,2,0,0.5)
# ψ2 = qnmfunctionnew(-2,2,2,0,0.0)
# ψ3 = qnmfunctionnew(-2,2,2,0,0.0000001)

ψm2even=qnmfunctionnew(-2,2,2,0,0.5)
ψp2even=qnmfunctionnew(2,2,2,0,0.5)
ψm2odd=qnmfunctionnew(-2,3,2,0,0.5)
ψp2odd=qnmfunctionnew(2,3,2,0,0.5)
ψm2oddalt=qnmfunctionnew(-2,5,3,0,0.5)
ψp2oddalt=qnmfunctionnew(2,5,3,0,0.5)

println("Check ψm2even conjugates:")
@show ψm2even.R(0.5)
@show ψm2even.R(0.5,isconjugate=true)
@show ψm2even.S(0.2)
@show ψm2even.S(0.2,isconjugate=true)
@show ψm2even(0.5,0.2)
@show ψm2even(0.5,0.2,isconjugate=true)
println("Check ψm2even isminus:")
@show ψm2even.R(0.5)
@show ψm2even.R(0.5,isminus=true)
@show ψp2even.S(0.2)
@show ψm2even.S(0.2,isminus=true)
@show ψm2even(0.5,.2)
@show ψm2even(0.5,.2,isminus=true)
println("Check ψm2oddalt isminus:")
@show ψm2oddalt.R(0.5)
@show ψm2oddalt.R(0.5,isminus=true)
@show ψm2oddalt.S(0.2)
@show ψm2oddalt.S(0.2,isminus=true)
@show ψm2oddalt(0.5,.2)
@show ψm2oddalt(0.5,.2,isminus=true)


# weight = let s = ψ.s , a= ψ.a
#     (r,z) -> sqrt(1-z^2)*(r^2+a^2-2*r)^s
# end

# ∂ωOpluss1 = OperatorSandwich(ψ,∂ωOplus,weight,ψ)
# ∂ωOplus1 = ∂ωOpluss1.Op;

# ∂ωOminuss1 = OperatorSandwich(ψ,∂ωOminus,weight,ψ)
# ∂ωOminus1 = ∂ωOminuss1.Op;

# Hpluss1 = OperatorSandwich(ψ,Hplus,weight,ψ)
# Hplus1 = Hpluss1.Op;

# Hminuss1 = OperatorSandwich(ψ,Hminus,weight,ψ)
# Hminus1 = Hminuss1.Op;

# Ipluss1 = OperatorSandwich(ψ,Iplus,weight,ψ)
# Iplus1 = Ipluss1.Op;

# Iminuss1 = OperatorSandwich(ψ,Iminus,weight,ψ)
# Iminus1 = Iminuss1.Op;


# ∂ωOpluss2 = OperatorSandwich(ψ2,∂ωOplus,weight,ψ2)
# ∂ωOplus2 = ∂ωOpluss2.Op;

# ∂ωOminuss2 = OperatorSandwich(ψ2,∂ωOminus,weight,ψ2)
# ∂ωOminus2 = ∂ωOminuss2.Op;

# Hpluss2 = OperatorSandwich(ψ2,Hplus,weight,ψ2)
# Hplus2 = Hpluss2.Op;

# Hminuss2 = OperatorSandwich(ψ2,Hminus,weight,ψ2)
# Hminus2 = Hminuss2.Op;

# Ipluss2 = OperatorSandwich(ψ2,Iplus,weight,ψ2)
# Iplus2 = Ipluss2.Op;

# Iminuss2 = OperatorSandwich(ψ2,Iminus,weight,ψ2)
# Iminus2 = Iminuss2.Op;


# ∂ωOpluss3 = OperatorSandwich(ψ3,∂ωOplus,weight,ψ3)
# ∂ωOplus3 = ∂ωOpluss3.Op;

# ∂ωOminuss3 = OperatorSandwich(ψ3,∂ωOminus,weight,ψ3)
# ∂ωOminus3 = ∂ωOminuss3.Op;

# Hpluss3 = OperatorSandwich(ψ3,Hplus,weight,ψ3)
# Hplus3 = Hpluss3.Op;

# Hminuss3 = OperatorSandwich(ψ3,Hminus,weight,ψ3)
# Hminus3 = Hminuss3.Op;

# Ipluss3 = OperatorSandwich(ψ3,Iplus,weight,ψ3)
# Iplus3 = Ipluss3.Op;

# Iminuss3 = OperatorSandwich(ψ3,Iminus,weight,ψ3)
# Iminus3 = Iminuss3.Op;

# @show ∂ωOplus1(1,0)
# @show ∂ωOminus1(1,0)


# @show Hplus1(0.5,.2)
# @show Hplus1(0.5,.2,isconjugate=true)
# @show Hplus1(0.5,.2,isminus=true)
# @show Hplus1(0.5,.2,LHSisminus=true)
# @show Hplus1(0.5,.2,isconjugate=true, isminus=true)
# @show Hplus1(0.5,.2,isconjugate=true,isminus=true,LHSisminus=true)
# @show Hplus1(0.5,.2,LHSisminus=true)

# @show Hminus1(1,0)
# @show Iplus1(1,0)
# @show Iminus1(1,0)

# @show ∂ωOplus2(1,0)
# @show ∂ωOminus2(1,0)
# @show Hplus2(1,0)
# @show Hminus2(1,0)
# @show Iplus2(1,0)
# @show Iminus2(1,0)

# @show ∂ωOplus3(1,0)
# @show ∂ωOminus3(1,0)
# @show Hplus3(1,0)
# @show Hminus3(1,0)
# @show Iplus3(1,0)
# @show Iminus3(1,0)

# ### Define the useful contours
# r₊ = ψ.R.r₊ ; r₋ = ψ.R.r₋ ; s = ψ.s ; Δr = 0.1*(r₊-r₋); ϵ = eps(0.1);

# point1 = r₊ + Δr - Δr*im
# point2 = r₊ - Δr - Δr*im

# radial1 = SemiInfiniteLine(point1 , point1 + Δr*im , false)
# angular = LineSegment(-1.0+100*ϵ , 1.0-100*ϵ , true) #to avoid the NaNs at the edges
# C1 = radial1 ⊗ angular

# radial2 = LineSegment(point1,point2,true)
# C2 = radial2 ⊗ angular

# radial3 = SemiInfiniteLine(point2 , point2 + Δr*im , true)
# C3 = radial3 ⊗ angular

# TheContour = C1⊕C2⊕C3

# println("Just integrals left")

# 𝒪plus= Integrate(∂ωOplus, TheContour)[1]
# 𝒪minus= conj(Integrate(∂ωOminus, TheContour,isminus=true,LHSisminus=true)[1])
# println("Done O's")
# ℋplus= Integrate(Hplus, TheContour)[1]
# println("Done Hplus")
# ℋminus= conj(Integrate(Hminus, TheContour,isminus=true,LHSisminus=true)[1])
# println("Done Hminus")
# ℐplus= Integrate(Iplus, TheContour,isconjugate=true,isminus=true)[1]
# println("Done Iplus")
# ℐminus= conj(Integrate(Iminus, TheContour,isconjugate=true,LHSisminus=true)[1])
# println("Done Iminus")

# ω2shift=Computeω2(𝒪plus,𝒪minus,ℋplus,ℋminus,ℐplus,ℐminus,ψ)
# @show ω2shift

# 𝒪plus2= Integrate(∂ωOplus2, TheContour)[1]
# 𝒪minus2= conj(Integrate(∂ωOminus2, TheContour,isminus=true,LHSisminus=true)[1])
# println("Done O's")
# ℋplus2= Integrate(Hplus2, TheContour)[1]
# println("Done Hplus")
# ℋminus2= conj(Integrate(Hminus2, TheContour,isminus=true,LHSisminus=true)[1])
# println("Done Hminus")
# ℐplus2= Integrate(Iplus2, TheContour,isconjugate=true,isminus=true)[1]
# println("Done Iplus")
# ℐminus2= conj(Integrate(Iminus2, TheContour,isconjugate=true,LHSisminus=true)[1])
# println("Done Iminus")

# 𝒪plus3=-0.11044200242553877 + 0.49269650825808775im
# 𝒪minus3=-0.001844270788297432 + 0.042774974720934486im
# ℋplus3=-2.445285430275201e-7 + 1.0519871516885872e-7im
# ℋminus3=1.1008476910844426e-7 + 2.925268041340887e-7im
# ℐplus3=3.274068559342357e-6 + 9.925224697278246e-6im
# ℐminus3=2.2486149550283164e-6 - 3.7324257307655107e-6im


# ω2shift3=Computeω2(𝒪plus3,𝒪minus3,ℋplus3,ℋminus3,ℐplus3,ℐminus3,ψ3)
# @show ω2shift3


# 𝒪plus3= Integrate(∂ωOplus3, TheContour)[1]
# 𝒪minus3= conj(Integrate(∂ωOminus3, TheContour,isminus=true,LHSisminus=true)[1])
# println("Done O's")
# ℋplus3= Integrate(Hplus3, TheContour)[1]
# println("Done Hplus")
# ℋminus3= conj(Integrate(Hminus3, TheContour,isminus=true,LHSisminus=true)[1])
# println("Done Hminus")
# ℐplus3= Integrate(Iplus3, TheContour,isconjugate=true,isminus=true)[1]
# println("Done Iplus")
# ℐminus3= conj(Integrate(Iminus3, TheContour,isconjugate=true,LHSisminus=true)[1])
# println("Done Iminus")

# println("R's")
# println(ψ.R(1))
# println(ψ.R(1;isconjugate=true))
# println(ψ.R(1;isminus=true))
# println(ψ.R(1;isconjugate=true,isminus=true))


# println("S's")
# println(ψ.S(0))
# println(ψ.S(0,isconjugate=true))
# println(ψ.S(0,isminus=true))

# println("Full funcs")
# println(ψ(1,0))
# println(ψ(1,0,isconjugate=true))
# println(ψ(1,0,isminus=true))
# println(ψ(1,0,isconjugate=true,isminus=true))

# println("ψ2.S")
# ψ2 = qnmfunctionnew(2,2,2,0,0.5)
# println("ψ2.S(0): ", ψ2.S(0))

# ∂ωOminus=OperatorShift(file2)
# ∂ωOplus=OperatorShift(file1)

# weight = let s = ψ.s , a= ψ.a
#     (r,z) -> sqrt(1-z^2)*(r^2+a^2-2*r)^s
# end

# ∂ωOminuss = OperatorSandwich(ψ,∂ωOminus,weight,ψ)
# ∂ωOminus = ∂ωOminuss.Op;

# ∂ωOpluss = OperatorSandwich(ψ,∂ωOplus,weight,ψ)
# ∂ωOplus = ∂ωOpluss.Op;

# @show ∂ωOminus(1.0,0.0)

# ### Define the useful contours
# r₊ = ψ.R.r₊ ; r₋ = ψ.R.r₋ ; s = ψ.s ; Δr = 0.1*(r₊-r₋); ϵ = eps(0.1);

# point1 = r₊ + Δr - Δr*im
# point2 = r₊ - Δr - Δr*im

# radial1 = SemiInfiniteLine(point1 , point1 + Δr*im , false)
# angular = LineSegment(-1.0+100*ϵ , 1.0-100*ϵ , true) #to avoid the NaNs at the edges
# C1 = radial1 ⊗ angular

# radial2 = LineSegment(point1,point2,true)
# C2 = radial2 ⊗ angular

# radial3 = SemiInfiniteLine(point2 , point2 + Δr*im , true)
# C3 = radial3 ⊗ angular

# TheContour = C1⊕C2⊕C3

# println(Integrate(∂ωOminus, TheContour))

# qnmfunctionnew(s,l,m,n,a; qnm=qnm)



# weight = let s = ψ.s , a= ψ.a
#     (r,z) -> sqrt(1-z^2)*(r^2+a^2-2*r)^s
# end

# println("omega:", ψ.ω)

# println("weight(1,0):", weight(1,0))
# println("ψ00(1,0)=",ψ(1,0))

# ∂ωOminuss = OperatorSandwich(ψ,∂ωOminus,weight,ψ)
# ∂ωOminus = ∂ωOminuss.Op;

# ∂ωOminus(1,0,isconjugate=true) |> println
# ∂ωOminus(1,0) |> println 
# ∂ωOminus(1,0,isconjugate=false) |> println 


# testfunc=qnmfunctionnew(-2,2,2,0,0.5)
# println(testfunc(0))
# println(testfunc(0, conjugate=true))



# # The below function prints out what each operator is from both file 1 and file2
# ∂ωOplus = OperatorShift(file1)
# ∂ωOminus= OperatorShift(file2)
# Hplus = OperatorShift(file3)
# Hminus = OperatorShift(file4)
# Iplus = OperatorShift(file5)
# Iminus = OperatorShift(file6)


# ### Get the integrand in the contour integral for the inner product computation
# # qnmfunctionnew(s,l,m,n,a, plusminus);
# ψplus = qnmfunctionnew(-2,2,2,0,0.5,"plus")
# ψminus = qnmfunctionnew(-2,2,2,0,0.5,"minus")

# weight = let s = ψplus.s , a= ψplus.a
#     (r,z) -> sqrt(1-z^2)*(r^2+a^2-2*r)^s
# end

# ∂ωOpluss = OperatorSandwich(ψplus,∂ωOplus,weight,ψplus)
# ∂ωOplus = ∂ωOpluss.Op;

# ∂ωOminuss = OperatorSandwich(ψminus,∂ωOminus,weight,ψminus)
# ∂ωOminus = ∂ωOminuss.Op;

# Hpluss = OperatorSandwich(ψplus,Hplus,weight,ψplus)
# Hplus = Hpluss.Op;

# Hminuss = OperatorSandwich(ψminus,Hminus,weight,ψminus)
# Hminus = Hminuss.Op;

# Ipluss = OperatorSandwich(ψplus,Iplus,weight,ψminus)
# Iplus = Ipluss.Op;

# Iminuss = OperatorSandwich(ψminus,Iminus,weight,ψplus)
# Iminus = Iminuss.Op;

# println("Second Print")

# ### Compute the expressions so that they compile the first time
# ∂ωOplus(2.2,0.1) |> println
# ∂ωOminus(2.2,0.1) |> println
# Hplus(2.2,0.1) |> println
# Hminus(2.2,0.1) |> println
# Iplus(2.2,0.1) |> println
# Iminus(2.2,0.1) |> println

# ### Define the useful contours
# r₊ = ψ.R.r₊ ; r₋ = ψ.R.r₋ ; s = ψ.s ; Δr = 0.1*(r₊-r₋); ϵ = eps(0.1);

# point1 = r₊ + Δr - Δr*im
# point2 = r₊ - Δr - Δr*im

# radial1 = SemiInfiniteLine(point1 , point1 + Δr*im , false)
# angular = LineSegment(-1.0+100*ϵ , 1.0-100*ϵ , true) #to avoid the NaNs at the edges
# C1 = radial1 ⊗ angular

# radial2 = LineSegment(point1,point2,true)
# C2 = radial2 ⊗ angular

# radial3 = SemiInfiniteLine(point2 , point2 + Δr*im , true)
# C3 = radial3 ⊗ angular

# TheContour = C1⊕C2⊕C3
