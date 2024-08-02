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

ψ = qnmfunctionnew(-2,2,2,0,0.5)

# testop=1.0+0.5im
# testom=1.0+0.5im
# testhp=1.0+0.5im
# testhm=1.0+0.5im
# testip=1.0+0.5im
# testim=1.0+0.5im

# @show Computeω2(testop,testom,testhp,testhm,testip,testim,ψ)

weight = let s = ψ.s , a= ψ.a
    (r,z) -> sqrt(1-z^2)*(r^2+a^2-2*r)^s
end

∂ωOplus = OperatorShift(file1)
# ∂ωOminus= OperatorShift(file2)
# Hplus = OperatorShift(file3)
# Hminus = OperatorShift(file4)
# Iplus = OperatorShift(file5)
# Iminus = OperatorShift(file6)

∂ωOpluss = OperatorSandwich(ψ,∂ωOplus,weight,ψ)
∂ωOplus = ∂ωOpluss.Op;

# ∂ωOminuss = OperatorSandwich(ψ,∂ωOminus,weight,ψ)
# ∂ωOminus = ∂ωOminuss.Op;

# Hpluss = OperatorSandwich(ψ,Hplus,weight,ψ)
# Hplus = Hpluss.Op;

# Hminuss = OperatorSandwich(ψ,Hminus,weight,ψ)
# Hminus = Hminuss.Op;

# Ipluss = OperatorSandwich(ψ,Iplus,weight,ψ)
# Iplus = Ipluss.Op;

# Iminuss = OperatorSandwich(ψ,Iminus,weight,ψ)
# Iminus = Iminuss.Op;

@show ∂ωOplus(1,0)
@show ∂ωOplus(1,0;LHSisconjugate=false)
@show ∂ωOplus(1,0;LHSisconjugate=true)
# @show ∂ωOminus(1,0)
# ∂ωOminus(1,0,isconjugate=true,isminus=true) |> println
# @show typeof(∂ωOminus)
# @show typeof(Hplus)
# @show Hplus(1,0,isconjugate=true)
# @show Hminus(1,0,isconjugate=false)
# @show Iplus(1,0,isconjugate=false)
# @show Iminus(1,0,isconjugate=true)

### Define the useful contours
r₊ = ψ.R.r₊ ; r₋ = ψ.R.r₋ ; s = ψ.s ; Δr = 0.1*(r₊-r₋); ϵ = eps(0.1);

point1 = r₊ + Δr - Δr*im
point2 = r₊ - Δr - Δr*im

radial1 = SemiInfiniteLine(point1 , point1 + Δr*im , false)
angular = LineSegment(-1.0+100*ϵ , 1.0-100*ϵ , true) #to avoid the NaNs at the edges
C1 = radial1 ⊗ angular

radial2 = LineSegment(point1,point2,true)
C2 = radial2 ⊗ angular

radial3 = SemiInfiniteLine(point2 , point2 + Δr*im , true)
C3 = radial3 ⊗ angular

TheContour = C1⊕C2⊕C3

println("Just integrals left")

@show Integrate(∂ωOplus, TheContour)
# @show Integrate(∂ωOminus, TheContour)
# @show Integrate(Hplus, TheContour)
# @show Integrate(Hminus, TheContour)
# @show Integrate(Iplus, TheContour)
# @show Integrate(Iminus, TheContour)

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

# ### Define Transformed Functions (that live on the compactified domain [0,1]⊗[0,1])
# δTₒ = TransformIntegrand(δT , thedomain)
# ∂ωTₒ = TransformIntegrand(∂ωT , thedomain)

# ### Compute the expressions so that they compile the first time
# δTₒ(0.1,0.1) |> println
# ∂ωTₒ(0.1,0.1) |> println

# #=
# ### Use expressions for the operators, noted in their respective csvs
# file1b = "OperatorShifts/KerrNewman0/TuekolskyShifts.csv"
# file2b = "OperatorShifts/KerrNewman0/TuekolskyFrequencyDerivative.csv"
# Tb = OperatorShift(file1b,file2b)

# ### Get the integrand in the contour integral for the inner product computation
# Tsb = OperatorSandwich(ψ,Tb,weight,ψ)
# δTb = Tsb.δT;
# ∂ωTb = Tsb.∂ωT;

# ### Compute the expressions so that they compile the first time
# δTb(2.2,0.1) |> println
# ∂ωTb(2.2,0.1) |> println

# ### Define Transformed Functions (that live on the compactified domain [0,1]⊗[0,1])
# δTbₒ = TransformIntegrand(δTb , thedomain)
# ∂ωTbₒ = TransformIntegrand(∂ωTb , thedomain)

# ### Compute the expressions so that they compile the first time
# δTbₒ(0.1,0.1) |> println
# ∂ωTbₒ(0.1,0.1) |> println
# =#
