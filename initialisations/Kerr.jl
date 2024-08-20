# Add the path to your local package to the LOAD_PATH
push!(LOAD_PATH, "C:/Users/dwuuu/Documents/GitHub/ContourIntegrals.jl")
push!(LOAD_PATH, "C:/Users/dwuuu/Documents/GitHub/SchwarzschildQNMs")
push!(LOAD_PATH, "C:/Users/dwuuu/Documents/GitHub/KerrQuasinormalModes.jl")

using ContourIntegrals
using KerrQNMShifts
using KerrQuasinormalModes


### Use expressions for the operators, noted in their respective csvs
file1 = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/SchwarzschildQNMShifts/OperatorShifts/Schwarzschild/dwOpluscoefficients.csv"
file2 = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/SchwarzschildQNMShifts/OperatorShifts/Schwarzschild/dwOminuscoefficients.csv"
file3 = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/SchwarzschildQNMShifts/OperatorShifts/Schwarzschild/Hpluscoefficients.csv"
file4 = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/SchwarzschildQNMShifts/OperatorShifts/Schwarzschild/Hminuscoefficients.csv"
file5 = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/SchwarzschildQNMShifts/OperatorShifts/Schwarzschild/Ipluscoefficients.csv"
file6 = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/SchwarzschildQNMShifts/OperatorShifts/Schwarzschild/Iminuscoefficients.csv"

filetest="C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/derivtestcoefficients.csv"

testshift = OperatorShift(filetest)

∂ωOplus = OperatorShift(file1)
# ∂ωOminus= OperatorShift(file2)
# Hplus = OperatorShift(file3)
# Hminus = OperatorShift(file4)
# Iplus = OperatorShift(file5)
# Iminus = OperatorShift(file6)

ψ1 = qnmfunctionnew(-2,5,3,0,0.5)
ψ2= qnmfunctionnew(-2,5,3,0,0.5,modesign="minus")

η1=ψ1.R.η
α1=ψ1.R.α
ξ1=ψ1.R.ξ
ζ1=ψ1.R.ζ
r₊1=ψ1.R.r₊
r₋1=ψ1.R.r₋

η2=ψ2.R.η
α2=ψ2.R.α
ξ2=ψ2.R.ξ
ζ2=ψ2.R.ζ
r₊2=ψ2.R.r₊
r₋2=ψ2.R.r₋

r0=10
z0=0.2

# @show L4(ψ1,r0,z0)/ψ1.R(r0)/ψ3.S(z0)
# ComputeDplus(ψ1)

# # Check the derivative functions
# δr=.000001
# δz=.00001
# @show ∂r(ψ1.R)(r)-(ψ1.R(r+δr)-ψ1.R(r))/δr
# @show ∂r(ψ1)(r,z)-(ψ1(r+δr,z)-ψ1(r,z))/δr
# @show ∂r(∂r(ψ1))(r,z)-(ψ1(r+δr,z)-2*ψ1(r,z)+ψ1(r-δr,z))/δr^2

# @show ∂θ(ψ1.S)(z)-(ψ1.S(z+δz)-ψ1.S(z))/δz*(-sqrt(1-z^2))
# @show ∂θ(ψ1)(r,z)-(ψ1(r,z+δz)-ψ1(r,z))/δz*(-sqrt(1-z^2))
# @show ∂θ(∂θ(ψ1))(r,z)-(ψ1(r,z+δz)-2*ψ1(r,z)+ψ1(r,z-δz))/δz^2*(1-z^2)-(ψ1(r,z+δz)-ψ1(r,z))/δz*(-z)

# @show ∂r(∂θ(ψ1))(r,z)-(ψ1(r+δr,z+δz)-ψ1(r+δr,z)-ψ1(r,z+δz)+ψ1(r,z))/(δr*δz)*(-sqrt(1-z^2))



@show c1=Complex(-1)^(α1-η1-ξ1)

# function c2(r)
#     im^(α1-η1-ξ1)*exp(-(1/2)*im*(π*(-α1+η1+ξ1)+2*(-α1+η1)*angle(r-r₋1)+2*ξ1*angle(r-r₊1)))*(r-r₋1)^(α1-η1)*((r-r₋1)^2)^(1/2*(-α1+η1))*(r-r₊1)^(-ξ1)*((r-r₊1)^2)^(ξ1/2)
# end
function c2(r)
    (im*(r₋1-r))^(η1-α1)*(im*(r₊1-r))^ξ1/((im*(r-r₋1))^(η1-α1)*(im*(r-r₊1))^ξ1)
end
@show c2(r0)

@show exp(im*π*(α1-η1-3*ξ1))*(r0-r₋1)^(α1-η1)*sqrt(Complex((r0-r₋1)^(-2*α1+2*η1)))*(-r0+r₊1)^(-ξ1)*sqrt(Complex((-r0+r₊1)^(2*ξ1)))
# @show ψ1(r0)
# @show ψ1(r0,z0)

# @show ψ1.R(r)
# @show ψ2.R(r)
# @show ψ1.R(r)*c
# @show conj(ψ2.R(r))


# weight = let s = ψ1.s , a= ψ1.a
#     (r,z) -> sqrt(1-z^2)*(r^2+a^2-2*r)^s
# end

# println("Past weight")

# testshifts11=OperatorSandwich(ψ1,testshift,weight,ψ1)
# testshift11 = testshifts11.Op;

# println("Past testshift11")

# testshifts21=OperatorSandwich(ψ2,testshift,weight,ψ1)
# testshift21 = testshifts21.Op;

# testshifts12=OperatorSandwich(ψ1,testshift,weight,ψ2)
# testshift12 = testshifts12.Op;

# testshifts22=OperatorSandwich(ψ2,testshift,weight,ψ2)
# testshift22 = testshifts22.Op;

# ∂ωOpluss1 = OperatorSandwich(ψ1,∂ωOplus,weight,ψ1)
# ∂ωOplus1 = ∂ωOpluss1.Op;

# ∂ωOminuss1 = OperatorSandwich(ψ1,∂ωOminus,weight,ψ1)
# ∂ωOminus1 = ∂ωOminuss1.Op;

# Hpluss1 = OperatorSandwich(ψ1,Hplus,weight,ψ1)
# Hplus1 = Hpluss1.Op;

# Hminuss1 = OperatorSandwich(ψ1,Hminus,weight,ψ1)
# Hminus1 = Hminuss1.Op;

# Ipluss1 = OperatorSandwich(ψ1,Iplus,weight,ψ1)
# Iplus1 = Ipluss1.Op;

# Iminuss1 = OperatorSandwich(ψ1,Iminus,weight,ψ1)
# Iminus1 = Iminuss1.Op;

# @show ψ1(0.5,.2)
# @show ψ1(0.5,.2,isconjugate=true)
# @show ψ1(0.5,.2,isminus=true)
# @show ψ1(0.5,.2,isminus=true,isconjugate=true)

# @show testshift11(r,z)
# @show testshift11(r,z,isconjugate=true)
# @show testshift21(r,z,isconjugate=true)
# @show testshift21(r,z)
# @show testshift12(r,z)
# @show testshift12(r,z,isconjugate=true)
# @show testshift22(r,z)
# @show testshift22(r,z,isconjugate=true)





# @show ∂ωOplus1(r,z)
# @show ∂ωOplus1(0.5,.2,isconjugate=true)
# @show ∂ωOplus1(0.5,.2,isminus=true)
# @show ∂ωOplus1(0.5,.2,LHSisminus=true)
# @show ∂ωOplus1(0.5,.2,isconjugate=true,isminus=true)
# @show ∂ωOplus1(0.5,.2,isconjugate=true,LHSisminus=true)
# @show ∂ωOplus1(0.5,.2,isminus=true,LHSisminus=true)



# @show ∂ωOminus1(0.5,.2,LHSisminus=true)
# @show Hplus1(0.5,.2)
# @show Hplus1(0.5,.2,isconjugate=true)
# @show Hplus1(0.5,.2,isminus=true)
# @show Hplus1(0.5,.2,LHSisminus=true)
# @show Hplus1(0.5,.2,isconjugate=true, isminus=true)
# @show Hplus1(0.5,.2,isconjugate=true,isminus=true,LHSisminus=true)
# @show Hplus1(0.5,.2,LHSisminus=true)


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


# ∂ωOminus=OperatorShift(file2)
# ∂ωOplus=OperatorShift(file1)

# ∂ωOminuss = OperatorSandwich(ψ,∂ωOminus,weight,ψ)
# ∂ωOminus = ∂ωOminuss.Op;

# ∂ωOpluss = OperatorSandwich(ψ,∂ωOplus,weight,ψ)
# ∂ωOplus = ∂ωOpluss.Op;

# @show ∂ωOminus(1.0,0.0)

# println(Integrate(∂ωOminus, TheContour))

# qnmfunctionnew(s,l,m,n,a; qnm=qnm)


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