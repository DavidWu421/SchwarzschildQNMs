# Add the path to your local package to the LOAD_PATH
push!(LOAD_PATH, "C:/Users/dwuuu/Documents/GitHub/ContourIntegrals.jl")
push!(LOAD_PATH, "C:/Users/dwuuu/Documents/GitHub/SchwarzschildQNMs")
push!(LOAD_PATH, "C:/Users/dwuuu/Documents/GitHub/KerrQuasinormalModes.jl")

using ContourIntegrals
using KerrQNMShifts
using KerrQuasinormalModes


### Use expressions for the operators, noted in their respective csvs
# filetest = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/SchwarzschildQNMShifts/OperatorShifts/Schwarzschild/Testcoefficients.csv"
file1 = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/SchwarzschildQNMShifts/OperatorShifts/Schwarzschild/dwOpluscoefficients.csv"
file2 = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/SchwarzschildQNMShifts/OperatorShifts/Schwarzschild/dwOminuscoefficients.csv"
file3 = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/SchwarzschildQNMShifts/OperatorShifts/Schwarzschild/Hpluscoefficients.csv"
file4 = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/SchwarzschildQNMShifts/OperatorShifts/Schwarzschild/Hminuscoefficients.csv"
file5 = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/SchwarzschildQNMShifts/OperatorShifts/Schwarzschild/Ipluscoefficients.csv"
file6 = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/SchwarzschildQNMShifts/OperatorShifts/Schwarzschild/Iminuscoefficients.csv"

# Test = OperatorShift(filetest)
∂ωOplus = OperatorShift(file1)
∂ωOminus= OperatorShift(file2)
Hplus = OperatorShift(file3)
Hminus = OperatorShift(file4)
Iplus = OperatorShift(file5)
Iminus = OperatorShift(file6)

ψ1 = qnmfunctionnew(-2,2,2,0,0.) 
ψ2= qnmfunctionnew(-2,2,2,0,0.,modesign="minus")
ψ3 = qnmfunctionnew(2,2,2,0,0.)
ψ4= qnmfunctionnew(2,2,2,0,0.,modesign="minus")

r₊=ψ1.R.r₊
r₋=ψ1.R.r₋

println("Check Radial minus mode equality")
@show ψ1.R(r₊+.1)-conj(ψ2.R(r₊+.1))
@show ψ1.R(r₊+.1+im)-conj(ψ2.R(r₊+.1+im))
@show ψ1.R(r₊-.1)-conj(ψ2.R(r₊-.1))
@show ψ1.R(r₊-.1+im)-conj(ψ2.R(r₊-.1+im))
println("Check Radial exponential decay")
@show ψ1.R(r₊+.1+100*im)
@show ψ2.R(r₊+.1+100*im)
@show ψ1.R(r₊-.1+100*im)
@show ψ2.R(r₊-.1+100*im)

println("Check Angular minus mode equality")
@show ψ3.S(.5)-conj(ψ2.S(.5))

α1=ψ1.R.α
η1=ψ1.R.η
ξ1=ψ1.R.ξ
ζ1=ψ1.R.ζ
α2=ψ2.R.α
η2=ψ2.R.η
ξ2=ψ2.R.ξ
ζ2=ψ2.R.ζ


# Compile ψ
ψ1(1,.5)
println("Past ψ compile")

weight = let s = ψ1.s , a= ψ1.a
    (r,z) -> sqrt(1-z^2)*(r^2+a^2-2*r)^s
end

println("Past weight")

∂ωOpluss1 = OperatorSandwich(ψ1,∂ωOplus,weight,ψ1)
∂ωOplus1 = ∂ωOpluss1.Op;

∂ωOminuss1 = OperatorSandwich(ψ2,∂ωOminus,weight,ψ2)
∂ωOminus1 = ∂ωOminuss1.Op;

println("Made O operators")

Hpluss1 = OperatorSandwich(ψ1,Hplus,weight,ψ1)
Hplus1 = Hpluss1.Op;

Hminuss1 = OperatorSandwich(ψ2,Hminus,weight,ψ2)
Hminus1 = Hminuss1.Op;

println("Made H operators")

Ipluss1 = OperatorSandwich(ψ1,Iplus,weight,ψ2)
Iplus1 = Ipluss1.Op;

Iminuss1 = OperatorSandwich(ψ2,Iminus,weight,ψ1)
Iminus1 = Iminuss1.Op;

println("Made I operators")

# Compile the operator sandwiches

@show ∂ωOplus1(1,.5)
@show ∂ωOminus1(1,.5)
@show Hplus1(1,.5)
@show Hminus1(1,.5)
@show Iplus1(1,.5,isconjugate=true)
@show Iminus1(1,.5,isconjugate=true)

## Define the useful contours
r₊ = ψ1.R.r₊ ; r₋ = ψ1.R.r₋ ; s = ψ1.s ; Δr = 0.1*(r₊-r₋); ϵ = eps(0.1);

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

𝒪plus= Integrate(∂ωOplus1, TheContour)[1]
𝒪minus= conj(Integrate(∂ωOminus1, TheContour)[1])
println(𝒪minus)
println("Done O's")
ℋplus= Integrate(Hplus1, TheContour)[1]
println("Done Hplus")
ℋminus= conj(Integrate(Hminus1, TheContour)[1])
println("Done Hminus")
ℐplus= Integrate(Iplus1, TheContour,isconjugate=true)[1]
println("Done Iplus")
ℐminus= conj(Integrate(Iminus1, TheContour,isconjugate=true)[1])
println("Done Iminus")

ω2shift=Computeω2(𝒪plus,𝒪minus,ℋplus,ℋminus,ℐplus,ℐminus,ψ1)
@show ω2shift