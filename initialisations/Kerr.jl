# Add the path to your local package to the LOAD_PATH
push!(LOAD_PATH, "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/ContourIntegrals.jl")
push!(LOAD_PATH, "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/SchwarzschildQNMShifts")

using ContourIntegrals
using KerrQNMShifts

### Use expressions for the operators, noted in their respective csvs
file1 = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/SchwarzschildQNMShifts/OperatorShifts/Schwarzschild/dwOpluscoefficients.csv"
file2 = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/SchwarzschildQNMShifts/OperatorShifts/Schwarzschild/dwOminuscoefficients.csv"
file3 = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/SchwarzschildQNMShifts/OperatorShifts/Schwarzschild/Hpluscoefficients.csv"
file4 = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/SchwarzschildQNMShifts/OperatorShifts/Schwarzschild/Hminuscoefficients.csv"
file5 = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/SchwarzschildQNMShifts/OperatorShifts/Schwarzschild/Ipluscoefficients.csv"
file6 = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/ComputeQNMs/SchwarzschildQNMShifts/OperatorShifts/Schwarzschild/Iminuscoefficients.csv"


# The below function prints out what each operator is from both file 1 and file2
∂ωOplus = OperatorShift(file1)
∂ωOminus= OperatorShift(file2)
Hplus = OperatorShift(file3)
Hminus = OperatorShift(file4)
Iplus = OperatorShift(file5)
Iminus = OperatorShift(file6)


### Get the integrand in the contour integral for the inner product computation
# qnmfunctionnew(s,l,m,n,a, plusminus);
ψplus = qnmfunctionnew(-2,2,2,0,0.5,"plus")
ψminus = qnmfunctionnew(-2,2,2,0,0.5,"minus")

weight = let s = ψplus.s , a= ψplus.a
    (r,z) -> sqrt(1-z^2)*(r^2+a^2-2*r)^s
end

∂ωOpluss = OperatorSandwich(ψplus,∂ωOplus,weight,ψplus)
∂ωOplus = ∂ωOpluss.Op;

∂ωOminuss = OperatorSandwich(ψminus,∂ωOminus,weight,ψminus)
∂ωOminus = ∂ωOminuss.Op;

Hpluss = OperatorSandwich(ψplus,Hplus,weight,ψplus)
Hplus = Hpluss.Op;

Hminuss = OperatorSandwich(ψminus,Hminus,weight,ψminus)
Hminus = Hminuss.Op;

Ipluss = OperatorSandwich(ψplus,Iplus,weight,ψminus)
Iplus = Ipluss.Op;

Iminuss = OperatorSandwich(ψminus,Iminus,weight,ψplus)
Iminus = Iminuss.Op;

println("Second Print")

### Compute the expressions so that they compile the first time
∂ωOplus(2.2,0.1) |> println
∂ωOminus(2.2,0.1) |> println
Hplus(2.2,0.1) |> println
Hminus(2.2,0.1) |> println
Iplus(2.2,0.1) |> println
Iminus(2.2,0.1) |> println

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
