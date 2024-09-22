push!(LOAD_PATH, "C:/Users/dwuuu/Documents/GitHub/ContourIntegrals.jl")
push!(LOAD_PATH, "C:/Users/dwuuu/Documents/GitHub/SchwarzschildQNMs")
push!(LOAD_PATH, "C:/Users/dwuuu/Documents/GitHub/KerrQuasinormalModes.jl")

using ContourIntegrals
using KerrQNMShifts
using KerrQuasinormalModes
using Test

@testset "TeukolskyOperator" begin
    
    ψ = qnmfunctionnew(-2,2,2,0,0.01)
    ψ0 = qnmfunctionnew(-2,2,2,0,0.)

    ω=ψ.ω
    ω0=ψ0.ω
    @show ω2=ω-ω0

    # Compile ψ
    ψ(1,.5)
    println("Past ψ compile")

    weight = let s = ψ0.s , a= ψ0.a
        (r,z) -> sqrt(1-z^2)*(r^2+a^2-2*r)^s
    end
    println("Past weight")

    ## Define the useful contours
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

    r₊0 = ψ0.R.r₊ ; r₋0 = ψ0.R.r₋ ; s0 = ψ0.s ; Δr0 = 0.1*(r₊0-r₋0); ϵ = eps(0.1);

    point10 = r₊0 + Δr0 - Δr0*im
    point20 = r₊0 - Δr0 - Δr0*im

    radial10 = SemiInfiniteLine(point10 , point10 + Δr0*im , false)
    angular0 = LineSegment(-1.0+100*ϵ , 1.0-100*ϵ , true) #to avoid the NaNs at the edges
    C10 = radial10 ⊗ angular0

    radial20 = LineSegment(point10,point20,true)
    C20 = radial20 ⊗ angular0

    radial30 = SemiInfiniteLine(point20 , point20 + Δr0*im , true)
    C30 = radial30 ⊗ angular0

    TheContour0 = C1⊕C2⊕C3

    println("Done Contours")

    ΣOplusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/WithSigmaOpluscoefficients.csv"
    Oplusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/Opluscoefficients.csv"
    Ominusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/Ominuscoefficients.csv"
    ΣdwOplusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/WithSigmadwOpluscoefficients.csv"
    dwOplusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/dwOpluscoefficients.csv"
    dwOminusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/dwOminuscoefficients.csv"

    ΣOplus = OperatorShift(ΣOplusfile)
    Σ∂ωOplus = OperatorShift(ΣdwOplusfile)
    Oplus = OperatorShift(Oplusfile)
    ∂ωOplus = OperatorShift(dwOplusfile)

    ΣOplusKerr = OperatorSandwich(ψ,ΣOplus,weight,ψ).Op
    ΣOplusSchw = OperatorSandwich(ψ0,ΣOplus,weight,ψ0).Op
    Σ∂ωOplusSchw = OperatorSandwich(ψ0,Σ∂ωOplus,weight,ψ0).Op
    OplusKerr = OperatorSandwich(ψ,Oplus,weight,ψ).Op
    OplusSchw = OperatorSandwich(ψ0,Oplus,weight,ψ0).Op
    ∂ωOplusSchw = OperatorSandwich(ψ0,∂ωOplus,weight,ψ0).Op

    println("Made Operators")

    @show ΣOplusKerr(3,0)
    @show ΣOplusSchw(3,0)
    @show Σ∂ωOplusSchw(3,0)

    @show OplusKerr(3,0)
    @show OplusSchw(3,0)
    @show ∂ωOplusSchw(3,0)

    println("Complied Operators")
   
    @show Σ∂ω𝒪plusSchw= Integrate(Σ∂ωOplusSchw, TheContour0)[1]
    println("Done first integral")
    @show Σ𝒪plusKerr= Integrate(ΣOplusKerr, TheContour)[1]
    @show Σ𝒪plusKerr
    @show ω2*Σ∂ω𝒪plusSchw
    println("Done WithSigma")

    @show ∂ω𝒪plusSchw= Integrate(∂ωOplusSchw, TheContour0)[1]
    println("Done first integral")
    @show 𝒪plusKerr= Integrate(OplusKerr, TheContour)[1]
    @show 𝒪plusKerr
    @show ω2*∂ω𝒪plusSchw
    

end

# # Random number generator
# function random_complex(length)
#     complex_numbers = []
#     for _ in 0:length
#         lower_bound =  r₋ - 10*im
#         upper_bound = 10 + 10*im
#         real_part = rand() * (real(upper_bound) - real(lower_bound)) + real(lower_bound)
#         imag_part = rand() * (imag(upper_bound) - imag(lower_bound)) + imag(lower_bound)
#         random = real_part + imag_part * im
#         push!(complex_numbers, random)
#     end
#     return complex_numbers
# end

# random_r = random_complex(10)

# random_z = rand(10)

# @testset "WavefunctionProperties" begin
#     ₋₂ψ₂₂₀₊ = qnmfunctionnew(-2,2,2,0,0.5) 
#     ₋₂ψ₂₂₀₋ = qnmfunctionnew(-2,2,2,0,0.5,modesign="minus")
#     ₊₂ψ₂₂₀₊ = qnmfunctionnew(2,2,2,0,0.5)
#     ₊₂ψ₂₂₀₋= qnmfunctionnew(2,2,2,0,0.5,modesign="minus")
#     r₊=₋₂ψ₂₂₀₊.R.r₊
#     r₋=₋₂ψ₂₂₀₊.R.r₋

#     println("Check Radial minus mode equality")
#     for i in 1:10
#         @assert abs(₋₂ψ₂₂₀₊(random_r[i])-conj(₋₂ψ₂₂₀₋(random_r[i])))<=10^(-8)
#     end

#     println("Check Radial exponential decay")
#     @assert abs(₋₂ψ₂₂₀₊(r₊+.1+100*im))  <= 10^(-8)
#     @assert abs(₋₂ψ₂₂₀₊(r₊-.1+100*im)) <= 10^(-8)
#     @assert abs(₋₂ψ₂₂₀₋(r₊+.1+100*im)) <= 10^(-8)
#     @assert abs(₋₂ψ₂₂₀₋(r₊-.1+100*im)) <= 10^(-8)

#     println("Check Angular minus mode equality")
#     for i in 1:10
#         @assert abs(₋₂ψ₂₂₀₊.S(random_z[i])-conj(₊₂ψ₂₂₀₋.S(random_z[i]))) <= 10^(-8)
#     end

#     println("Check derivatives")
#     δ=10^(-6)
#     for i in 1:10
#         @assert abs(∂r(₋₂ψ₂₂₀₊)(random_r[i])-(₋₂ψ₂₂₀₊(random_r[i]+δ)-₋₂ψ₂₂₀₊(random_r[i]))/δ^2))<=10^(-8)
#         @assert abs(∂θ(₋₂ψ₂₂₀₊)(random_r[i],random_z[i])-(₋₂ψ₂₂₀₊(random_r[i],random_z[i]+δ)-₋₂ψ₂₂₀₊(random_r[i],random_z[i]))/δ^2))<=10^(-8)

#         # Second-order derivative test
#         finite_diff_2nd_r = (₋₂ψ₂₂₀₊(random_r[i]+δ) - 2*₋₂ψ₂₂₀₊(random_r[i]) + ₋₂ψ₂₂₀₊(random_r[i]-δ)) / δ^2
#         @assert abs(∂r(∂r(₋₂ψ₂₂₀₊))(random_r[i]) - finite_diff_2nd_r) <= 10^(-8)

#         finite_diff_2nd_θ = (₋₂ψ₂₂₀₊(random_r[i],random_z[i]+δ) - 2*₋₂ψ₂₂₀₊(random_r[i],random_z[i]) + ₋₂ψ₂₂₀₊(random_r[i],random_z[i]-δ)) / δ^2
#         @assert abs(∂r(∂r(₋₂ψ₂₂₀₊))(random_r[i],random_z[i]) - finite_diff_2nd_θ) <= 10^(-8)

#         # Mixed derivative test
#         finite_diff_rθ = (₋₂ψ₂₂₀₊(random_r[i]+δ, random_z[i]+δ) - ₋₂ψ₂₂₀₊(random_r[i]+δ, random_z[i]-δ) - ₋₂ψ₂₂₀₊(random_r[i]-δ, random_z[i]+δ) + ₋₂ψ₂₂₀₊(random_r[i]-δ, random_z[i]-δ)) / (4*δ^2)
#          @assert abs(∂r(∂θ(₋₂ψ₂₂₀₊))(random_r[i],random_z[i]) - finite_diff_rθ) <= 10^(-8)

#     end
  
# end
