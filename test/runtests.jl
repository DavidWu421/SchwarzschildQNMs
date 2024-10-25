push!(LOAD_PATH, "C:/Users/dwuuu/Documents/GitHub/ContourIntegrals.jl")
push!(LOAD_PATH, "C:/Users/dwuuu/Documents/GitHub/SchwarzschildQNMs")
push!(LOAD_PATH, "C:/Users/dwuuu/Documents/GitHub/KerrQuasinormalModes.jl")

using ContourIntegrals
using KerrQNMShifts
using KerrQuasinormalModes
using Test

function compute_derivative_matrix(ψ, r, θ)
    # Initialize a 5x5 matrix to hold the derivatives
    derivatives = Matrix{ComplexF64}(undef, 5, 5)
    for i in 0:4
        for j in 0:4
            temp_ψ = ψ
            for _ in 1:i
                temp_ψ = ∂r(temp_ψ)
            end
            for _ in 1:j
                temp_ψ = ∂θ(temp_ψ)
            end
            derivatives[i+1, j+1] = temp_ψ(r, θ)
        end
    end
    
    return derivatives
end

function matrix_to_mathematica_format(matrix)
    # Convert the matrix to a string in Mathematica's List format
    mat_str = string("{{", join([join([string(el) for el in row], ", ") for row in eachrow(matrix)], "}, {"), "}}")
    return mat_str
end

# @testset "TeukolskyOperator" begin
    
    ψ = qnmfunctionnew(-2,2,2,0,0.01)
    ψ0 = qnmfunctionnew(-2,2,2,0,0.)

    ψm = qnmfunctionnew(-2,2,2,0,0.01,modesign="minus")
    ψ0m = qnmfunctionnew(-2,2,2,0,0.,modesign="minus")

    ψtest = qnmfunctionnew(2,5,3,0,0.)
    ψtest.ω=ψ0.ω
    ψtest.a=ψ0.a
    ψtest.m=ψ0.m
    ψtest.s=ψ0.s


    ω=ψ.ω
    ω0=ψ0.ω
    @show ω2=ω-ω0

    # Compile ψ
    ψ(1,.5)
    println("Past ψ compile")

    Σ = let a= ψ.a
        (r,z) -> Complex(r)^2+a^2*z^2
    end
    Δ = let a= ψ.a
        (r,z) -> Complex(r)^2+a^2-2*r
    end
    ζ = let a= ψ.a
        (r,z) -> r-im*a*z
    end

    weight = let s = ψ.s , a= ψ.a
        (r,z) ->2^(5+s)*ζ(r,z)^(8+2*s)*Σ(r,z)*sqrt(1-z^2)/((Δ(r,z))^2)
    end


    Σ0 = let a= ψ0.a
        (r,z) -> Complex(r)^2+a^2*z^2
    end
    Δ0 = let a= ψ0.a
        (r,z) -> Complex(r)^2+a^2-2*r
    end
    ζ0 = let a= ψ0.a
        (r,z) -> r-im*a*z
    end
    weight0 = let s = ψ0.s , a= ψ0.a
        (r,z) ->2^(5+s)*ζ0(r,z)^(8+2*s)*Σ0(r,z)*sqrt(1-z^2)/((Δ0(r,z))^2)
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

    KerrOplusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/KerrOcoefficients.csv"
    Oplusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/Ocoefficients.csv"
    dwOplusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/dwOcoefficients.csv"

    KerrOplus = OperatorShift(KerrOplusfile)
    KerrOminus = OperatorShift(KerrOplusfile)
    Oplus = OperatorShift(Oplusfile)
    ∂ωOplus = OperatorShift(dwOplusfile)
    Ominus = OperatorShift(Oplusfile)
    ∂ωOminus = OperatorShift(dwOplusfile)

    println("Made operator shifts")

    # KerrOplusKerr = OperatorSandwich(ψ,KerrOplus,weight,ψ).Op
    # KerrOminusKerr = OperatorSandwich(ψm,KerrOminus,weight,ψm).Op
    HermiticityTest1 = OperatorSandwich(ψ0,Oplus,weight0,ψtest).Op
    HermiticityTest2 = OperatorSandwich(ψtest,Oplus,weight0,ψ0).Op

    # println("Made the first set of operators")
    
    # OplusKerr = OperatorSandwich(ψ0,Oplus,weight0,ψ).Op
    # OplusSchw = OperatorSandwich(ψ0,Oplus,weight0,ψ0).Op
    # ∂ωOplusSchw = OperatorSandwich(ψ0,∂ωOplus,weight0,ψ0).Op
    # OminusKerr = OperatorSandwich(ψ0m,Ominus,weight0,ψm).Op
    # OminusSchw = OperatorSandwich(ψ0m,Ominus,weight0,ψ0m).Op
    # ∂ωOminusSchw = OperatorSandwich(ψ0m,∂ωOminus,weight0,ψ0m).Op

    println("Made Operators")

    @show HermiticityTest1(3,0.5)
    @show HermiticityTest2(3,0.5)

    # @show KerrOplusKerr(3,0.5)
    # @show KerrOminusKerr(3,0.5)

    # @show OplusKerr(3,0.5)
    # @show OplusSchw(3,0.5)
    # @show ∂ωOplusSchw(3,0.5)

    # @show OminusKerr(3,0.5)
    # @show OminusSchw(3,0.5)
    # @show ∂ωOminusSchw(3,0.5)

    println("Complied Operators")

    HermiticityTestResult1 = Integrate(HermiticityTest1,TheContour0)[1]
    @show HermiticityTestResult1
    HermiticityTestResult2 = Integrate(HermiticityTest2,TheContour0)[1]
    @show HermiticityTestResult2
   
    # ∂ω𝒪plusSchw= Integrate(∂ωOplusSchw, TheContour0)[1]
    # 𝒪plusKerr= Integrate(OplusKerr, TheContour)[1]
    # @show 𝒪plusKerr
    # @show ∂ω𝒪plusSchw
    # @show ω2*∂ω𝒪plusSchw

    # ∂ω𝒪minusSchw= Integrate(∂ωOminusSchw, TheContour0)[1]
    # 𝒪minusKerr= Integrate(OminusKerr, TheContour)[1]
    # @show 𝒪minusKerr
    # @show ∂ω𝒪minusSchw
    # @show -conj(ω2)*∂ω𝒪minusSchw
    
# end

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

#         @assert isapprox(∂r(₋₂ψ₂₂₀₊)(random_r[i]), (₋₂ψ₂₂₀₊(random_r[i]+δ)- ₋₂ψ₂₂₀₊(random_r[i]))/δ, rtol=1e-5)
#         @assert isapprox(∂θ(₋₂ψ₂₂₀₊)(random_r[i],random_z[i]), -sqrt(1-random_z[i]^2)* (₋₂ψ₂₂₀₊(random_r[i],random_z[i]+δ)- ₋₂ψ₂₂₀₊(random_r[i],random_z[i]))/δ, rtol=1e-5)
        
#         # Second-order derivative test
#         finite_diff_2nd_r = (₋₂ψ₂₂₀₊(random_r[i]+δ) - 2*(₋₂ψ₂₂₀₊(random_r[i])) + ₋₂ψ₂₂₀₊(random_r[i]-δ)) / δ^2
#         @assert isapprox(∂r(∂r(₋₂ψ₂₂₀₊))(random_r[i]), finite_diff_2nd_r, rtol=1e-1)

#         finite_diff_2nd_θ = (1-random_z[i]^2)*(₋₂ψ₂₂₀₊(random_r[i],random_z[i]+δ) - 2*(₋₂ψ₂₂₀₊(random_r[i],random_z[i])) + ₋₂ψ₂₂₀₊(random_r[i],random_z[i]-δ)) / δ^2 - random_z[i]*(₋₂ψ₂₂₀₊(random_r[i],random_z[i]+δ)- ₋₂ψ₂₂₀₊(random_r[i],random_z[i]))/δ
#         @assert isapprox(∂θ(∂θ(₋₂ψ₂₂₀₊))(random_r[i],random_z[i]), finite_diff_2nd_θ, rtol=1e-2)

#         # Mixed derivative test
#         finite_diff_rθ = -sqrt(1-random_z[i]^2)*(₋₂ψ₂₂₀₊(random_r[i]+δ, random_z[i]+δ) - ₋₂ψ₂₂₀₊(random_r[i]+δ, random_z[i]-δ) - ₋₂ψ₂₂₀₊(random_r[i]-δ, random_z[i]+δ) + ₋₂ψ₂₂₀₊(random_r[i]-δ, random_z[i]-δ)) / (4*δ^2)
#         @assert isapprox(∂r(∂θ(₋₂ψ₂₂₀₊))(random_r[i],random_z[i]), finite_diff_rθ, rtol=1e-3)
#         @assert isapprox(∂θ(∂r(₋₂ψ₂₂₀₊))(random_r[i],random_z[i]), finite_diff_rθ, rtol=1e-3)


#         #Minus mode derivatives

#         @assert isapprox(∂r(₋₂ψ₂₂₀₋)(random_r[i]), (₋₂ψ₂₂₀₋(random_r[i]+δ)- ₋₂ψ₂₂₀₋(random_r[i]))/δ, rtol=1e-5)
#         @assert isapprox(∂θ(₋₂ψ₂₂₀₋)(random_r[i],random_z[i]), -sqrt(1-random_z[i]^2)* (₋₂ψ₂₂₀₋(random_r[i],random_z[i]+δ)- ₋₂ψ₂₂₀₋(random_r[i],random_z[i]))/δ, rtol=1e-4)
        
#         # Second-order derivative test
#         finite_diff_2nd_r = (₋₂ψ₂₂₀₋(random_r[i]+δ) - 2*(₋₂ψ₂₂₀₋(random_r[i])) + ₋₂ψ₂₂₀₋(random_r[i]-δ)) / δ^2
#         @assert isapprox(∂r(∂r(₋₂ψ₂₂₀₋))(random_r[i]), finite_diff_2nd_r, rtol=1e-1)

#         finite_diff_2nd_θ = (1-random_z[i]^2)*(₋₂ψ₂₂₀₋(random_r[i],random_z[i]+δ) - 2*(₋₂ψ₂₂₀₋(random_r[i],random_z[i])) + ₋₂ψ₂₂₀₋(random_r[i],random_z[i]-δ)) / δ^2 - random_z[i]*(₋₂ψ₂₂₀₋(random_r[i],random_z[i]+δ)- ₋₂ψ₂₂₀₋(random_r[i],random_z[i]))/δ
#         @assert isapprox(∂θ(∂θ(₋₂ψ₂₂₀₋))(random_r[i],random_z[i]), finite_diff_2nd_θ, rtol=1e-2)

#         # Mixed derivative test
#         finite_diff_rθ = -sqrt(1-random_z[i]^2)*(₋₂ψ₂₂₀₋(random_r[i]+δ, random_z[i]+δ) - ₋₂ψ₂₂₀₋(random_r[i]+δ, random_z[i]-δ) - ₋₂ψ₂₂₀₋(random_r[i]-δ, random_z[i]+δ) + ₋₂ψ₂₂₀₋(random_r[i]-δ, random_z[i]-δ)) / (4*δ^2)
#         @assert isapprox(∂r(∂θ(₋₂ψ₂₂₀₋))(random_r[i],random_z[i]), finite_diff_rθ, rtol=1e-3)
#         @assert isapprox(∂θ(∂r(₋₂ψ₂₂₀₋))(random_r[i],random_z[i]), finite_diff_rθ, rtol=1e-3)

         

#     end
  
# end
