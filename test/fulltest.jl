using ContourIntegrals
using KerrQNMShifts
using KerrQuasinormalModes
using Test
using Zygote

println("Done usings")

@testset "HermiticityTest" begin
    println("Started HermiticityTest: ")
    pert_a=.1
    
    ψ = qnmfunctionnew(-2,2,2,0,0.)

    ψm = qnmfunctionnew(-2,2,2,0,0.,modesign="minus")

    ψ1 = qnmfunctionnew(-2,4,2,0,0.1)
    ψm1 = qnmfunctionnew(-2,4,2,0,0.1,modesign="minus")
    ψ2 = qnmfunctionnew(-2,4,2,0,0.)
    ψm2 = qnmfunctionnew(-2,4,2,0,0.,modesign="minus")
    ψ3 = qnmfunctionnew(-2,5,2,0,0.)
    ψm3 = qnmfunctionnew(-2,5,2,0,0.,modesign="minus")
    ψ4 = qnmfunctionnew(-2,5,3,0,0.)
    ψm4 = qnmfunctionnew(-2,5,3,0,0.,modesign="minus")

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

    weightplus = let a= ψ.a
        (r,z) ->8*ζ(r,z)^(4)*Σ(r,z)/((Δ(r,z))^2)
    end

    weightminus = let a= ψ.a
        (r,z) ->8*ζ(conj(r),z)^(4)*Σ(conj(r),z)/((Δ(conj(r),z))^2)
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
    println("Done Contours")

    Oplusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/HermTestOpluscoefficients.csv"
    Ominusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/HermTestOminuscoefficients.csv"

    Oplus=OperatorShift(Oplusfile)
    Ominus=OperatorShift(Ominusfile)

    println("Made operator shifts")
    
    
    OplusSchw1 = OperatorSandwich(ψ,Oplus,weightplus,ψ1).Op
    OminusSchw1 = OperatorSandwich(ψm,Ominus,weightminus,ψ1).Op
    OplusSchw2 = OperatorSandwich(ψ,Oplus,weightplus,ψ2).Op
    OminusSchw2 = OperatorSandwich(ψm,Ominus,weightminus,ψ2).Op
    OplusSchw3 = OperatorSandwich(ψ,Oplus,weightplus,ψ3).Op
    OminusSchw3 = OperatorSandwich(ψm,Ominus,weightminus,ψ3).Op
    OplusSchw4 = OperatorSandwich(ψ,Oplus,weightplus,ψ4).Op
    OminusSchw4 = OperatorSandwich(ψm,Ominus,weightminus,ψ4).Op

    # OminusSchw0check = OperatorSandwich(ψmweird,Ominus,weightminus,ψm).Op

    println("Made Operators")

    # @show OminusSchw0check(8+im,.6,pertparam=pert_a)

    @show OplusSchw1(8+im,0.6,pertparam=pert_a)
    @show OminusSchw1(8+im,0.6,pertparam=pert_a)
    @show OplusSchw2(8+im,0.6,pertparam=pert_a)
    @show OminusSchw2(8+im,0.6,pertparam=pert_a)
    @show OplusSchw3(8+im,0.6,pertparam=pert_a)
    @show OminusSchw3(8+im,0.6,pertparam=pert_a)
    @show OplusSchw4(8+im,0.6,pertparam=pert_a)
    @show OminusSchw4(8+im,0.6,pertparam=pert_a)

    println("Complied Operators")

    𝒪plusSchw1= Integrate(OplusSchw1, TheContour,pertparam=pert_a,abstol=1e-6)[1]
    𝒪minusSchw1= Integrate(OminusSchw1, TheContour,pertparam=pert_a,abstol=1e-6)[1]
    𝒪plusSchw2= Integrate(OplusSchw2, TheContour,pertparam=pert_a,abstol=1e-6)[1]
    𝒪minusSchw2= Integrate(OminusSchw2, TheContour,pertparam=pert_a,abstol=1e-6)[1]
    𝒪plusSchw3= Integrate(OplusSchw3, TheContour,pertparam=pert_a,abstol=1e-6)[1]
    𝒪minusSchw3= Integrate(OminusSchw3, TheContour,pertparam=pert_a,abstol=1e-6)[1]
    𝒪plusSchw4= Integrate(OplusSchw4, TheContour,pertparam=pert_a,abstol=1e-6)[1]
    𝒪minusSchw4= Integrate(OminusSchw4, TheContour,pertparam=pert_a,abstol=1e-6)[1]

    @show 𝒪plusSchw1
    @show 𝒪minusSchw1
    @show 𝒪plusSchw2
    @show 𝒪minusSchw2
    @show 𝒪plusSchw3
    @show 𝒪minusSchw3
    @show 𝒪plusSchw4
    @show 𝒪minusSchw4
end

@testset "KerrAsPerturb" begin

    println("Started KerrAsPertub Test: ")
    pert_a=.01
    
    ψ = qnmfunctionnew(-2,2,2,0,0.)
    ψconj = qnmfunctionnew(-2,2,2,0,0.,is_conjugate=true)

    ψm = qnmfunctionnew(-2,2,2,0,0.,modesign="minus")
    ψmconj = qnmfunctionnew(-2,2,2,0,0.,modesign="minus",is_conjugate=true)

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

    weightplus = let a= ψ.a
        (r,z) ->8*ζ(r,z)^(4)*Σ(r,z)/((Δ(r,z))^2)
    end

    weightminus = let a= ψ.a
        (r,z) ->8*ζ(conj(r),z)^(4)*Σ(conj(r),z)/((Δ(conj(r),z))^2)
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
    println("Done Contours")


    dwOplusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/dwOpluscoefficients.csv"
    dwOminusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/dwOminuscoefficients.csv"
    Hplusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/Hpluscoefficients.csv"
    Hminusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/Hminuscoefficients.csv"
    Iplusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/Ipluscoefficients.csv"
    Iminusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/Iminuscoefficients.csv"

    ∂ωOplus = OperatorShift(dwOplusfile)
    ∂ωOminus = OperatorShift(dwOminusfile)
    Hplus = OperatorShift(Hplusfile)
    Hminus = OperatorShift(Hminusfile)
    Iplus = OperatorShift(Iplusfile)
    Iminus = OperatorShift(Iminusfile)

    println("Made operator shifts")
    
    ∂ωOplusSchw = OperatorSandwich(ψ,∂ωOplus,weightplus,ψ).Op
    ∂ωOminusSchw = OperatorSandwich(ψm,∂ωOminus,weightminus,ψm).Op
    HplusSchw = OperatorSandwich(ψ,Hplus,weightplus,ψ).Op
    HminusSchw = OperatorSandwich(ψm,Hminus,weightminus,ψm).Op
    IplusSchw = OperatorSandwich(ψ,Iplus,weightplus,ψmconj).Op
    IminusSchw = OperatorSandwich(ψm,Iminus,weightminus,ψconj).Op

    println("Made Operators")

    @show ∂ωOplusSchw(8+im,0.6,pertparam=pert_a)
    @show ∂ωOminusSchw(8+im,0.6,pertparam=pert_a)
    @show HplusSchw(8+im,0.6,pertparam=pert_a)
    @show HminusSchw(8+im,0.6,pertparam=pert_a)
    @show IplusSchw(8+im,0.6,pertparam=pert_a)
    @show IminusSchw(8+im,0.6,pertparam=pert_a)

    println("Complied Operators")

    ∂ω𝒪plusSchw= Integrate(∂ωOplusSchw, TheContour,pertparam=pert_a,abstol=1e-6)[1]
    ∂ω𝒪minusSchw= conj(Integrate(∂ωOminusSchw, TheContour,pertparam=pert_a,abstol=1e-6)[1])
    ℋplusSchw= Integrate(HplusSchw, TheContour,pertparam=pert_a,abstol=1e-6)[1]
    ℋminusSchw= conj(Integrate(HminusSchw, TheContour,pertparam=pert_a,abstol=1e-6)[1])
    ℐplusSchw= Integrate(IplusSchw, TheContour,pertparam=pert_a,abstol=1e-6)[1]
    ℐminusSchw= conj(Integrate(IminusSchw, TheContour,pertparam=pert_a,abstol=1e-6)[1])

    @show ∂ω𝒪plusSchw
    @show ∂ω𝒪minusSchw
    @show ℋplusSchw
    @show ℋminusSchw
    @show ℐplusSchw
    @show ℐminusSchw 
    
    ω2s=Computeω2(∂ω𝒪plusSchw,∂ω𝒪minusSchw,ℋplusSchw,ℋminusSchw,ℐplusSchw,ℐminusSchw,ψ)
    @show ω2s

end

@testset "QNM functions" begin
    # Random number generator
    function random_complex(length)
        complex_numbers = []
        for _ in 0:length
            lower_bound =  r₋ - 10*im
            upper_bound = 10 + 10*im
            real_part = rand() * (real(upper_bound) - real(lower_bound)) + real(lower_bound)
            imag_part = rand() * (imag(upper_bound) - imag(lower_bound)) + imag(lower_bound)
            random = real_part + imag_part * im
            push!(complex_numbers, random)
        end
        return complex_numbers
    end

    ₋₂ψ₂₂₀₊ = qnmfunctionnew(-2,2,2,0,0.5) 
    ₋₂ψ₂₂₀₋ = qnmfunctionnew(-2,2,2,0,0.5,modesign="minus")
    ₊₂ψ₂₂₀₊ = qnmfunctionnew(2,2,2,0,0.5)
    ₊₂ψ₂₂₀₋= qnmfunctionnew(2,2,2,0,0.5,modesign="minus")
    ₋₂ψ₂₂₀₊conj = qnmfunctionnew(-2,2,2,0,0.5,is_conjugate=true) 
    ₋₂ψ₂₂₀₋conj = qnmfunctionnew(-2,2,2,0,0.5,modesign="minus",is_conjugate=true)
    ₊₂ψ₂₂₀₊conj = qnmfunctionnew(2,2,2,0,0.5,is_conjugate=true)
    ₊₂ψ₂₂₀₋conj= qnmfunctionnew(2,2,2,0,0.5,modesign="minus",is_conjugate=true)
    r₊=₋₂ψ₂₂₀₊.R.r₊
    r₋=₋₂ψ₂₂₀₊.R.r₋

    random_r = random_complex(10)
    random_z = rand(10)

    println("Check Radial minus mode equality")
    for i in 1:10
        @assert abs(₋₂ψ₂₂₀₊(random_r[i])-conj(₋₂ψ₂₂₀₋(random_r[i])))<=10^(-8)
        @assert abs(₋₂ψ₂₂₀₊conj(random_r[i])-conj(₋₂ψ₂₂₀₋conj(random_r[i])))<=10^(-8)
        @assert abs(₊₂ψ₂₂₀₊(random_r[i])-conj(₊₂ψ₂₂₀₋(random_r[i])))<=10^(-8)
        @assert abs(₊₂ψ₂₂₀₊conj(random_r[i])-conj(₊₂ψ₂₂₀₋conj(random_r[i])))<=10^(-8)
    end

    println("Check Radial exponential decay")
    @assert abs(₋₂ψ₂₂₀₊(r₊+.1+100*im)) <= 10^(-8)
    @assert abs(₊₂ψ₂₂₀₊(r₊+.1+100*im)) <= 10^(-8)
    @assert abs(₋₂ψ₂₂₀₋(r₊+.1+100*im)) <= 10^(-8)
    @assert abs(₊₂ψ₂₂₀₋(r₊+.1+100*im)) <= 10^(-8)
    @assert abs(₋₂ψ₂₂₀₊conj(r₊+.1+100*im)) <= 10^(-8)
    @assert abs(₊₂ψ₂₂₀₊conj(r₊+.1+100*im)) <= 10^(-8)
    @assert abs(₋₂ψ₂₂₀₋conj(r₊+.1+100*im)) <= 10^(-8)
    @assert abs(₊₂ψ₂₂₀₋conj(r₊+.1+100*im)) <= 10^(-8)

    println("Check Angular minus mode equality")
    for i in 1:10
        @assert abs(₋₂ψ₂₂₀₊.S(random_z[i])-conj(₊₂ψ₂₂₀₋.S(random_z[i]))) <= 10^(-8)
        @assert abs(₋₂ψ₂₂₀₊conj.S(random_z[i])-conj(₊₂ψ₂₂₀₋conj.S(random_z[i]))) <= 10^(-8)
    end

    println("Check derivatives")
    δ=10^(-6)
    for i in 1:10
        @assert isapprox(∂r(₋₂ψ₂₂₀₊)(random_r[i]), (₋₂ψ₂₂₀₊(random_r[i]+δ)- ₋₂ψ₂₂₀₊(random_r[i]))/δ, rtol=1e-5)
        @assert isapprox(∂r(₋₂ψ₂₂₀₊conj)(random_r[i]), (₋₂ψ₂₂₀₊conj(random_r[i]+δ)- ₋₂ψ₂₂₀₊conj(random_r[i]))/δ, rtol=1e-5)
        @assert isapprox(∂r(₋₂ψ₂₂₀₋)(random_r[i]), (₋₂ψ₂₂₀₋(random_r[i]+δ)- ₋₂ψ₂₂₀₋(random_r[i]))/δ, rtol=1e-5)
        @assert isapprox(∂r(₋₂ψ₂₂₀₋conj)(random_r[i]), (₋₂ψ₂₂₀₋conj(random_r[i]+δ)- ₋₂ψ₂₂₀₋conj(random_r[i]))/δ, rtol=1e-5)
        @assert isapprox(∂r(₊₂ψ₂₂₀₊)(random_r[i]), (₊₂ψ₂₂₀₊(random_r[i]+δ)- ₊₂ψ₂₂₀₊(random_r[i]))/δ, rtol=1e-5)
        @assert isapprox(∂r(₊₂ψ₂₂₀₊conj)(random_r[i]), (₊₂ψ₂₂₀₊conj(random_r[i]+δ)- ₊₂ψ₂₂₀₊conj(random_r[i]))/δ, rtol=1e-5)
        @assert isapprox(∂r(₊₂ψ₂₂₀₋)(random_r[i]), (₊₂ψ₂₂₀₋(random_r[i]+δ)- ₊₂ψ₂₂₀₋(random_r[i]))/δ, rtol=1e-5)
        @assert isapprox(∂r(₊₂ψ₂₂₀₋conj)(random_r[i]), (₊₂ψ₂₂₀₋conj(random_r[i]+δ)- ₊₂ψ₂₂₀₋conj(random_r[i]))/δ, rtol=1e-5)

        @assert isapprox(∂θ(₋₂ψ₂₂₀₊)(random_r[i],random_z[i]), -sqrt(1-random_z[i]^2)* (₋₂ψ₂₂₀₊(random_r[i],random_z[i]+δ)- ₋₂ψ₂₂₀₊(random_r[i],random_z[i]))/δ, rtol=1e-5)
        @assert isapprox(∂θ(₋₂ψ₂₂₀₊conj)(random_r[i],random_z[i]), -sqrt(1-random_z[i]^2)* (₋₂ψ₂₂₀₊conj(random_r[i],random_z[i]+δ)- ₋₂ψ₂₂₀₊conj(random_r[i],random_z[i]))/δ, rtol=1e-5)
        @assert isapprox(∂θ(₋₂ψ₂₂₀₋)(random_r[i],random_z[i]), -sqrt(1-random_z[i]^2)* (₋₂ψ₂₂₀₋(random_r[i],random_z[i]+δ)- ₋₂ψ₂₂₀₋(random_r[i],random_z[i]))/δ, rtol=1e-5)
        @assert isapprox(∂θ(₋₂ψ₂₂₀₋conj)(random_r[i],random_z[i]), -sqrt(1-random_z[i]^2)* (₋₂ψ₂₂₀₋conj(random_r[i],random_z[i]+δ)- ₋₂ψ₂₂₀₋conj(random_r[i],random_z[i]))/δ, rtol=1e-5)
        @assert isapprox(∂θ(₊₂ψ₂₂₀₊)(random_r[i],random_z[i]), -sqrt(1-random_z[i]^2)* (₊₂ψ₂₂₀₊(random_r[i],random_z[i]+δ)- ₊₂ψ₂₂₀₊(random_r[i],random_z[i]))/δ, rtol=1e-5)
        @assert isapprox(∂θ(₊₂ψ₂₂₀₊conj)(random_r[i],random_z[i]), -sqrt(1-random_z[i]^2)* (₊₂ψ₂₂₀₊conj(random_r[i],random_z[i]+δ)- ₊₂ψ₂₂₀₊conj(random_r[i],random_z[i]))/δ, rtol=1e-5)
        @assert isapprox(∂θ(₊₂ψ₂₂₀₋)(random_r[i],random_z[i]), -sqrt(1-random_z[i]^2)* (₊₂ψ₂₂₀₋(random_r[i],random_z[i]+δ)- ₊₂ψ₂₂₀₋(random_r[i],random_z[i]))/δ, rtol=1e-5)
        @assert isapprox(∂θ(₊₂ψ₂₂₀₋conj)(random_r[i],random_z[i]), -sqrt(1-random_z[i]^2)* (₊₂ψ₂₂₀₋conj(random_r[i],random_z[i]+δ)- ₊₂ψ₂₂₀₋conj(random_r[i],random_z[i]))/δ, rtol=1e-5)
    end
end