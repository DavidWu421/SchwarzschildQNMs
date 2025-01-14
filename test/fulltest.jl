using ContourIntegrals
using KerrQNMShifts
using KerrQuasinormalModes
using Test
using Zygote

println("Done usings")

@testset "FullTest" begin
    pert_a=.1
    
    ψ = qnmfunctionnew(-2,2,2,0,0.)
    ψconj = qnmfunctionnew(-2,2,2,0,0.,is_conjugate=true)

    ψm = qnmfunctionnew(-2,2,2,0,0.,is_minus=true)
    ψmconj = qnmfunctionnew(-2,2,2,0,0.,is_minus=true,is_conjugate=true)

    ψup = qnmfunctionnew(-2,5,3,0,0.)
    ψdown = qnmfunctionnew(-2,5,3,0,0.,is_minus=true)


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

    weight = let a= ψ.a
        (r,z) ->8*ζ(r,z)^(4)*Σ(r,z)/((Δ(r,z))^2)
    end

    println("Past weight")

    ## Define the useful contours
    r₊ = ψ.R.r₊ ; r₋ = ψ.R.r₋ ; s = ψ.s ; Δr = 0.1*(r₊-r₋); ϵ = eps(0.1);

    #The upwards pointing contour
    point1up = r₊ + Δr - Δr*im
    point2up = r₊ - Δr - Δr*im

    radial1up = SemiInfiniteLine(point1up , point1up + Δr*im , false)
    angular = LineSegment(-1.0+100*ϵ , 1.0-100*ϵ , true) #to avoid the NaNs at the edges
    C1up = radial1up ⊗ angular

    radial2up = LineSegment(point1up,point2up,true)
    C2up = radial2up ⊗ angular

    radial3up = SemiInfiniteLine(point2up , point2up + Δr*im , true)
    C3up = radial3up ⊗ angular

    TheContourup = C1up⊕C2up⊕C3up

    #The downwards pointing contour
    point1down = r₊ + Δr + Δr*im
    point2down = r₊ - Δr + Δr*im

    radial1down = SemiInfiniteLine(point1down , point1down - Δr*im , false)
    angular = LineSegment(-1.0+100*ϵ , 1.0-100*ϵ , true) #to avoid the NaNs at the edges
    C1down = radial1down ⊗ angular

    radial2down = LineSegment(point1down,point2down,true)
    C2down = radial2down ⊗ angular

    radial3down = SemiInfiniteLine(point2down , point2down - Δr*im , true)
    C3down = radial3down ⊗ angular

    TheContourdown = C1down⊕C2down⊕C3down
    println("Done Contours")


    Ofile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/FullTest/Ocoefficients.csv"
    dwOfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/FullTest/dwOcoefficients.csv"
    Hfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/FullTest/Hcoefficients.csv"
    Ifile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/FullTest/Icoefficients.csv"

    O = OperatorShift(Ofile)
    ∂ωO = OperatorShift(dwOfile)
    H = OperatorShift(Hfile)
    I = OperatorShift(Ifile)

    println("Made operator shifts")
    
    OSchw1 = OperatorSandwich(ψ,O,weight,ψ).Op
    OSchw2 = OperatorSandwich(ψm,O,weight,ψm).Op
    ∂ωOplusSchw = OperatorSandwich(ψ,∂ωO,weight,ψ).Op
    ∂ωOminusSchw = OperatorSandwich(ψm,∂ωO,weight,ψm).Op
    HplusSchw = OperatorSandwich(ψ,H,weight,ψ).Op
    HminusSchw = OperatorSandwich(ψm,H,weight,ψm).Op
    IplusSchw = OperatorSandwich(ψ,I,weight,ψmconj).Op
    IminusSchw = OperatorSandwich(ψm,I,weight,ψconj).Op

    println("Made Operators")

    @show OSchw1(8+im,0.6,pertparam=pert_a)
    @show OSchw1(8+im,0.6,pertparam=pert_a)
    @show ∂ωOplusSchw(8+im,0.6,pertparam=pert_a)
    @show ∂ωOminusSchw(8+im,0.6,pertparam=pert_a)
    @show HplusSchw(8+im,0.6,pertparam=pert_a)
    @show HminusSchw(8+im,0.6,pertparam=pert_a)
    @show IplusSchw(8+im,0.6,pertparam=pert_a)
    @show IminusSchw(8+im,0.6,pertparam=pert_a)

    println("Complied Operators")

    𝒪plusSchw= Integrate(OplusSchw, TheContourup,pertparam=pert_a,abstol=1e-6)[1]
    𝒪minusSchw= Integrate(OminusSchw, TheContourdown,pertparam=pert_a,abstol=1e-6)[1]
    ∂ω𝒪plusSchw= Integrate(∂ωOplusSchw, TheContourup,pertparam=pert_a,abstol=1e-6)[1]
    ∂ω𝒪minusSchw= conj(Integrate(∂ωOminusSchw, TheContourdown,pertparam=pert_a,abstol=1e-6)[1])
    ℋplusSchw= Integrate(HplusSchw, TheContourup,pertparam=pert_a,abstol=1e-6)[1]
    ℋminusSchw= conj(Integrate(HminusSchw, TheContourdown,pertparam=pert_a,abstol=1e-6)[1])
    ℐplusSchw= Integrate(IplusSchw, TheContourup,pertparam=pert_a,abstol=1e-6)[1]
    ℐminusSchw= conj(Integrate(IminusSchw, TheContourdown,pertparam=pert_a,abstol=1e-6)[1])

    @show ∂ω𝒪plusSchw
    @show ∂ω𝒪minusSchw
    @show ℋplusSchw
    @show ℋminusSchw
    @show ℐplusSchw
    @show ℐminusSchw 
    
    ω2s=Computeω2(∂ω𝒪plusSchw,∂ω𝒪minusSchw,ℋplusSchw,ℋminusSchw,ℐplusSchw,ℐminusSchw,ψ)
    @show ω2s

end


@testset "Whole thing" begin
    println("Started KerrAsPertub Test: ")
    pert_a=.1
    
    ψ = qnmfunctionnew(-2,2,2,0,0.)
    ψconj = qnmfunctionnew(-2,2,2,0,0.,is_conjugate=true)

    ψm = qnmfunctionnew(-2,2,2,0,0.,is_minus=true)
    ψmconj = qnmfunctionnew(-2,2,2,0,0.,is_minus=true,is_conjugate=true)

    ψup = qnmfunctionnew(-2,5,3,0,0.)
    ψdown = qnmfunctionnew(-2,5,3,0,0.,is_minus=true)


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

    weight = let a= ψ.a
        (r,z) ->8*ζ(r,z)^(4)*Σ(r,z)/((Δ(r,z))^2)
    end

    println("Past weight")

    ## Define the useful contours
    r₊ = ψ.R.r₊ ; r₋ = ψ.R.r₋ ; s = ψ.s ; Δr = 0.1*(r₊-r₋); ϵ = eps(0.1);

    #The upwards pointing contour
    point1up = r₊ + Δr - Δr*im
    point2up = r₊ - Δr - Δr*im

    radial1up = SemiInfiniteLine(point1up , point1up + Δr*im , false)
    angular = LineSegment(-1.0+100*ϵ , 1.0-100*ϵ , true) #to avoid the NaNs at the edges
    C1up = radial1up ⊗ angular

    radial2up = LineSegment(point1up,point2up,true)
    C2up = radial2up ⊗ angular

    radial3up = SemiInfiniteLine(point2up , point2up + Δr*im , true)
    C3up = radial3up ⊗ angular

    TheContourup = C1up⊕C2up⊕C3up

    #The downwards pointing contour
    point1down = r₊ + Δr + Δr*im
    point2down = r₊ - Δr + Δr*im

    radial1down = SemiInfiniteLine(point1down , point1down - Δr*im , false)
    angular = LineSegment(-1.0+100*ϵ , 1.0-100*ϵ , true) #to avoid the NaNs at the edges
    C1down = radial1down ⊗ angular

    radial2down = LineSegment(point1down,point2down,true)
    C2down = radial2down ⊗ angular

    radial3down = SemiInfiniteLine(point2down , point2down - Δr*im , true)
    C3down = radial3down ⊗ angular

    TheContourdown = C1down⊕C2down⊕C3down
    println("Done Contours")


    Oplusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/HermTestOpluscoefficients.csv"
    Ominusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/HermTestOminuscoefficients2.csv"
    dwOfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/dwOpluscoefficients.csv"
    Hfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/Hpluscoefficients.csv"
    Ifile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/Ipluscoefficients.csv"

    Oplus = OperatorShift(Oplusfile)
    Ominus  = OperatorShift(Ominusfile)
    ∂ωO = OperatorShift(dwOfile)
    H = OperatorShift(Hfile)
    I = OperatorShift(Ifile)

    println("Made operator shifts")
    
    OplusSchw = OperatorSandwich(ψ,Oplus,weight,ψup).Op
    OminusSchw = OperatorSandwich(ψm,Ominus,weight,ψdown).Op
    ∂ωOplusSchw = OperatorSandwich(ψ,∂ωO,weight,ψ).Op
    ∂ωOminusSchw = OperatorSandwich(ψm,∂ωO,weight,ψm).Op
    HplusSchw = OperatorSandwich(ψ,H,weight,ψ).Op
    HminusSchw = OperatorSandwich(ψm,H,weight,ψm).Op
    IplusSchw = OperatorSandwich(ψ,I,weight,ψmconj).Op
    IminusSchw = OperatorSandwich(ψm,I,weight,ψconj).Op

    println("Made Operators")

    @show OplusSchw(8+im,0.6,pertparam=pert_a)
    @show OminusSchw(8+im,0.6,pertparam=pert_a)
    @show ∂ωOplusSchw(8+im,0.6,pertparam=pert_a)
    @show ∂ωOminusSchw(8+im,0.6,pertparam=pert_a)
    @show HplusSchw(8+im,0.6,pertparam=pert_a)
    @show HminusSchw(8+im,0.6,pertparam=pert_a)
    @show IplusSchw(8+im,0.6,pertparam=pert_a)
    @show IminusSchw(8+im,0.6,pertparam=pert_a)

    println("Complied Operators")

    𝒪plusSchw= Integrate(OplusSchw, TheContourup,pertparam=pert_a,abstol=1e-6)[1]
    𝒪minusSchw= Integrate(OminusSchw, TheContourdown,pertparam=pert_a,abstol=1e-6)[1]
    ∂ω𝒪plusSchw= Integrate(∂ωOplusSchw, TheContourup,pertparam=pert_a,abstol=1e-6)[1]
    ∂ω𝒪minusSchw= conj(Integrate(∂ωOminusSchw, TheContourdown,pertparam=pert_a,abstol=1e-6)[1])
    ℋplusSchw= Integrate(HplusSchw, TheContourup,pertparam=pert_a,abstol=1e-6)[1]
    ℋminusSchw= conj(Integrate(HminusSchw, TheContourdown,pertparam=pert_a,abstol=1e-6)[1])
    ℐplusSchw= Integrate(IplusSchw, TheContourup,pertparam=pert_a,abstol=1e-6)[1]
    ℐminusSchw= conj(Integrate(IminusSchw, TheContourdown,pertparam=pert_a,abstol=1e-6)[1])

    @show ∂ω𝒪plusSchw
    @show ∂ω𝒪minusSchw
    @show ℋplusSchw
    @show ℋminusSchw
    @show ℐplusSchw
    @show ℐminusSchw 
    
    ω2s=Computeω2(∂ω𝒪plusSchw,∂ω𝒪minusSchw,ℋplusSchw,ℋminusSchw,ℐplusSchw,ℐminusSchw,ψ)
    @show ω2s

end


@testset "derivative test" begin
    pert_a=.1
    
    ψ = qnmfunctionnew(-2,2,2,0,0.1)
    ψm = qnmfunctionnew(-2,2,2,0,0.1,is_minus=true)
    ψsm = qnmfunctionnew(2,2,2,0,0.1,is_minus=true)
    ψconj = qnmfunctionnew(-2,2,2,0,0.1,is_conjugate=true)
    ψmconj = qnmfunctionnew(-2,2,2,0,0.1,is_minus=true,is_conjugate=true)

    ψ(3)
    ψm(3)
    ψconj(3)
    ψmconj(3)
    
    ψ(3+im)
    ψm(3+im)
    ψconj(3+im)
    ψmconj(3+im)

    ψ.S(.4)
    ψsm.S(.4)
    ψm.S(.4)
    ψconj.S(.4)
    ψmconj.S(.4)



    ∂r(ψ)(3+im)-(ψ(3.001+1.001im)-ψ(3+im))/(.001+.001im)
    ∂r(ψ)(3)-(ψ(3.001)-ψ(3))/(.001)

    ∂r(ψm)(3+im)-(ψm(3.001+1.001im)-ψm(3+im))/(.001+.001im)
    ∂r(ψm)(3)-(ψm(3.001)-ψm(3))/(.001)

    ∂r(ψconj)(3+im)-(ψconj(3.001+1.001im)-ψconj(3+im))/(.001+.001im)
    ∂r(ψconj)(3)-(ψconj(3.001)-ψconj(3))/(.001)

    ∂r(ψmconj)(3+im)-(ψmconj(3.001+1.001im)-ψmconj(3+im))/(.001+.001im)
    ∂r(ψmconj)(3)-(ψmconj(3.001)-ψmconj(3))/(.001)


    ∂θ(ψ.S)(.4)-(ψ.S(.4001)-ψ.S(.4))/.0001*(-sqrt(1-.4^2))

    ∂θ(ψm.S)(.4)-(ψm.S(.4001)-ψm.S(.4))/.0001*(-sqrt(1-.4^2))

    ∂θ(ψconj.S)(.4)-(ψconj.S(.4001)-ψconj.S(.4))/.0001*(-sqrt(1-.4^2))

    ∂θ(ψmconj.S)(.4)-(ψmconj.S(.4001)-ψmconj.S(.4))/.0001*(-sqrt(1-.4^2))
end

@testtset "DirectionFlip" begin
    pert_a=.1
    
    ψ = qnmfunctionnew(-2,2,2,0,0.)

    ψm = qnmfunctionnew(-2,2,2,0,0.,is_minus=true)

    ψ1 = qnmfunctionnew(-2,4,2,0,0.1)
    ψm1 = qnmfunctionnew(-2,4,2,0,0.1,is_minus=true)

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

    weight = let a= ψ.a
        (r,z) ->8*ζ(r,z)^(4)*Σ(r,z)/((Δ(r,z))^2)
    end

    println("Past weight")

    ## Define the useful contours
    r₊ = ψ.R.r₊ ; r₋ = ψ.R.r₋ ; s = ψ.s ; Δr = 0.1*(r₊-r₋); ϵ = eps(0.1);

    #The upwards pointing contour
    point1up = r₊ + Δr - Δr*im
    point2up = r₊ - Δr - Δr*im

    radial1up = SemiInfiniteLine(point1up , point1up + Δr*im , false)
    angular = LineSegment(-1.0+100*ϵ , 1.0-100*ϵ , true) #to avoid the NaNs at the edges
    C1up = radial1up ⊗ angular

    radial2up = LineSegment(point1up,point2up,true)
    C2up = radial2up ⊗ angular

    radial3up = SemiInfiniteLine(point2up , point2up + Δr*im , true)
    C3up = radial3up ⊗ angular

    TheContourup = C1up⊕C2up⊕C3up

    #The downwards pointing contour
    point1down = r₊ + Δr + Δr*im
    point2down = r₊ - Δr + Δr*im

    radial1down = SemiInfiniteLine(point1down , point1down - Δr*im , false)
    angular = LineSegment(-1.0+100*ϵ , 1.0-100*ϵ , true) #to avoid the NaNs at the edges
    C1down = radial1down ⊗ angular

    radial2down = LineSegment(point1down,point2down,true)
    C2down = radial2down ⊗ angular

    radial3down = SemiInfiniteLine(point2down , point2down - Δr*im , true)
    C3down = radial3down ⊗ angular

    TheContourdown = C1down⊕C2down⊕C3down
    println("Done Contours")

    Oplusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/HermTestOpluscoefficients.csv"
    Ominusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/HermTestOminuscoefficients2.csv"

    Oplus=OperatorShift(Oplusfile)
    Ominus=OperatorShift(Ominusfile)

    println("Made operator shifts")
    
    OplusSchw0check = OperatorSandwich(ψ1,Oplus,weight,ψ).Op
    OminusSchw0check = OperatorSandwich(ψm1,Ominus,weight,ψm).Op

    
    OplusSchw1 = OperatorSandwich(ψ,Oplus,weight,ψ1).Op
    OminusSchw1 = OperatorSandwich(ψm,Ominus,weight,ψm1).Op

    
    println("Made Operators")

    @show OplusSchw0check(8+im,.6,pertparam=pert_a)
    @show OminusSchw0check(8+im,.6,pertparam=pert_a)

    @show OplusSchw1(8+im,0.6,pertparam=pert_a)
    @show OminusSchw1(8+im,0.6,pertparam=pert_a)

    println("Complied Operators")

    𝒪plusSchw1= Integrate(OplusSchw1, TheContourup,pertparam=pert_a,abstol=1e-6)[1]
    𝒪minusSchw1= Integrate(OminusSchw1, TheContourdown,pertparam=pert_a,abstol=1e-6)[1]

    @show 𝒪plusSchw1
    @show 𝒪minusSchw1
end

@testset "HermiticityTest" begin
    println("Started HermiticityTest: ")
    pert_a=.1
    
    ψ = qnmfunctionnew(-2,2,2,0,0.)

    ψm = qnmfunctionnew(-2,2,2,0,0.,is_minus=true)

    ψ1 = qnmfunctionnew(-2,4,2,0,0.1)
    ψm1 = qnmfunctionnew(-2,4,2,0,0.1,is_minus=true)
    ψ2 = qnmfunctionnew(-2,4,2,0,0.)
    ψm2 = qnmfunctionnew(-2,4,2,0,0.,is_minus=true)
    ψ3 = qnmfunctionnew(-2,5,2,0,0.)
    ψm3 = qnmfunctionnew(-2,5,2,0,0.,is_minus=true)

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
    
    OplusSchw0check = OperatorSandwich(ψm,Oplus,weightplus,ψ).Op
    OminusSchw0check = OperatorSandwich(ψ,Ominus,weightminus,ψm).Op

    
    OplusSchw1 = OperatorSandwich(ψ,Oplus,weightplus,ψ1).Op
    OminusSchw1 = OperatorSandwich(ψm,Ominus,weightminus,ψm1).Op
    OplusSchw2 = OperatorSandwich(ψ,Oplus,weightplus,ψ2).Op
    OminusSchw2 = OperatorSandwich(ψm,Ominus,weightminus,ψm2).Op
    OplusSchw3 = OperatorSandwich(ψ,Oplus,weightplus,ψ3).Op
    OminusSchw3 = OperatorSandwich(ψm,Ominus,weightminus,ψm3).Op

    
    println("Made Operators")

    @show OplusSchw0check(8+im,.6,pertparam=pert_a)
    @show OminusSchw0check(8+im,.6,pertparam=pert_a)

    @show OplusSchw1(8+im,0.6,pertparam=pert_a)
    @show OminusSchw1(8+im,0.6,pertparam=pert_a)
    @show OplusSchw2(8+im,0.6,pertparam=pert_a)
    @show OminusSchw2(8+im,0.6,pertparam=pert_a)
    @show OplusSchw3(8+im,0.6,pertparam=pert_a)
    @show OminusSchw3(8+im,0.6,pertparam=pert_a)

    println("Complied Operators")

    𝒪plusSchw1= Integrate(OplusSchw1, TheContour,pertparam=pert_a,abstol=1e-6)[1]
    𝒪minusSchw1= Integrate(OminusSchw1, TheContour,pertparam=pert_a,abstol=1e-6)[1]
    𝒪plusSchw2= Integrate(OplusSchw2, TheContour,pertparam=pert_a,abstol=1e-6)[1]
    𝒪minusSchw2= Integrate(OminusSchw2, TheContour,pertparam=pert_a,abstol=1e-6)[1]
    𝒪plusSchw3= Integrate(OplusSchw3, TheContour,pertparam=pert_a,abstol=1e-6)[1]
    𝒪minusSchw3= Integrate(OminusSchw3, TheContour,pertparam=pert_a,abstol=1e-6)[1]

    @show 𝒪plusSchw1
    @show 𝒪minusSchw1
    @show 𝒪plusSchw2
    @show 𝒪minusSchw2
    @show 𝒪plusSchw3
    @show 𝒪minusSchw3
end

@testset "KerrAsPerturb" begin

    println("Started KerrAsPertub Test: ")
    pert_a=.1
    
    ψ = qnmfunctionnew(-2,2,2,0,0.)
    ψconj = qnmfunctionnew(-2,2,2,0,0.,is_conjugate=true)

    ψm = qnmfunctionnew(-2,2,2,0,0.,is_minus=true)
    ψmconj = qnmfunctionnew(-2,2,2,0,0.,is_minus=true,is_conjugate=true)

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
    ₋₂ψ₂₂₀₋ = qnmfunctionnew(-2,2,2,0,0.5,is_minus=true)
    ₊₂ψ₂₂₀₊ = qnmfunctionnew(2,2,2,0,0.5)
    ₊₂ψ₂₂₀₋= qnmfunctionnew(2,2,2,0,0.5,is_minus=true)
    ₋₂ψ₂₂₀₊conj = qnmfunctionnew(-2,2,2,0,0.5,is_conjugate=true) 
    ₋₂ψ₂₂₀₋conj = qnmfunctionnew(-2,2,2,0,0.5,is_minus=true,is_conjugate=true)
    ₊₂ψ₂₂₀₊conj = qnmfunctionnew(2,2,2,0,0.5,is_conjugate=true)
    ₊₂ψ₂₂₀₋conj= qnmfunctionnew(2,2,2,0,0.5,is_minus=true,is_conjugate=true)
    r₊=₋₂ψ₂₂₀₊.R.r₊
    r₋=₋₂ψ₂₂₀₊.R.r₋

    # random_r = random_complex(10)
    random_r= rand(1:10,10)
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
    @assert abs(₋₂ψ₂₂₀₋(r₊+.1-100*im)) <= 10^(-8)
    @assert abs(₊₂ψ₂₂₀₋(r₊+.1-100*im)) <= 10^(-8)
    @assert abs(₋₂ψ₂₂₀₊conj(r₊+.1-100*im)) <= 10^(-8)
    @assert abs(₊₂ψ₂₂₀₊conj(r₊+.1-100*im)) <= 10^(-8)
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
        # @assert isapprox(∂r(₋₂ψ₂₂₀₊)(random_r[i]), (₋₂ψ₂₂₀₊(random_r[i]+δ)- ₋₂ψ₂₂₀₊(random_r[i]))/δ, rtol=1e-5)
        # @assert isapprox(∂r(₋₂ψ₂₂₀₊conj)(random_r[i]), (₋₂ψ₂₂₀₊conj(random_r[i]+δ)- ₋₂ψ₂₂₀₊conj(random_r[i]))/δ, rtol=1e-5)
        # @assert isapprox(∂r(₋₂ψ₂₂₀₋)(random_r[i]), (₋₂ψ₂₂₀₋(random_r[i]+δ)- ₋₂ψ₂₂₀₋(random_r[i]))/δ, rtol=1e-5)
        # @assert isapprox(∂r(₋₂ψ₂₂₀₋conj)(random_r[i]), (₋₂ψ₂₂₀₋conj(random_r[i]+δ)- ₋₂ψ₂₂₀₋conj(random_r[i]))/δ, rtol=1e-5)
        # @assert isapprox(∂r(₊₂ψ₂₂₀₊)(random_r[i]), (₊₂ψ₂₂₀₊(random_r[i]+δ)- ₊₂ψ₂₂₀₊(random_r[i]))/δ, rtol=1e-5)
        # @assert isapprox(∂r(₊₂ψ₂₂₀₊conj)(random_r[i]), (₊₂ψ₂₂₀₊conj(random_r[i]+δ)- ₊₂ψ₂₂₀₊conj(random_r[i]))/δ, rtol=1e-5)
        # @assert isapprox(∂r(₊₂ψ₂₂₀₋)(random_r[i]), (₊₂ψ₂₂₀₋(random_r[i]+δ)- ₊₂ψ₂₂₀₋(random_r[i]))/δ, rtol=1e-5)
        # @assert isapprox(∂r(₊₂ψ₂₂₀₋conj)(random_r[i]), (₊₂ψ₂₂₀₋conj(random_r[i]+δ)- ₊₂ψ₂₂₀₋conj(random_r[i]))/δ, rtol=1e-5)

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