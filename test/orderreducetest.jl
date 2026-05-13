using ContourIntegrals
using KerrQNMShifts
using KerrQuasinormalModes
using Test
using CSV, DataFrames

println("Done usings")

@testset "ScalarKerrPertToSchw" begin   #CHECKED AND GOOD

    # dwOfile = "./OperatorShifts/OrderReducedShifts/dwscalarOcoefficients.csv"
    # dϵOfile = "./OperatorShifts/OrderReducedShifts/δscalarOcoefficients.csv"
    dwOfile = "./OperatorShifts/OrderReducedShifts/NPdwscalarOcoefficients.csv"
    dϵOfile = "./OperatorShifts/OrderReducedShifts/NPδscalarOcoefficients.csv"

    ∂ϵO = OperatorShift(dϵOfile)
    ∂ωO = OperatorShift(dwOfile)

    println("Made operator shifts")
    
    ψ = qnmfunctionnew(0,10,10,0,0.)

    # Compile ψ
    ψ(1,.5)
    println("Past ψ compile")

    ## Define the useful contours
    r₊ = ψ.R.r₊ ; r₋ = ψ.R.r₋ ; s = ψ.s ; Δr = 0.1*(r₊-r₋); ϵ = eps(0.1);

    # Define the Weight
    Σ = let a= ψ.a
        (r,z) -> Complex(r)^2+a^2*z^2
    end
    Δ = let a= ψ.a
        (r,z) -> Complex(r)^2+a^2-2*r
    end
    weight = let a= ψ.a
        (r,z) ->Σ(r,z)
    end

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

    # Make the operators
    ∂ϵOm2 = OperatorSandwich(ψ,∂ϵO,weight,ψ).Op
    ∂ωOm2 = OperatorSandwich(ψ,∂ωO,weight,ψ).Op
    println("Made Operators")

    ∂ϵ𝒪m2 = Integrate(∂ϵOm2, TheContourup,abstol=1e-10)[1]
    ∂ω𝒪m2 = Integrate(∂ωOm2, TheContourup,abstol=1e-10)[1]

    @show δω = -(∂ϵ𝒪m2/∂ω𝒪m2)
end
