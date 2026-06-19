using ContourIntegrals
using KerrQNMShifts
using KerrQuasinormalModes
using Test
using CSV, DataFrames

println("Done usings")

@testset "FullJPTake2" begin 

    dwOfile = "./OperatorShifts/FullNPTest/FullG2Take2/dwOOps/dwOcoefficients.csv"

    H1Afile = "./OperatorShifts/FullNPTest/FullG2Take2/HOps/H1A.csv"
    H1Bfile = "./OperatorShifts/FullNPTest/FullG2Take2/HOps/H1B.csv"
    H2Afile = "./OperatorShifts/FullNPTest/FullG2Take2/HOps/H2A.csv"
    H2Bfile = "./OperatorShifts/FullNPTest/FullG2Take2/HOps/H2B.csv"
    H3Afile = "./OperatorShifts/FullNPTest/FullG2Take2/HOps/H3A.csv"
    H3Bfile = "./OperatorShifts/FullNPTest/FullG2Take2/HOps/H3B.csv"
    H4Afile = "./OperatorShifts/FullNPTest/FullG2Take2/HOps/H4A.csv"
    H4Bfile = "./OperatorShifts/FullNPTest/FullG2Take2/HOps/H4B.csv"
    H10file = "./OperatorShifts/FullNPTest/FullG2Take2/HOps/H10.csv"

    I1Afile = "./OperatorShifts/FullNPTest/FullG2Take2/IOps/I1A.csv"
    I1Bfile = "./OperatorShifts/FullNPTest/FullG2Take2/IOps/I1B.csv"
    I2Afile = "./OperatorShifts/FullNPTest/FullG2Take2/IOps/I2A.csv"
    I2Bfile = "./OperatorShifts/FullNPTest/FullG2Take2/IOps/I2B.csv"
    I3Afile = "./OperatorShifts/FullNPTest/FullG2Take2/IOps/I3A.csv"
    I3Bfile = "./OperatorShifts/FullNPTest/FullG2Take2/IOps/I3B.csv"
    I4Afile = "./OperatorShifts/FullNPTest/FullG2Take2/IOps/I4A.csv"
    I4Bfile = "./OperatorShifts/FullNPTest/FullG2Take2/IOps/I4B.csv"
    I10file = "./OperatorShifts/FullNPTest/FullG2Take2/IOps/I10.csv"

    dwO = OperatorShift(dwOfile)

    H1A = OperatorShift(H1Afile)
    H1B = OperatorShift(H1Bfile)
    H2A = OperatorShift(H2Afile)
    H2B = OperatorShift(H2Bfile)
    H3A = OperatorShift(H3Afile)
    H3B = OperatorShift(H3Bfile)
    H4A = OperatorShift(H4Afile)
    H4B = OperatorShift(H4Bfile)
    H10 = OperatorShift(H10file)

    I1A = OperatorShift(I1Afile)
    I1B = OperatorShift(I1Bfile)
    I2A = OperatorShift(I2Afile)
    I2B = OperatorShift(I2Bfile)
    I3A = OperatorShift(I3Afile)
    I3B = OperatorShift(I3Bfile)
    I4A = OperatorShift(I4Afile)
    I4B = OperatorShift(I4Bfile)
    I10 = OperatorShift(I10file)

    println("Made operator shifts")
    
    ψ = qnmfunctionnew(-2,2,2,0,10^(-8))
    ψm = qnmfunctionnew(-2,2,2,0,10^(-8), is_minus=true)
    ψconj =  qnmfunctionnew(-2,2,2,0,10^(-8), is_conjugate=true)
    ψmconj = qnmfunctionnew(-2,2,2,0,10^(-8),is_conjugate=true,is_minus=true)

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
    ζ = let a= ψ.a
        (r,z) -> r-im*a*z
    end

    weight = let a= ψ.a
        (r,z) ->8*ζ(r,z)^(4)*Σ(r,z)/((Δ(r,z))^2)
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

    # Make the operators
    dwOplus = OperatorSandwich(ψ,dwO,weight,ψ).Op
    dwOminus = OperatorSandwich(ψm,dwO,weight,ψm).Op

    H1Aplus = OperatorSandwich(ψ,H1A,weight,ψ).Op
    H1Bplus = OperatorSandwich(ψ,H1B,weight,ψ).Op
    H2Aplus = OperatorSandwich(ψ,H2A,weight,ψ).Op
    H2Bplus = OperatorSandwich(ψ,H2B,weight,ψ).Op
    H3Aplus = OperatorSandwich(ψ,H3A,weight,ψ).Op
    H3Bplus = OperatorSandwich(ψ,H3B,weight,ψ).Op
    H4Aplus = OperatorSandwich(ψ,H4A,weight,ψ).Op
    H4Bplus = OperatorSandwich(ψ,H4B,weight,ψ).Op
    H10plus = OperatorSandwich(ψ,H10,weight,ψ).Op

    H1Aminus = OperatorSandwich(ψm,H1A,weight,ψm).Op
    H1Bminus = OperatorSandwich(ψm,H1B,weight,ψm).Op
    H2Aminus = OperatorSandwich(ψm,H2A,weight,ψm).Op
    H2Bminus = OperatorSandwich(ψm,H2B,weight,ψm).Op
    H3Aminus = OperatorSandwich(ψm,H3A,weight,ψm).Op
    H3Bminus = OperatorSandwich(ψm,H3B,weight,ψm).Op
    H4Aminus = OperatorSandwich(ψm,H4A,weight,ψm).Op
    H4Bminus = OperatorSandwich(ψm,H4B,weight,ψm).Op
    H10minus = OperatorSandwich(ψm,H10,weight,ψm).Op

    I1Aplus = OperatorSandwich(ψ,I1A,weight,ψmconj).Op
    I1Bplus = OperatorSandwich(ψ,I1B,weight,ψmconj).Op
    I2Aplus = OperatorSandwich(ψ,I2A,weight,ψmconj).Op
    I2Bplus = OperatorSandwich(ψ,I2B,weight,ψmconj).Op
    I3Aplus = OperatorSandwich(ψ,I3A,weight,ψmconj).Op
    I3Bplus = OperatorSandwich(ψ,I3B,weight,ψmconj).Op
    I4Aplus = OperatorSandwich(ψ,I4A,weight,ψmconj).Op
    I4Bplus = OperatorSandwich(ψ,I4B,weight,ψmconj).Op
    I10plus = OperatorSandwich(ψ,H10,weight,ψmconj).Op

    I1Aminus = OperatorSandwich(ψm,I1A,weight,ψconj).Op
    I1Bminus = OperatorSandwich(ψm,I1B,weight,ψconj).Op
    I2Aminus = OperatorSandwich(ψm,I2A,weight,ψconj).Op
    I2Bminus = OperatorSandwich(ψm,I2B,weight,ψconj).Op
    I3Aminus = OperatorSandwich(ψm,I3A,weight,ψconj).Op
    I3Bminus = OperatorSandwich(ψm,I3B,weight,ψconj).Op
    I4Aminus = OperatorSandwich(ψm,I4A,weight,ψconj).Op
    I4Bminus = OperatorSandwich(ψm,I4B,weight,ψconj).Op
    I10minus = OperatorSandwich(ψm,H10,weight,ψconj).Op

    println("Made Operators")

    # Upwards-facing contours

    ∂ω𝒪plus = Integrate(dwOplus, TheContourup,abstol=1e-4)[1]

    ℋ1Aplus = Integrate(H1Aplus, TheContourup,abstol=1e-4)[1]
    ℋ1Bplus = Integrate(H1Bplus, TheContourup,abstol=1e-4)[1]
    ℋ2Aplus = Integrate(H2Aplus, TheContourup,abstol=1e-4)[1]
    ℋ2Bplus = Integrate(H2Bplus, TheContourup,abstol=1e-4)[1]
    ℋ3Aplus = Integrate(H3Aplus, TheContourup,abstol=1e-4)[1]
    ℋ3Bplus = Integrate(H3Bplus, TheContourup,abstol=1e-4)[1]
    ℋ4Aplus = Integrate(H4Aplus, TheContourup,abstol=1e-4)[1]
    ℋ4Bplus = Integrate(H4Bplus, TheContourup,abstol=1e-4)[1]
    ℋ10plus = Integrate(H10plus, TheContourup,abstol=1e-4)[1]

    ℐ1Aplus = Integrate(I1Aplus, TheContourup,abstol=1e-4)[1]
    ℐ1Bplus = Integrate(I1Bplus, TheContourup,abstol=1e-4)[1]
    ℐ2Aplus = Integrate(I2Aplus, TheContourup,abstol=1e-4)[1]
    ℐ2Bplus = Integrate(I2Bplus, TheContourup,abstol=1e-4)[1]
    ℐ3Aplus = Integrate(I3Aplus, TheContourup,abstol=1e-4)[1]
    ℐ3Bplus = Integrate(I3Bplus, TheContourup,abstol=1e-4)[1]
    ℐ4Aplus = Integrate(I4Aplus, TheContourup,abstol=1e-4)[1]
    ℐ4Bplus = Integrate(I4Bplus, TheContourup,abstol=1e-4)[1]
    ℐ10plus = Integrate(I10plus, TheContourup,abstol=1e-4)[1]

    # Downwards-facing contours

    ∂ω𝒪minus = Integrate(dwOminus, TheContourdown,abstol=1e-4)[1]

    ℋ1Aminus = Integrate(H1Aminus, TheContourdown,abstol=1e-4)[1]
    ℋ1Bminus = Integrate(H1Bminus, TheContourdown,abstol=1e-4)[1]
    ℋ2Aminus = Integrate(H2Aminus, TheContourdown,abstol=1e-4)[1]
    ℋ2Bminus = Integrate(H2Bminus, TheContourdown,abstol=1e-4)[1]
    ℋ3Aminus = Integrate(H3Aminus, TheContourdown,abstol=1e-4)[1]
    ℋ3Bminus = Integrate(H3Bminus, TheContourdown,abstol=1e-4)[1]
    ℋ4Aminus = Integrate(H4Aminus, TheContourdown,abstol=1e-4)[1]
    ℋ4Bminus = Integrate(H4Bminus, TheContourdown,abstol=1e-4)[1]
    ℋ10minus = Integrate(H10minus, TheContourdown,abstol=1e-4)[1]

    ℐ1Aminus = Integrate(I1Aminus, TheContourdown,abstol=1e-4)[1]
    ℐ1Bminus = Integrate(I1Bminus, TheContourdown,abstol=1e-4)[1]
    ℐ2Aminus = Integrate(I2Aminus, TheContourdown,abstol=1e-4)[1]
    ℐ2Bminus = Integrate(I2Bminus, TheContourdown,abstol=1e-4)[1]
    ℐ3Aminus = Integrate(I3Aminus, TheContourdown,abstol=1e-4)[1]
    ℐ3Bminus = Integrate(I3Bminus, TheContourdown,abstol=1e-4)[1]
    ℐ4Aminus = Integrate(I4Aminus, TheContourdown,abstol=1e-4)[1]
    ℐ4Bminus = Integrate(I4Bminus, TheContourdown,abstol=1e-4)[1]
    ℐ10minus = Integrate(I10minus, TheContourdown,abstol=1e-4)[1]

    ℋplus = ℋ1Aplus+ℋ1Bplus+ℋ2Aplus+ℋ2Bplus+ℋ3Aplus+ℋ3Bplus+ℋ4Aplus+ℋ4Bplus+ℋ10plus
    ℐplus = ℐ1Aplus+ℐ1Bplus+ℐ2Aplus+ℐ2Bplus+ℐ3Aplus+ℐ3Bplus+ℐ4Aplus+ℐ4Bplus+ℐ10plus

    ℋminus = ℋ1Aminus+ℋ1Bminus+ℋ2Aminus+ℋ2Bminus+ℋ3Aminus+ℋ3Bminus+ℋ4Aminus+ℋ4Bminus+ℋ10minus
    ℐminus = ℐ1Aminus+ℐ1Bminus+ℐ2Aminus+ℐ2Bminus+ℐ3Aminus+ℐ3Bminus+ℐ4Aminus+ℐ4Bminus+ℐ10minus

    ω2s=Computeω2(∂ω𝒪plus,∂ω𝒪minus,ℋplus,ℋminus,ℐplus,ℐminus,ψ)
    @show ω2s
    
end

@testset "FullButWrongJP" begin 

    dwOfile = "./OperatorShifts/FullNPTest/FullG2/dwOOps/dwOcoefficients.csv"

    H1Afile = "./OperatorShifts/FullNPTest/FullG2/HOps/H1A.csv"
    H1Bfile = "./OperatorShifts/FullNPTest/FullG2/HOps/H1B.csv"
    H2Afile = "./OperatorShifts/FullNPTest/FullG2/HOps/H2A.csv"
    H2Bfile = "./OperatorShifts/FullNPTest/FullG2/HOps/H2B.csv"
    H3Afile = "./OperatorShifts/FullNPTest/FullG2/HOps/H3A.csv"
    H3Bfile = "./OperatorShifts/FullNPTest/FullG2/HOps/H3B.csv"
    H4Afile = "./OperatorShifts/FullNPTest/FullG2/HOps/H4A.csv"
    H4Bfile = "./OperatorShifts/FullNPTest/FullG2/HOps/H4B.csv"
    H5Afile = "./OperatorShifts/FullNPTest/FullG2/HOps/H5A.csv"
    H5Bfile = "./OperatorShifts/FullNPTest/FullG2/HOps/H5B.csv"
    H6Afile = "./OperatorShifts/FullNPTest/FullG2/HOps/H6A.csv"
    H6Bfile = "./OperatorShifts/FullNPTest/FullG2/HOps/H6B.csv"
    H7Afile = "./OperatorShifts/FullNPTest/FullG2/HOps/H7A.csv"
    H7Bfile = "./OperatorShifts/FullNPTest/FullG2/HOps/H7B.csv"
    H8Afile = "./OperatorShifts/FullNPTest/FullG2/HOps/H8A.csv"
    H8Bfile = "./OperatorShifts/FullNPTest/FullG2/HOps/H8B.csv"
    H9file = "./OperatorShifts/FullNPTest/FullG2/HOps/H9.csv"
    H10file = "./OperatorShifts/FullNPTest/FullG2/HOps/H10.csv"

    I1Afile = "./OperatorShifts/FullNPTest/FullG2/IOps/I1A.csv"
    I1Bfile = "./OperatorShifts/FullNPTest/FullG2/IOps/I1B.csv"
    I2Afile = "./OperatorShifts/FullNPTest/FullG2/IOps/I2A.csv"
    I2Bfile = "./OperatorShifts/FullNPTest/FullG2/IOps/I2B.csv"
    I3Afile = "./OperatorShifts/FullNPTest/FullG2/IOps/I3A.csv"
    I3Bfile = "./OperatorShifts/FullNPTest/FullG2/IOps/I3B.csv"
    I4Afile = "./OperatorShifts/FullNPTest/FullG2/IOps/I4A.csv"
    I4Bfile = "./OperatorShifts/FullNPTest/FullG2/IOps/I4B.csv"
    I5Afile = "./OperatorShifts/FullNPTest/FullG2/IOps/I5A.csv"
    I5Bfile = "./OperatorShifts/FullNPTest/FullG2/IOps/I5B.csv"
    I6Afile = "./OperatorShifts/FullNPTest/FullG2/IOps/I6A.csv"
    I6Bfile = "./OperatorShifts/FullNPTest/FullG2/IOps/I6B.csv"
    I7Afile = "./OperatorShifts/FullNPTest/FullG2/IOps/I7A.csv"
    I7Bfile = "./OperatorShifts/FullNPTest/FullG2/IOps/I7B.csv"
    I8Afile = "./OperatorShifts/FullNPTest/FullG2/IOps/I8A.csv"
    I8Bfile = "./OperatorShifts/FullNPTest/FullG2/IOps/I8B.csv"
    I9file = "./OperatorShifts/FullNPTest/FullG2/IOps/I9.csv"
    I10file = "./OperatorShifts/FullNPTest/FullG2/IOps/I10.csv"

    dwO = OperatorShift(dwOfile)

    H1A = OperatorShift(H1Afile)
    H1B = OperatorShift(H1Bfile)
    H2A = OperatorShift(H2Afile)
    H2B = OperatorShift(H2Bfile)
    H3A = OperatorShift(H3Afile)
    H3B = OperatorShift(H3Bfile)
    H4A = OperatorShift(H4Afile)
    H4B = OperatorShift(H4Bfile)
    # H5A = OperatorShift(H5Afile)
    # H5B = OperatorShift(H5Bfile)
    # H6A = OperatorShift(H6Afile)
    # H6B = OperatorShift(H6Bfile)
    # H7A = OperatorShift(H7Afile)
    # H7B = OperatorShift(H7Bfile)
    # H8A = OperatorShift(H8Afile)
    # H8B = OperatorShift(H8Bfile)
    # H9 = OperatorShift(H9file)
    H10 = OperatorShift(H10file)

    I1A = OperatorShift(I1Afile)
    I1B = OperatorShift(I1Bfile)
    I2A = OperatorShift(I2Afile)
    I2B = OperatorShift(I2Bfile)
    I3A = OperatorShift(I3Afile)
    I3B = OperatorShift(I3Bfile)
    I4A = OperatorShift(I4Afile)
    I4B = OperatorShift(I4Bfile)
    # I5A = OperatorShift(I5Afile)
    # I5B = OperatorShift(I5Bfile)
    # I6A = OperatorShift(I6Afile)
    # I6B = OperatorShift(I6Bfile)
    # I7A = OperatorShift(I7Afile)
    # I7B = OperatorShift(I7Bfile)
    # I8A = OperatorShift(I8Afile)
    # I8B = OperatorShift(I8Bfile)
    # I9 = OperatorShift(I9file)
    I10 = OperatorShift(I10file)

    println("Made operator shifts")
    
    ψ = qnmfunctionnew(-2,2,2,0,0.1)
    ψm = qnmfunctionnew(-2,2,2,0,0.1, is_minus=true)
    ψconj =  qnmfunctionnew(-2,2,2,0,0.1, is_conjugate=true)
    ψmconj = qnmfunctionnew(-2,2,2,0,0.1,is_conjugate=true,is_minus=true)

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

    # Make the operators
    dwOplus = OperatorSandwich(ψ,dwO,weight,ψ).Op
    dwOminus = OperatorSandwich(ψm,dwO,weight,ψm).Op

    H1Aplus = OperatorSandwich(ψ,H1A,weight,ψ).Op
    H1Bplus = OperatorSandwich(ψ,H1B,weight,ψ).Op
    H2Aplus = OperatorSandwich(ψ,H2A,weight,ψ).Op
    H2Bplus = OperatorSandwich(ψ,H2B,weight,ψ).Op
    H3Aplus = OperatorSandwich(ψ,H3A,weight,ψ).Op
    H3Bplus = OperatorSandwich(ψ,H3B,weight,ψ).Op
    H4Aplus = OperatorSandwich(ψ,H4A,weight,ψ).Op
    H4Bplus = OperatorSandwich(ψ,H4B,weight,ψ).Op
    # H5Aplus = OperatorSandwich(ψ,H5A,weight,ψ).Op
    # H5Bplus = OperatorSandwich(ψ,H5B,weight,ψ).Op
    # H6Aplus = OperatorSandwich(ψ,H6A,weight,ψ).Op
    # H6Bplus = OperatorSandwich(ψ,H6B,weight,ψ).Op
    # H7Aplus = OperatorSandwich(ψ,H7A,weight,ψ).Op
    # H7Bplus = OperatorSandwich(ψ,H7B,weight,ψ).Op
    # H8Aplus = OperatorSandwich(ψ,H8A,weight,ψ).Op
    # H8Bplus = OperatorSandwich(ψ,H8B,weight,ψ).Op
    # H9plus = OperatorSandwich(ψ,H9,weight,ψ).Op
    H10plus = OperatorSandwich(ψ,H10,weight,ψ).Op

    H1Aminus = OperatorSandwich(ψm,H1A,weight,ψm).Op
    H1Bminus = OperatorSandwich(ψm,H1B,weight,ψm).Op
    H2Aminus = OperatorSandwich(ψm,H2A,weight,ψm).Op
    H2Bminus = OperatorSandwich(ψm,H2B,weight,ψm).Op
    H3Aminus = OperatorSandwich(ψm,H3A,weight,ψm).Op
    H3Bminus = OperatorSandwich(ψm,H3B,weight,ψm).Op
    H4Aminus = OperatorSandwich(ψm,H4A,weight,ψm).Op
    H4Bminus = OperatorSandwich(ψm,H4B,weight,ψm).Op
    # H5Aminus = OperatorSandwich(ψm,H5A,weight,ψm).Op
    # H5Bminus = OperatorSandwich(ψm,H5B,weight,ψm).Op
    # H6Aminus = OperatorSandwich(ψm,H6A,weight,ψm).Op
    # H6Bminus = OperatorSandwich(ψm,H6B,weight,ψm).Op
    # H7Aminus = OperatorSandwich(ψm,H7A,weight,ψm).Op
    # H7Bminus = OperatorSandwich(ψm,H7B,weight,ψm).Op
    # H8Aminus = OperatorSandwich(ψm,H8A,weight,ψm).Op
    # H8Bminus = OperatorSandwich(ψm,H8B,weight,ψm).Op
    # H9minus = OperatorSandwich(ψm,H9,weight,ψm).Op
    H10minus = OperatorSandwich(ψm,H10,weight,ψm).Op

    I1Aplus = OperatorSandwich(ψ,I1A,weight,ψmconj).Op
    I1Bplus = OperatorSandwich(ψ,I1B,weight,ψmconj).Op
    I2Aplus = OperatorSandwich(ψ,I2A,weight,ψmconj).Op
    I2Bplus = OperatorSandwich(ψ,I2B,weight,ψmconj).Op
    I3Aplus = OperatorSandwich(ψ,I3A,weight,ψmconj).Op
    I3Bplus = OperatorSandwich(ψ,I3B,weight,ψmconj).Op
    I4Aplus = OperatorSandwich(ψ,I4A,weight,ψmconj).Op
    I4Bplus = OperatorSandwich(ψ,I4B,weight,ψmconj).Op
    # I5Aplus = OperatorSandwich(ψ,I5A,weight,ψmconj).Op
    # I5Bplus = OperatorSandwich(ψ,I5B,weight,ψmconj).Op
    # I6Aplus = OperatorSandwich(ψ,I6A,weight,ψmconj).Op
    # I6Bplus = OperatorSandwich(ψ,I6B,weight,ψmconj).Op
    # I7Aplus = OperatorSandwich(ψ,I7A,weight,ψmconj).Op
    # I7Bplus = OperatorSandwich(ψ,I7B,weight,ψmconj).Op
    # I8Aplus = OperatorSandwich(ψ,I8A,weight,ψmconj).Op
    # I8Bplus = OperatorSandwich(ψ,I8B,weight,ψmconj).Op
    # I9plus = OperatorSandwich(ψ,H9,weight,ψmconj).Op
    I10plus = OperatorSandwich(ψ,H10,weight,ψmconj).Op

    I1Aminus = OperatorSandwich(ψm,I1A,weight,ψconj).Op
    I1Bminus = OperatorSandwich(ψm,I1B,weight,ψconj).Op
    I2Aminus = OperatorSandwich(ψm,I2A,weight,ψconj).Op
    I2Bminus = OperatorSandwich(ψm,I2B,weight,ψconj).Op
    I3Aminus = OperatorSandwich(ψm,I3A,weight,ψconj).Op
    I3Bminus = OperatorSandwich(ψm,I3B,weight,ψconj).Op
    I4Aminus = OperatorSandwich(ψm,I4A,weight,ψconj).Op
    I4Bminus = OperatorSandwich(ψm,I4B,weight,ψconj).Op
    # I5Aminus = OperatorSandwich(ψm,I5A,weight,ψconj).Op
    # I5Bminus = OperatorSandwich(ψm,I5B,weight,ψconj).Op
    # I6Aminus = OperatorSandwich(ψm,I6A,weight,ψconj).Op
    # I6Bminus = OperatorSandwich(ψm,I6B,weight,ψconj).Op
    # I7Aminus = OperatorSandwich(ψm,I7A,weight,ψconj).Op
    # I7Bminus = OperatorSandwich(ψm,I7B,weight,ψconj).Op
    # I8Aminus = OperatorSandwich(ψm,I8A,weight,ψconj).Op
    # I8Bminus = OperatorSandwich(ψm,I8B,weight,ψconj).Op
    # I9minus = OperatorSandwich(ψm,H9,weight,ψconj).Op
    I10minus = OperatorSandwich(ψm,H10,weight,ψconj).Op

    println("Made Operators")

    # Upwards-facing contours

    ∂ω𝒪plus = Integrate(dwOplus, TheContourup,abstol=1e-4)[1]

    ℋ1Aplus = Integrate(H1Aplus, TheContourup,abstol=1e-4)[1]
    ℋ1Bplus = Integrate(H1Bplus, TheContourup,abstol=1e-4)[1]
    ℋ2Aplus = Integrate(H2Aplus, TheContourup,abstol=1e-4)[1]
    ℋ2Bplus = Integrate(H2Bplus, TheContourup,abstol=1e-4)[1]
    ℋ3Aplus = Integrate(H3Aplus, TheContourup,abstol=1e-4)[1]
    ℋ3Bplus = Integrate(H3Bplus, TheContourup,abstol=1e-4)[1]
    ℋ4Aplus = Integrate(H4Aplus, TheContourup,abstol=1e-4)[1]
    ℋ4Bplus = Integrate(H4Bplus, TheContourup,abstol=1e-4)[1]
    # ℋ5Aplus = Integrate(H5Aplus, TheContourup,abstol=1e-4)[1]
    # ℋ5Bplus = Integrate(H5Bplus, TheContourup,abstol=1e-4)[1]
    # ℋ6Aplus = Integrate(H6Aplus, TheContourup,abstol=1e-4)[1]
    # ℋ6Bplus = Integrate(H6Bplus, TheContourup,abstol=1e-4)[1]
    # ℋ7Aplus = Integrate(H7Aplus, TheContourup,abstol=1e-4)[1]
    # ℋ7Bplus = Integrate(H7Bplus, TheContourup,abstol=1e-4)[1]
    # ℋ8Aplus = Integrate(H8Aplus, TheContourup,abstol=1e-4)[1]
    # ℋ8Bplus = Integrate(H8Bplus, TheContourup,abstol=1e-4)[1]
    # ℋ9plus = Integrate(H9plus, TheContourup,abstol=1e-4)[1]
    ℋ10plus = Integrate(H10plus, TheContourup,abstol=1e-4)[1]

    ℐ1Aplus = Integrate(I1Aplus, TheContourup,abstol=1e-4)[1]
    ℐ1Bplus = Integrate(I1Bplus, TheContourup,abstol=1e-4)[1]
    ℐ2Aplus = Integrate(I2Aplus, TheContourup,abstol=1e-4)[1]
    ℐ2Bplus = Integrate(I2Bplus, TheContourup,abstol=1e-4)[1]
    ℐ3Aplus = Integrate(I3Aplus, TheContourup,abstol=1e-4)[1]
    ℐ3Bplus = Integrate(I3Bplus, TheContourup,abstol=1e-4)[1]
    ℐ4Aplus = Integrate(I4Aplus, TheContourup,abstol=1e-4)[1]
    ℐ4Bplus = Integrate(I4Bplus, TheContourup,abstol=1e-4)[1]
    # ℐ5Aplus = Integrate(I5Aplus, TheContourup,abstol=1e-4)[1]
    # ℐ5Bplus = Integrate(I5Bplus, TheContourup,abstol=1e-4)[1]
    # ℐ6Aplus = Integrate(I6Aplus, TheContourup,abstol=1e-4)[1]
    # ℐ6Bplus = Integrate(I6Bplus, TheContourup,abstol=1e-4)[1]
    # ℐ7Aplus = Integrate(I7Aplus, TheContourup,abstol=1e-4)[1]
    # ℐ7Bplus = Integrate(I7Bplus, TheContourup,abstol=1e-4)[1]
    # ℐ8Aplus = Integrate(I8Aplus, TheContourup,abstol=1e-4)[1]
    # ℐ8Bplus = Integrate(I8Bplus, TheContourup,abstol=1e-4)[1]
    # ℐ9plus = Integrate(I9plus, TheContourup,abstol=1e-4)[1]
    ℐ10plus = Integrate(I10plus, TheContourup,abstol=1e-4)[1]

    # Downwards-facing contours

    ∂ω𝒪minus = Integrate(dwOminus, TheContourdown,abstol=1e-4)[1]

    ℋ1Aminus = Integrate(H1Aminus, TheContourdown,abstol=1e-4)[1]
    ℋ1Bminus = Integrate(H1Bminus, TheContourdown,abstol=1e-4)[1]
    ℋ2Aminus = Integrate(H2Aminus, TheContourdown,abstol=1e-4)[1]
    ℋ2Bminus = Integrate(H2Bminus, TheContourdown,abstol=1e-4)[1]
    ℋ3Aminus = Integrate(H3Aminus, TheContourdown,abstol=1e-4)[1]
    ℋ3Bminus = Integrate(H3Bminus, TheContourdown,abstol=1e-4)[1]
    ℋ4Aminus = Integrate(H4Aminus, TheContourdown,abstol=1e-4)[1]
    ℋ4Bminus = Integrate(H4Bminus, TheContourdown,abstol=1e-4)[1]
    # ℋ5Aminus = Integrate(H5Aminus, TheContourdown,abstol=1e-4)[1]
    # ℋ5Bminus = Integrate(H5Bminus, TheContourdown,abstol=1e-4)[1]
    # ℋ6Aminus = Integrate(H6Aminus, TheContourdown,abstol=1e-4)[1]
    # ℋ6Bminus = Integrate(H6Bminus, TheContourdown,abstol=1e-4)[1]
    # ℋ7Aminus = Integrate(H7Aminus, TheContourdown,abstol=1e-4)[1]
    # ℋ7Bminus = Integrate(H7Bminus, TheContourdown,abstol=1e-4)[1]
    # ℋ8Aminus = Integrate(H8Aminus, TheContourdown,abstol=1e-4)[1]
    # ℋ8Bminus = Integrate(H8Bminus, TheContourdown,abstol=1e-4)[1]
    # ℋ9minus = Integrate(H9minus, TheContourdown,abstol=1e-4)[1]
    ℋ10minus = Integrate(H10minus, TheContourdown,abstol=1e-4)[1]

    ℐ1Aminus = Integrate(I1Aminus, TheContourdown,abstol=1e-4)[1]
    ℐ1Bminus = Integrate(I1Bminus, TheContourdown,abstol=1e-4)[1]
    ℐ2Aminus = Integrate(I2Aminus, TheContourdown,abstol=1e-4)[1]
    ℐ2Bminus = Integrate(I2Bminus, TheContourdown,abstol=1e-4)[1]
    ℐ3Aminus = Integrate(I3Aminus, TheContourdown,abstol=1e-4)[1]
    ℐ3Bminus = Integrate(I3Bminus, TheContourdown,abstol=1e-4)[1]
    ℐ4Aminus = Integrate(I4Aminus, TheContourdown,abstol=1e-4)[1]
    ℐ4Bminus = Integrate(I4Bminus, TheContourdown,abstol=1e-4)[1]
    # ℐ5Aminus = Integrate(I5Aminus, TheContourdown,abstol=1e-4)[1]
    # ℐ5Bminus = Integrate(I5Bminus, TheContourdown,abstol=1e-4)[1]
    # ℐ6Aminus = Integrate(I6Aminus, TheContourdown,abstol=1e-4)[1]
    # ℐ6Bminus = Integrate(I6Bminus, TheContourdown,abstol=1e-4)[1]
    # ℐ7Aminus = Integrate(I7Aminus, TheContourdown,abstol=1e-4)[1]
    # ℐ7Bminus = Integrate(I7Bminus, TheContourdown,abstol=1e-4)[1]
    # ℐ8Aminus = Integrate(I8Aminus, TheContourdown,abstol=1e-4)[1]
    # ℐ8Bminus = Integrate(I8Bminus, TheContourdown,abstol=1e-4)[1]
    # ℐ9minus = Integrate(I9minus, TheContourdown,abstol=1e-4)[1]
    ℐ10minus = Integrate(I10minus, TheContourdown,abstol=1e-4)[1]

    # ℋplus = ℋ1Aplus+ℋ1Bplus+ℋ2Aplus+ℋ2Bplus+ℋ3Aplus+ℋ3Bplus+ℋ4Aplus+ℋ4Bplus+ℋ5Aplus+ℋ5Bplus+ℋ6Aplus+ℋ6Bplus+ℋ7Aplus+ℋ7Bplus+ℋ8Aplus+ℋ8Bplus+ℋ9plus+ℋ10plus
    # ℐplus = ℐ1Aplus+ℐ1Bplus+ℐ2Aplus+ℐ2Bplus+ℐ3Aplus+ℐ3Bplus+ℐ4Aplus+ℐ4Bplus+ℐ5Aplus+ℐ5Bplus+ℐ6Aplus+ℐ6Bplus+ℐ7Aplus+ℐ7Bplus+ℐ8Aplus+ℐ8Bplus+ℐ9plus+ℐ10plus

    # ℋminus = ℋ1Aminus+ℋ1Bminus+ℋ2Aminus+ℋ2Bminus+ℋ3Aminus+ℋ3Bminus+ℋ4Aminus+ℋ4Bminus+ℋ5Aminus+ℋ5Bminus+ℋ6Aminus+ℋ6Bminus+ℋ7Aminus+ℋ7Bminus+ℋ8Aminus+ℋ8Bminus+ℋ9minus+ℋ10minus
    # ℐminus = ℐ1Aminus+ℐ1Bminus+ℐ2Aminus+ℐ2Bminus+ℐ3Aminus+ℐ3Bminus+ℐ4Aminus+ℐ4Bminus+ℐ5Aminus+ℐ5Bminus+ℐ6Aminus+ℐ6Bminus+ℐ7Aminus+ℐ7Bminus+ℐ8Aminus+ℐ8Bminus+ℐ9minus+ℐ10minus

    ℋplus = ℋ1Aplus+ℋ1Bplus+ℋ2Aplus+ℋ2Bplus+ℋ3Aplus+ℋ3Bplus+ℋ4Aplus+ℋ4Bplus+ℋ10plus
    ℐplus = ℐ1Aplus+ℐ1Bplus+ℐ2Aplus+ℐ2Bplus+ℐ3Aplus+ℐ3Bplus+ℐ4Aplus+ℐ4Bplus+ℐ10plus

    ℋminus = ℋ1Aminus+ℋ1Bminus+ℋ2Aminus+ℋ2Bminus+ℋ3Aminus+ℋ3Bminus+ℋ4Aminus+ℋ4Bminus+ℋ10minus
    ℐminus = ℐ1Aminus+ℐ1Bminus+ℐ2Aminus+ℐ2Bminus+ℐ3Aminus+ℐ3Bminus+ℐ4Aminus+ℐ4Bminus+ℐ10minus

    ω2s=Computeω2(∂ω𝒪plus,∂ω𝒪minus,ℋplus,ℋminus,ℐplus,ℐminus,ψ)
    @show ω2s
    
end


@testset "PartOfFullJP" begin   

    # dwOfile = "./OperatorShifts/OrderReducedShifts/dwscalarOcoefficients.csv"
    # dϵOfile = "./OperatorShifts/OrderReducedShifts/δscalarOcoefficients.csv"
    # dwOfile = "./OperatorShifts/OrderReducedShifts/NPdwscalarOcoefficients.csv"
    # dϵOfile = "./OperatorShifts/OrderReducedShifts/NPδscalarOcoefficients.csv"
    dwOfile = "./OperatorShifts/FullNPTest/PartialOps/dwOcoefficients.csv"
    ℋfile = "./OperatorShifts/FullNPTest/PartialOps/PartialHcoefficients.csv"
    ℐfile = "./OperatorShifts/FullNPTest/PartialOps/PartialIcoefficients.csv"

    ∂ωO = OperatorShift(dwOfile)
    H = OperatorShift(ℋfile)
    I = OperatorShift(ℐfile)

    println("Made operator shifts")
    
    ψ = qnmfunctionnew(-2,2,2,0,0.1)
    ψmconj = qnmfunctionnew(-2,2,2,0,0.1,is_conjugate=true,is_minus=true)

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
    ∂ωOm2 = OperatorSandwich(ψ,∂ωO,weight,ψ).Op
    Hm2 = OperatorSandwich(ψ,H,weight,ψ).Op
    Im2 = OperatorSandwich(ψ,I,weight,ψmconj).Op
    println("Made Operators")

    ∂ω𝒪m2 = Integrate(∂ωOm2, TheContourup,abstol=1e-4)[1]
    ℋm2 = Integrate(Hm2, TheContourup,abstol=1e-4)[1]
    ℐm2 = Integrate(Im2, TheContourup,abstol=1e-4)[1]

    # @show δω = -(∂ϵ𝒪m2/∂ω𝒪m2)
end
