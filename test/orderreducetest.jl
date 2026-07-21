using ContourIntegrals
using KerrQNMShifts
using KerrQuasinormalModes
using Test
using CSV, DataFrames

println("Done usings")

@testset "JPPertToKerrManySpins" begin 

    dwOfile = "./OperatorShifts/FullNPTest/FullG2Take5/dwOOps/dwOcoefficients.csv"

    H1Afile = "./OperatorShifts/FullNPTest/FullG2Take5/HOps/H1A.csv"
    H1Bfile = "./OperatorShifts/FullNPTest/FullG2Take5/HOps/H1B.csv"
    H2Afile = "./OperatorShifts/FullNPTest/FullG2Take5/HOps/H2A.csv"
    H2Bfile = "./OperatorShifts/FullNPTest/FullG2Take5/HOps/H2B.csv"
    H3Afile = "./OperatorShifts/FullNPTest/FullG2Take5/HOps/H3A.csv"
    H3Bfile = "./OperatorShifts/FullNPTest/FullG2Take5/HOps/H3B.csv"
    H4Afile = "./OperatorShifts/FullNPTest/FullG2Take5/HOps/H4A.csv"
    H4Bfile = "./OperatorShifts/FullNPTest/FullG2Take5/HOps/H4B.csv"
    H10file = "./OperatorShifts/FullNPTest/FullG2Take5/HOps/H10.csv"

    I1Afile = "./OperatorShifts/FullNPTest/FullG2Take5/IOps/I1A.csv"
    I1Bfile = "./OperatorShifts/FullNPTest/FullG2Take5/IOps/I1B.csv"
    I2Afile = "./OperatorShifts/FullNPTest/FullG2Take5/IOps/I2A.csv"
    I2Bfile = "./OperatorShifts/FullNPTest/FullG2Take5/IOps/I2B.csv"
    I3Afile = "./OperatorShifts/FullNPTest/FullG2Take5/IOps/I3A.csv"
    I3Bfile = "./OperatorShifts/FullNPTest/FullG2Take5/IOps/I3B.csv"
    I4Afile = "./OperatorShifts/FullNPTest/FullG2Take5/IOps/I4A.csv"
    I4Bfile = "./OperatorShifts/FullNPTest/FullG2Take5/IOps/I4B.csv"
    I10file = "./OperatorShifts/FullNPTest/FullG2Take5/IOps/I10.csv"

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

    freqpertsmatrix =DataFrame(l = Int[], m = Int[], n = Int[], a = Float64[] , δωE = ComplexF64[], γE = ComplexF64[], δωO = ComplexF64[], γO = ComplexF64[])

    l=3
    m=1
    n=0

    filename = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/FullNP/QNMPlots/JPShifts" * string(l, m, n) * ".csv"
    
    for a in collect(0.0+1e-8:0.01:0.8)

        println(l,", ",m,", ",n,", ",a)

        ψ = qnmfunctionnew(-2,l,m,n,a)
        ψm = qnmfunctionnew(-2,l,m,n,a, is_minus=true)
        ψconj =  qnmfunctionnew(-2,l,m,n,a, is_conjugate=true)
        ψmconj = qnmfunctionnew(-2,l,m,n,a,is_conjugate=true,is_minus=true)

        # Compile ψ
        ψ(1,.5)
        # println("Past ψ compile")

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
        I10plus = OperatorSandwich(ψ,I10,weight,ψmconj).Op

        I1Aminus = OperatorSandwich(ψm,I1A,weight,ψconj).Op
        I1Bminus = OperatorSandwich(ψm,I1B,weight,ψconj).Op
        I2Aminus = OperatorSandwich(ψm,I2A,weight,ψconj).Op
        I2Bminus = OperatorSandwich(ψm,I2B,weight,ψconj).Op
        I3Aminus = OperatorSandwich(ψm,I3A,weight,ψconj).Op
        I3Bminus = OperatorSandwich(ψm,I3B,weight,ψconj).Op
        I4Aminus = OperatorSandwich(ψm,I4A,weight,ψconj).Op
        I4Bminus = OperatorSandwich(ψm,I4B,weight,ψconj).Op
        I10minus = OperatorSandwich(ψm,I10,weight,ψconj).Op

        println("Made Operators")

        # Upwards-facing contours

        ℋ1Aplus = Integrate(H1Aplus, TheContourup,abstol=1e-5)[1]
        ℋ1Bplus = Integrate(H1Bplus, TheContourup,abstol=1e-5)[1]
        ℋ2Aplus = Integrate(H2Aplus, TheContourup,abstol=1e-5)[1]
        ℋ2Bplus = Integrate(H2Bplus, TheContourup,abstol=1e-5)[1]
        ℋ3Aplus = Integrate(H3Aplus, TheContourup,abstol=1e-5)[1]
        ℋ3Bplus = Integrate(H3Bplus, TheContourup,abstol=1e-5)[1]
        ℋ4Aplus = Integrate(H4Aplus, TheContourup,abstol=1e-5)[1]
        ℋ4Bplus = Integrate(H4Bplus, TheContourup,abstol=1e-5)[1]
        ℋ10plus = Integrate(H10plus, TheContourup,abstol=1e-5)[1]
        println("Done ℋplus")

        ℐ1Aplus = Integrate(I1Aplus, TheContourup,abstol=1e-5)[1]
        ℐ1Bplus = Integrate(I1Bplus, TheContourup,abstol=1e-5)[1]
        ℐ2Aplus = Integrate(I2Aplus, TheContourup,abstol=1e-5)[1]
        ℐ2Bplus = Integrate(I2Bplus, TheContourup,abstol=1e-5)[1]
        ℐ3Aplus = Integrate(I3Aplus, TheContourup,abstol=1e-5)[1]
        ℐ3Bplus = Integrate(I3Bplus, TheContourup,abstol=1e-5)[1]
        ℐ4Aplus = Integrate(I4Aplus, TheContourup,abstol=1e-5)[1]
        ℐ4Bplus = Integrate(I4Bplus, TheContourup,abstol=1e-5)[1]
        ℐ10plus = Integrate(I10plus, TheContourup,abstol=1e-5)[1]
        println("Done ℐplus")

        ∂ω𝒪plus = Integrate(dwOplus, TheContourup,abstol=1e-5)[1]
        ℋplus = ℋ1Aplus+ℋ1Bplus+ℋ2Aplus+ℋ2Bplus+ℋ3Aplus+ℋ3Bplus+ℋ4Aplus+ℋ4Bplus+ℋ10plus
        ℐplus = ℐ1Aplus+ℐ1Bplus+ℐ2Aplus+ℐ2Bplus+ℐ3Aplus+ℐ3Bplus+ℐ4Aplus+ℐ4Bplus+ℐ10plus

        # Downwards-facing contours

        ℋ1Aminus = conj(Integrate(H1Aminus, TheContourdown,abstol=1e-5)[1])
        ℋ1Bminus = conj(Integrate(H1Bminus, TheContourdown,abstol=1e-5)[1])
        ℋ2Aminus = conj(Integrate(H2Aminus, TheContourdown,abstol=1e-5)[1])
        ℋ2Bminus = conj(Integrate(H2Bminus, TheContourdown,abstol=1e-5)[1])
        ℋ3Aminus = conj(Integrate(H3Aminus, TheContourdown,abstol=1e-5)[1])
        ℋ3Bminus = conj(Integrate(H3Bminus, TheContourdown,abstol=1e-5)[1])
        ℋ4Aminus = conj(Integrate(H4Aminus, TheContourdown,abstol=1e-5)[1])
        ℋ4Bminus = conj(Integrate(H4Bminus, TheContourdown,abstol=1e-5)[1])
        ℋ10minus = conj(Integrate(H10minus, TheContourdown,abstol=1e-5)[1])
        println("Done ℋminus")

        ℐ1Aminus = conj(Integrate(I1Aminus, TheContourdown,abstol=1e-5)[1])
        ℐ1Bminus = conj(Integrate(I1Bminus, TheContourdown,abstol=1e-5)[1])
        ℐ2Aminus = conj(Integrate(I2Aminus, TheContourdown,abstol=1e-5)[1])
        ℐ2Bminus = conj(Integrate(I2Bminus, TheContourdown,abstol=1e-5)[1])
        ℐ3Aminus = conj(Integrate(I3Aminus, TheContourdown,abstol=1e-5)[1])
        ℐ3Bminus = conj(Integrate(I3Bminus, TheContourdown,abstol=1e-5)[1])
        ℐ4Aminus = conj(Integrate(I4Aminus, TheContourdown,abstol=1e-5)[1])
        ℐ4Bminus = conj(Integrate(I4Bminus, TheContourdown,abstol=1e-5)[1])
        ℐ10minus = conj(Integrate(I10minus, TheContourdown,abstol=1e-5)[1])
        println("Done ℐminus")

        ∂ω𝒪minus = conj(Integrate(dwOminus, TheContourdown,abstol=1e-5)[1])
        ℋminus = ℋ1Aminus+ℋ1Bminus+ℋ2Aminus+ℋ2Bminus+ℋ3Aminus+ℋ3Bminus+ℋ4Aminus+ℋ4Bminus+ℋ10minus
        ℐminus = ℐ1Aminus+ℐ1Bminus+ℐ2Aminus+ℐ2Bminus+ℐ3Aminus+ℐ3Bminus+ℐ4Aminus+ℐ4Bminus+ℐ10minus

        δωs=Computeω2(∂ω𝒪plus,∂ω𝒪minus,ℋplus,ℋminus,ℐplus,ℐminus,ψ)
        @show δωs

        push!(freqpertsmatrix, (l, m, n, a, δωs[1], δωs[3], δωs[2], δωs[4]))
        
        CSV.write(filename, freqpertsmatrix)
    end

end

@testset "JPPertToKerr" begin 

    dwOfile = "./OperatorShifts/FullNPTest/FullG2Take5/dwOOps/dwOcoefficients.csv"

    H1Afile = "./OperatorShifts/FullNPTest/FullG2Take5/HOps/H1A.csv"
    H1Bfile = "./OperatorShifts/FullNPTest/FullG2Take5/HOps/H1B.csv"
    H2Afile = "./OperatorShifts/FullNPTest/FullG2Take5/HOps/H2A.csv"
    H2Bfile = "./OperatorShifts/FullNPTest/FullG2Take5/HOps/H2B.csv"
    H3Afile = "./OperatorShifts/FullNPTest/FullG2Take5/HOps/H3A.csv"
    H3Bfile = "./OperatorShifts/FullNPTest/FullG2Take5/HOps/H3B.csv"
    H4Afile = "./OperatorShifts/FullNPTest/FullG2Take5/HOps/H4A.csv"
    H4Bfile = "./OperatorShifts/FullNPTest/FullG2Take5/HOps/H4B.csv"
    H10file = "./OperatorShifts/FullNPTest/FullG2Take5/HOps/H10.csv"

    I1Afile = "./OperatorShifts/FullNPTest/FullG2Take5/IOps/I1A.csv"
    I1Bfile = "./OperatorShifts/FullNPTest/FullG2Take5/IOps/I1B.csv"
    I2Afile = "./OperatorShifts/FullNPTest/FullG2Take5/IOps/I2A.csv"
    I2Bfile = "./OperatorShifts/FullNPTest/FullG2Take5/IOps/I2B.csv"
    I3Afile = "./OperatorShifts/FullNPTest/FullG2Take5/IOps/I3A.csv"
    I3Bfile = "./OperatorShifts/FullNPTest/FullG2Take5/IOps/I3B.csv"
    I4Afile = "./OperatorShifts/FullNPTest/FullG2Take5/IOps/I4A.csv"
    I4Bfile = "./OperatorShifts/FullNPTest/FullG2Take5/IOps/I4B.csv"
    I10file = "./OperatorShifts/FullNPTest/FullG2Take5/IOps/I10.csv"

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

    freqpertsmatrix =DataFrame(l = Int[], m = Int[], n = Int[], a = Float64[] , δωE = ComplexF64[], γE = ComplexF64[], δωO = ComplexF64[], γO = ComplexF64[])

    l=2
    m=2
    n=0
    a=0.01

    filename = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/FullNP/QNMPlots/LinearizedKerr/KerrShifts" * string(l, m, n) * ".csv"
    
    println(l,", ",m,", ",n,", ",a)

    ψ = qnmfunctionnew(-2,l,m,n,a)
    ψm = qnmfunctionnew(-2,l,m,n,a, is_minus=true)
    ψconj =  qnmfunctionnew(-2,l,m,n,a, is_conjugate=true)
    ψmconj = qnmfunctionnew(-2,l,m,n,a,is_conjugate=true,is_minus=true)

    # Compile ψ
    ψ(1,.5)
    # println("Past ψ compile")

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
    I10plus = OperatorSandwich(ψ,I10,weight,ψmconj).Op

    I1Aminus = OperatorSandwich(ψm,I1A,weight,ψconj).Op
    I1Bminus = OperatorSandwich(ψm,I1B,weight,ψconj).Op
    I2Aminus = OperatorSandwich(ψm,I2A,weight,ψconj).Op
    I2Bminus = OperatorSandwich(ψm,I2B,weight,ψconj).Op
    I3Aminus = OperatorSandwich(ψm,I3A,weight,ψconj).Op
    I3Bminus = OperatorSandwich(ψm,I3B,weight,ψconj).Op
    I4Aminus = OperatorSandwich(ψm,I4A,weight,ψconj).Op
    I4Bminus = OperatorSandwich(ψm,I4B,weight,ψconj).Op
    I10minus = OperatorSandwich(ψm,I10,weight,ψconj).Op

    println("Made Operators")

        # Upwards-facing contours

    ℋ1Aplus = Integrate(H1Aplus, TheContourup,abstol=1e-5)[1]
    ℋ1Bplus = Integrate(H1Bplus, TheContourup,abstol=1e-5)[1]
    ℋ2Aplus = Integrate(H2Aplus, TheContourup,abstol=1e-5)[1]
    ℋ2Bplus = Integrate(H2Bplus, TheContourup,abstol=1e-5)[1]
    ℋ3Aplus = Integrate(H3Aplus, TheContourup,abstol=1e-5)[1]
    ℋ3Bplus = Integrate(H3Bplus, TheContourup,abstol=1e-5)[1]
    ℋ4Aplus = Integrate(H4Aplus, TheContourup,abstol=1e-5)[1]
    ℋ4Bplus = Integrate(H4Bplus, TheContourup,abstol=1e-5)[1]
    ℋ10plus = Integrate(H10plus, TheContourup,abstol=1e-5)[1]
    println("Done ℋplus")

    ℐ1Aplus = Integrate(I1Aplus, TheContourup,abstol=1e-5)[1]
    ℐ1Bplus = Integrate(I1Bplus, TheContourup,abstol=1e-5)[1]
    ℐ2Aplus = Integrate(I2Aplus, TheContourup,abstol=1e-5)[1]
    ℐ2Bplus = Integrate(I2Bplus, TheContourup,abstol=1e-5)[1]
    ℐ3Aplus = Integrate(I3Aplus, TheContourup,abstol=1e-5)[1]
    ℐ3Bplus = Integrate(I3Bplus, TheContourup,abstol=1e-5)[1]
    ℐ4Aplus = Integrate(I4Aplus, TheContourup,abstol=1e-5)[1]
    ℐ4Bplus = Integrate(I4Bplus, TheContourup,abstol=1e-5)[1]
    ℐ10plus = Integrate(I10plus, TheContourup,abstol=1e-5)[1]
    println("Done ℐplus")

    ∂ω𝒪plus = Integrate(dwOplus, TheContourup,abstol=1e-5)[1]
    ℋplus = ℋ1Aplus+ℋ1Bplus+ℋ2Aplus+ℋ2Bplus+ℋ3Aplus+ℋ3Bplus+ℋ4Aplus+ℋ4Bplus+ℋ10plus
    ℐplus = ℐ1Aplus+ℐ1Bplus+ℐ2Aplus+ℐ2Bplus+ℐ3Aplus+ℐ3Bplus+ℐ4Aplus+ℐ4Bplus+ℐ10plus

        # Downwards-facing contours

    ℋ1Aminus = conj(Integrate(H1Aminus, TheContourdown,abstol=1e-5)[1])
    ℋ1Bminus = conj(Integrate(H1Bminus, TheContourdown,abstol=1e-5)[1])
    ℋ2Aminus = conj(Integrate(H2Aminus, TheContourdown,abstol=1e-5)[1])
    ℋ2Bminus = conj(Integrate(H2Bminus, TheContourdown,abstol=1e-5)[1])
    ℋ3Aminus = conj(Integrate(H3Aminus, TheContourdown,abstol=1e-5)[1])
    ℋ3Bminus = conj(Integrate(H3Bminus, TheContourdown,abstol=1e-5)[1])
    ℋ4Aminus = conj(Integrate(H4Aminus, TheContourdown,abstol=1e-5)[1])
    ℋ4Bminus = conj(Integrate(H4Bminus, TheContourdown,abstol=1e-5)[1])
    ℋ10minus = conj(Integrate(H10minus, TheContourdown,abstol=1e-5)[1])
    println("Done ℋminus")

    ℐ1Aminus = conj(Integrate(I1Aminus, TheContourdown,abstol=1e-5)[1])
    ℐ1Bminus = conj(Integrate(I1Bminus, TheContourdown,abstol=1e-5)[1])
    ℐ2Aminus = conj(Integrate(I2Aminus, TheContourdown,abstol=1e-5)[1])
    ℐ2Bminus = conj(Integrate(I2Bminus, TheContourdown,abstol=1e-5)[1])
    ℐ3Aminus = conj(Integrate(I3Aminus, TheContourdown,abstol=1e-5)[1])
    ℐ3Bminus = conj(Integrate(I3Bminus, TheContourdown,abstol=1e-5)[1])
    ℐ4Aminus = conj(Integrate(I4Aminus, TheContourdown,abstol=1e-5)[1])
    ℐ4Bminus = conj(Integrate(I4Bminus, TheContourdown,abstol=1e-5)[1])
    ℐ10minus = conj(Integrate(I10minus, TheContourdown,abstol=1e-5)[1])
    println("Done ℐminus")

    ∂ω𝒪minus = conj(Integrate(dwOminus, TheContourdown,abstol=1e-5)[1])
    ℋminus = ℋ1Aminus+ℋ1Bminus+ℋ2Aminus+ℋ2Bminus+ℋ3Aminus+ℋ3Bminus+ℋ4Aminus+ℋ4Bminus+ℋ10minus
    ℐminus = ℐ1Aminus+ℐ1Bminus+ℐ2Aminus+ℐ2Bminus+ℐ3Aminus+ℐ3Bminus+ℐ4Aminus+ℐ4Bminus+ℐ10minus

    δωs=Computeω2(∂ω𝒪plus,∂ω𝒪minus,ℋplus,ℋminus,ℐplus,ℐminus,ψ)
    @show δωs

    push!(freqpertsmatrix, (l, m, n, a, δωs[1], δωs[3], δωs[2], δωs[4]))
        
    CSV.write(filename, freqpertsmatrix)

end

@testset "KerrPertToSchw" begin 

    dwOfile = "./OperatorShifts/FullNPTest/KerrPertToSchw/dwOOps/dwOcoefficients.csv"

    H1Afile = "./OperatorShifts/FullNPTest/KerrPertToSchw/HOps/H1A.csv"
    H1Bfile = "./OperatorShifts/FullNPTest/KerrPertToSchw/HOps/H1B.csv"
    H2Afile = "./OperatorShifts/FullNPTest/KerrPertToSchw/HOps/H2A.csv"
    H2Bfile = "./OperatorShifts/FullNPTest/KerrPertToSchw/HOps/H2B.csv"
    H3Afile = "./OperatorShifts/FullNPTest/KerrPertToSchw/HOps/H3A.csv"
    H3Bfile = "./OperatorShifts/FullNPTest/KerrPertToSchw/HOps/H3B.csv"
    H4Afile = "./OperatorShifts/FullNPTest/KerrPertToSchw/HOps/H4A.csv"
    H4Bfile = "./OperatorShifts/FullNPTest/KerrPertToSchw/HOps/H4B.csv"
    H10file = "./OperatorShifts/FullNPTest/KerrPertToSchw/HOps/H10.csv"

    I1Afile = "./OperatorShifts/FullNPTest/KerrPertToSchw/IOps/I1A.csv"
    I1Bfile = "./OperatorShifts/FullNPTest/KerrPertToSchw/IOps/I1B.csv"
    I2Afile = "./OperatorShifts/FullNPTest/KerrPertToSchw/IOps/I2A.csv"
    I2Bfile = "./OperatorShifts/FullNPTest/KerrPertToSchw/IOps/I2B.csv"
    I3Afile = "./OperatorShifts/FullNPTest/KerrPertToSchw/IOps/I3A.csv"
    I3Bfile = "./OperatorShifts/FullNPTest/KerrPertToSchw/IOps/I3B.csv"
    I4Afile = "./OperatorShifts/FullNPTest/KerrPertToSchw/IOps/I4A.csv"
    I4Bfile = "./OperatorShifts/FullNPTest/KerrPertToSchw/IOps/I4B.csv"
    I10file = "./OperatorShifts/FullNPTest/KerrPertToSchw/IOps/I10.csv"

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

    l=3
    m=1
    n=0
    a=10^(-8)
    
    println(l,", ",m,", ",n,", ",a)

    ψ = qnmfunctionnew(-2,l,m,n,a)
    ψm = qnmfunctionnew(-2,l,m,n,a, is_minus=true)
    ψconj =  qnmfunctionnew(-2,l,m,n,a, is_conjugate=true)
    ψmconj = qnmfunctionnew(-2,l,m,n,a,is_conjugate=true,is_minus=true)

    # Compile ψ
    ψ(1,.5)
    # println("Past ψ compile")

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
    I10plus = OperatorSandwich(ψ,I10,weight,ψmconj).Op

    I1Aminus = OperatorSandwich(ψm,I1A,weight,ψconj).Op
    I1Bminus = OperatorSandwich(ψm,I1B,weight,ψconj).Op
    I2Aminus = OperatorSandwich(ψm,I2A,weight,ψconj).Op
    I2Bminus = OperatorSandwich(ψm,I2B,weight,ψconj).Op
    I3Aminus = OperatorSandwich(ψm,I3A,weight,ψconj).Op
    I3Bminus = OperatorSandwich(ψm,I3B,weight,ψconj).Op
    I4Aminus = OperatorSandwich(ψm,I4A,weight,ψconj).Op
    I4Bminus = OperatorSandwich(ψm,I4B,weight,ψconj).Op
    I10minus = OperatorSandwich(ψm,I10,weight,ψconj).Op

    println("Made Operators")

        # Upwards-facing contours

    ℋ1Aplus = Integrate(H1Aplus, TheContourup,abstol=1e-5)[1]
    ℋ1Bplus = Integrate(H1Bplus, TheContourup,abstol=1e-5)[1]
    ℋ2Aplus = Integrate(H2Aplus, TheContourup,abstol=1e-5)[1]
    ℋ2Bplus = Integrate(H2Bplus, TheContourup,abstol=1e-5)[1]
    ℋ3Aplus = Integrate(H3Aplus, TheContourup,abstol=1e-5)[1]
    ℋ3Bplus = Integrate(H3Bplus, TheContourup,abstol=1e-5)[1]
    ℋ4Aplus = Integrate(H4Aplus, TheContourup,abstol=1e-5)[1]
    ℋ4Bplus = Integrate(H4Bplus, TheContourup,abstol=1e-5)[1]
    ℋ10plus = Integrate(H10plus, TheContourup,abstol=1e-5)[1]
    println("Done ℋplus")

    ℐ1Aplus = Integrate(I1Aplus, TheContourup,abstol=1e-5)[1]
    ℐ1Bplus = Integrate(I1Bplus, TheContourup,abstol=1e-5)[1]
    ℐ2Aplus = Integrate(I2Aplus, TheContourup,abstol=1e-5)[1]
    ℐ2Bplus = Integrate(I2Bplus, TheContourup,abstol=1e-5)[1]
    ℐ3Aplus = Integrate(I3Aplus, TheContourup,abstol=1e-5)[1]
    ℐ3Bplus = Integrate(I3Bplus, TheContourup,abstol=1e-5)[1]
    ℐ4Aplus = Integrate(I4Aplus, TheContourup,abstol=1e-5)[1]
    ℐ4Bplus = Integrate(I4Bplus, TheContourup,abstol=1e-5)[1]
    ℐ10plus = Integrate(I10plus, TheContourup,abstol=1e-5)[1]
    println("Done ℐplus")

    ∂ω𝒪plus = Integrate(dwOplus, TheContourup,abstol=1e-5)[1]
    ℋplus = ℋ1Aplus+ℋ1Bplus+ℋ2Aplus+ℋ2Bplus+ℋ3Aplus+ℋ3Bplus+ℋ4Aplus+ℋ4Bplus+ℋ10plus
    ℐplus = ℐ1Aplus+ℐ1Bplus+ℐ2Aplus+ℐ2Bplus+ℐ3Aplus+ℐ3Bplus+ℐ4Aplus+ℐ4Bplus+ℐ10plus

        # Downwards-facing contours

    ℋ1Aminus = conj(Integrate(H1Aminus, TheContourdown,abstol=1e-5)[1])
    ℋ1Bminus = conj(Integrate(H1Bminus, TheContourdown,abstol=1e-5)[1])
    ℋ2Aminus = conj(Integrate(H2Aminus, TheContourdown,abstol=1e-5)[1])
    ℋ2Bminus = conj(Integrate(H2Bminus, TheContourdown,abstol=1e-5)[1])
    ℋ3Aminus = conj(Integrate(H3Aminus, TheContourdown,abstol=1e-5)[1])
    ℋ3Bminus = conj(Integrate(H3Bminus, TheContourdown,abstol=1e-5)[1])
    ℋ4Aminus = conj(Integrate(H4Aminus, TheContourdown,abstol=1e-5)[1])
    ℋ4Bminus = conj(Integrate(H4Bminus, TheContourdown,abstol=1e-5)[1])
    ℋ10minus = conj(Integrate(H10minus, TheContourdown,abstol=1e-5)[1])
    println("Done ℋminus")

    ℐ1Aminus = conj(Integrate(I1Aminus, TheContourdown,abstol=1e-5)[1])
    ℐ1Bminus = conj(Integrate(I1Bminus, TheContourdown,abstol=1e-5)[1])
    ℐ2Aminus = conj(Integrate(I2Aminus, TheContourdown,abstol=1e-5)[1])
    ℐ2Bminus = conj(Integrate(I2Bminus, TheContourdown,abstol=1e-5)[1])
    ℐ3Aminus = conj(Integrate(I3Aminus, TheContourdown,abstol=1e-5)[1])
    ℐ3Bminus = conj(Integrate(I3Bminus, TheContourdown,abstol=1e-5)[1])
    ℐ4Aminus = conj(Integrate(I4Aminus, TheContourdown,abstol=1e-5)[1])
    ℐ4Bminus = conj(Integrate(I4Bminus, TheContourdown,abstol=1e-5)[1])
    ℐ10minus = conj(Integrate(I10minus, TheContourdown,abstol=1e-5)[1])
    println("Done ℐminus")

    ∂ω𝒪minus = conj(Integrate(dwOminus, TheContourdown,abstol=1e-5)[1])
    ℋminus = ℋ1Aminus+ℋ1Bminus+ℋ2Aminus+ℋ2Bminus+ℋ3Aminus+ℋ3Bminus+ℋ4Aminus+ℋ4Bminus+ℋ10minus
    ℐminus = ℐ1Aminus+ℐ1Bminus+ℐ2Aminus+ℐ2Bminus+ℐ3Aminus+ℐ3Bminus+ℐ4Aminus+ℐ4Bminus+ℐ10minus

    δωs=Computeω2(∂ω𝒪plus,∂ω𝒪minus,ℋplus,ℋminus,ℐplus,ℐminus,ψ)
    @show δωs

end

@testset "EasyKerrPertToSchw" begin 

    dwOfile = "./OperatorShifts/EasyKerrPertToSchw/dwOcoefficients.csv"

    δOfile = "./OperatorShifts/EasyKerrPertToSchw/δOcoefficients.csv"

    dwO = OperatorShift(dwOfile)
    δO = OperatorShift(δOfile)


    println("Made operator shifts")

    l=4
    m=1
    n=0
    a=10^(-8)
    s=-1

    
    println(l,", ",m,", ",n,", ",a)

    ψ = qnmfunctionnew(s,l,m,n,a)
    ψm = qnmfunctionnew(s,l,m,n,a, is_minus=true)
    ψconj =  qnmfunctionnew(s,l,m,n,a, is_conjugate=true)
    ψmconj = qnmfunctionnew(s,l,m,n,a,is_conjugate=true,is_minus=true)

    ψ.S(.3)
    ψm.S(-.3)

    # Compile ψ
    ψ(1,.5)
    # println("Past ψ compile")

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

    δOplus = OperatorSandwich(ψ,δO,weight,ψ).Op
    δOminus = OperatorSandwich(ψm,δO,weight,ψm).Op

    println("Made Operators")

    # Upwards-facing contours

    δ𝒪plus = Integrate(δOplus, TheContourup,abstol=1e-5)[1]

    ∂ω𝒪plus = Integrate(dwOplus, TheContourup,abstol=1e-5)[1]

    # Downwards-facing contours

    δ𝒪minus = conj(Integrate(δOminus, TheContourdown,abstol=1e-5)[1])

    ∂ω𝒪minus = conj(Integrate(dwOminus, TheContourdown,abstol=1e-5)[1])

    @show -δ𝒪plus/∂ω𝒪plus
    @show δ𝒪minus/∂ω𝒪minus

end

@testset "EasyConjKerrPertToSchw" begin 

    dwOfile = "./OperatorShifts/EasyKerrPertToSchw/ConjdwOcoefficients.csv"

    δOfile = "./OperatorShifts/EasyKerrPertToSchw/ConjδOcoefficients.csv"

    dwO = OperatorShift(dwOfile)
    δO = OperatorShift(δOfile)


    println("Made operator shifts")

    l=3
    m=1
    n=0
    a=10^(-8)

    
    println(l,", ",m,", ",n,", ",a)

    ψ = qnmfunctionnew(-2,l,m,n,a)
    ψm = qnmfunctionnew(-2,l,m,n,a, is_minus=true)
    ψconj =  qnmfunctionnew(-2,l,m,n,a, is_conjugate=true)
    ψmconj = qnmfunctionnew(-2,l,m,n,a,is_conjugate=true,is_minus=true)

    # Compile ψ
    ψ(1,.5)
    # println("Past ψ compile")

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
    dwOplus = OperatorSandwich(ψconj,dwO,weight,ψconj).Op
    dwOminus = OperatorSandwich(ψmconj,dwO,weight,ψmconj).Op

    δOplus = OperatorSandwich(ψconj,δO,weight,ψconj).Op
    δOminus = OperatorSandwich(ψmconj,δO,weight,ψmconj).Op

    println("Made Operators")

    # Upwards-facing contours

    δ𝒪plus = conj(Integrate(δOplus, TheContourdown,abstol=1e-5)[1])

    ∂ω𝒪plus = conj(Integrate(dwOplus, TheContourdown,abstol=1e-5)[1])

    # Downwards-facing contours

    δ𝒪minus = Integrate(δOminus, TheContourup,abstol=1e-5)[1]

    ∂ω𝒪minus = Integrate(dwOminus, TheContourup,abstol=1e-5)[1]

    @show -δ𝒪plus/∂ω𝒪plus
    @show δ𝒪minus/∂ω𝒪minus

end

@testset "EasyScalarKerrPertToSchw" begin 

    dwOfile = "./OperatorShifts/EasyKerrPertToSchw/ScalardwOcoefficients.csv"

    δOfile = "./OperatorShifts/EasyKerrPertToSchw/ScalarδOcoefficients.csv"

    dwO = OperatorShift(dwOfile)
    δO = OperatorShift(δOfile)


    println("Made operator shifts")

    l=2
    m=2
    n=0
    a=10^(-8)

    
    println(l,", ",m,", ",n,", ",a)

    ψ = qnmfunctionnew(0,l,m,n,a)
    ψm = qnmfunctionnew(0,l,m,n,a, is_minus=true)
    ψconj =  qnmfunctionnew(0,l,m,n,a, is_conjugate=true)
    ψmconj = qnmfunctionnew(0,l,m,n,a,is_conjugate=true,is_minus=true)

    # Compile ψ
    ψ(1,.5)
    # println("Past ψ compile")

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

    δOplus = OperatorSandwich(ψ,δO,weight,ψ).Op
    δOminus = OperatorSandwich(ψm,δO,weight,ψm).Op

    println("Made Operators")

    # Upwards-facing contours

    δ𝒪plus = Integrate(δOplus, TheContourup,abstol=1e-5)[1]

    ∂ω𝒪plus = Integrate(dwOplus, TheContourup,abstol=1e-5)[1]

    # Downwards-facing contours

    δ𝒪minus = conj(Integrate(δOminus, TheContourdown,abstol=1e-5)[1])

    ∂ω𝒪minus = conj(Integrate(dwOminus, TheContourdown,abstol=1e-5)[1])

    @show -δ𝒪plus/∂ω𝒪plus
    @show δ𝒪minus/∂ω𝒪minus

end
