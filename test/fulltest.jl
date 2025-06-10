using ContourIntegrals
using KerrQNMShifts
using KerrQuasinormalModes
using Test
using CSV, DataFrames

println("Done usings")

@testset "ScalarKerrPertToSchwMany" begin   #CHECKED AND GOOD

    dwOfile = "./OperatorShifts/KerrPertToSchw/dwscalarOcoefficients.csv"
    dϵOfile = "./OperatorShifts/KerrPertToSchw/δscalarOcoefficients.csv"
    ∂ϵO = OperatorShift(dϵOfile)
    ∂ωO = OperatorShift(dwOfile)

    println("Made operator shifts")

    freqpertsmatrix =DataFrame(l = Int[], m = Int[], n = Int[], a = Float64[] , δω = ComplexF64[])
    lmax=8

    for l in 2:lmax
        m=l
        for n in 0:2
            for a in [0.0]
                ψ = qnmfunctionnew(0,l,m,n,a)
                # Compile ψ
                ψ(1,.5)
                println("Past ψ compile")
                ## Define the useful contours
                r₊ = ψ.R.r₊ ; r₋ = ψ.R.r₋ ; s = ψ.s ; Δr = 0.1*(r₊-r₋); ϵ = eps(0.1);
                # Define the Weight
                Σ = let a= ψ.a
                    (r,z) -> Complex(r)^2+a^2*z^2
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
                ∂ϵ𝒪m2 = Integrate(∂ϵOm2, TheContourup,abstol=1e-6)[1]
                ∂ω𝒪m2 = Integrate(∂ωOm2, TheContourup,abstol=1e-6)[1]
                @show δω = -(∂ϵ𝒪m2/∂ω𝒪m2)
                push!(freqpertsmatrix, (l, m, n, a, δω))
            end
        end
    end
    CSV.write("/home/dgw763/Documents/LinearizedSpin/ScalarPerturbations/freqpertsmatrix.csv", freqpertsmatrix)
end

@testset "ScalarKerrPertToSchw" begin   #CHECKED AND GOOD

    dwOfile = "./OperatorShifts/KerrPertToSchw/dwscalarOcoefficients.csv"
    dϵOfile = "./OperatorShifts/KerrPertToSchw/δscalarOcoefficients.csv"

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

    ∂ϵ𝒪m2 = Integrate(∂ϵOm2, TheContourup,abstol=1e-6)[1]
    ∂ω𝒪m2 = Integrate(∂ωOm2, TheContourup,abstol=1e-6)[1]

    @show δω = -(∂ϵ𝒪m2/∂ω𝒪m2)
end

@testset "ScalarLinearizedJP" begin   #CHECKED AND GOOD

    dwOfile = "./OperatorShifts/ScalarLinearizedJP/dwscalarOcoefficients.csv"
    dϵOfile = "./OperatorShifts/ScalarLinearizedJP/δscalarOcoefficients.csv"

    ∂ϵO = OperatorShift(dϵOfile)
    ∂ωO = OperatorShift(dwOfile)

    println("Made operator shifts")
    
    ψ = qnmfunctionnew(0,2,2,0,0.1)

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

    ∂ϵ𝒪m2 = Integrate(∂ϵOm2, TheContourup,abstol=1e-6)[1]
    ∂ω𝒪m2 = Integrate(∂ωOm2, TheContourup,abstol=1e-6)[1]

    @show δω = -(∂ϵ𝒪m2/∂ω𝒪m2)
end

@testset "ScalarLinearizedJPMany" begin   #CHECKED AND GOOD

    dwOfile = "./OperatorShifts/ScalarLinearizedJP/dwscalarOcoefficients.csv"
    dϵOfile = "./OperatorShifts/ScalarLinearizedJP/δscalarOcoefficients.csv"

    ∂ϵO = OperatorShift(dϵOfile)
    ∂ωO = OperatorShift(dwOfile)

    println("Made operator shifts")

    freqpertsmatrix =DataFrame(l = Int[], m = Int[], n = Int[], a = Float64[] , δω = ComplexF64[])
    lmax=8

    for l in 2:lmax
        m=l
        for n in 0:2
            for a in [0.0]
                ψ = qnmfunctionnew(0,l,m,n,a)
                # Compile ψ
                ψ(1,.5)
                println("Past ψ compile")
                ## Define the useful contours
                r₊ = ψ.R.r₊ ; r₋ = ψ.R.r₋ ; s = ψ.s ; Δr = 0.1*(r₊-r₋); ϵ = eps(0.1);
                # Define the Weight
                Σ = let a= ψ.a
                    (r,z) -> Complex(r)^2+a^2*z^2
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
                ∂ϵ𝒪m2 = Integrate(∂ϵOm2, TheContourup,abstol=1e-6)[1]
                ∂ω𝒪m2 = Integrate(∂ωOm2, TheContourup,abstol=1e-6)[1]
                @show δω = -(∂ϵ𝒪m2/∂ω𝒪m2)
                push!(freqpertsmatrix, (l, m, n, a, δω))
            end
        end
    end
    CSV.write("/home/dgw763/Documents/LinearizedSpin/ScalarPerturbations/freqpertsmatrix.csv", freqpertsmatrix)
end


@testset "KerrNewmanpertm2Fitting" begin

    dqHfile = "./OperatorShifts/KerrNewman/sm2KNTeukolskyChargeDeriv.csv"
    dwHfile = "./OperatorShifts/KerrNewman/sm2KNTeukolskyFreqDeriv.csv"

    ∂qH = OperatorShift(dqHfile)
    ∂ωH = OperatorShift(dwHfile)

    println("Made operator shifts")

    modelist = [[2,2,0],[2,2,1],[2,2,2],[2,1,0],[2,1,1],[2,0,0],[2,0,1],[3,2,0],[3,2,1],[3,3,0],[3,3,1],[4,4,0],[4,2,0]]
    for elem in modelist
        l = elem[1]
        m = elem[2]
        n = elem[3]

        freqpertsmatrix =DataFrame(l = Int[], m = Int[], n = Int[], a = Float64[] , δω = ComplexF64[])
        filenamestart="/home/dgw763/Documents/KerrNewman/Fitting/freqpertsmatrix"
        filenameend= ".csv"
        lstr=string(l)
        mstr=string(m)
        nstr=string(n)
        filename=filenamestart*lstr*mstr*nstr*filenameend
        println(filename)

        for a in range(0, stop=0.99, step=0.01)
        ψ = qnmfunctionnew(-2,l,m,n,a)

            # Compile ψ
            ψ(1,.5)
            println("Past ψ compile")

            ## Define the useful contours
            r₊ = ψ.R.r₊ ; r₋ = ψ.R.r₋ ; s = ψ.s ; Δr = 0.1*(r₊-r₋); ϵ = eps(0.1);

            # Define the Weight from Mark 2014
            weight = let a= ψ.a
                (r,z) ->(r-r₊)^s * (r-r₋)^s 
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
            ∂qHm2 = OperatorSandwich(ψ,∂qH,weight,ψ).Op
            ∂ωHm2 = OperatorSandwich(ψ,∂ωH,weight,ψ).Op
            println("Made Operators")

            ∂qℋm2 = Integrate(∂qHm2, TheContourup,abstol=1e-6)[1]
            ∂ωℋm2 = Integrate(∂ωHm2, TheContourup,abstol=1e-6)[1]

            @show δω = -(∂qℋm2/∂ωℋm2)

            push!(freqpertsmatrix, (l, m, n, a, δω))
        end
        CSV.write(filename, freqpertsmatrix)
    end
end

@testset "Eikonal" begin

    dwOfile = "/home/dgw763/Documents/LinearizedSpin/Operators/dwOcoefficients.csv"
    Hlinfile = "/home/dgw763/Documents/LinearizedSpin/Operators/Hlincoefficients.csv"
    Ilinfile = "/home/dgw763/Documents/LinearizedSpin/Operators/Ilincoefficients.csv"

    ∂ωO = OperatorShift(dwOfile)
    Hlin = OperatorShift(Hlinfile)
    Ilin = OperatorShift(Ilinfile)

    println("Made operator shifts")

    freqpertsmatrix =DataFrame(l = Int[], m = Int[], n = Int[], a = Float64[] , ωup = ComplexF64[], γup = ComplexF64[], ωdown = ComplexF64[], γdown = ComplexF64[])
    lmax=8

    for l in 2:lmax
        m=l
        n=2
        for a in [0.]
            ψ = qnmfunctionnew(-2,l,m,n,a)
            ψconj = qnmfunctionnew(-2,l,m,n,a,is_conjugate=true)
            ψm = qnmfunctionnew(-2,l,m,n,a,is_minus=true)
            ψmconj = qnmfunctionnew(-2,l,m,n,a,is_minus=true,is_conjugate=true)
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
            ∂ωOplusSchw = OperatorSandwich(ψ,∂ωO,weight,ψ).Op
            ∂ωOminusSchw = OperatorSandwich(ψm,∂ωO,weight,ψm).Op
            println("Done ∂ωO's")
            HlinplusSchw = OperatorSandwich(ψ,Hlin,weight,ψ).Op
            HlinminusSchw = OperatorSandwich(ψm,Hlin,weight,ψm).Op
            println("Done Hlin's")
            IlinplusSchw = OperatorSandwich(ψ,Ilin,weight,ψmconj).Op
            IlinminusSchw = OperatorSandwich(ψm,Ilin,weight,ψconj).Op
            println("Done I's")
            println("Made Operators")

            pert_ϵ=1
            ∂ω𝒪plusSchw= Integrate(∂ωOplusSchw, TheContourup,pertparam=pert_ϵ,abstol=1e-6)[1]
            ∂ω𝒪minusSchw= conj(Integrate(∂ωOminusSchw, TheContourdown,pertparam=pert_ϵ,abstol=1e-6)[1])
            ℋlinplusSchw= Integrate(HlinplusSchw, TheContourup,pertparam=pert_ϵ,abstol=1e-6)[1]
            ℋlinminusSchw= conj(Integrate(HlinminusSchw, TheContourdown,pertparam=pert_ϵ,abstol=1e-6)[1])
            ℐlinplusSchw= Integrate(IlinplusSchw, TheContourup,pertparam=pert_ϵ,abstol=1e-6)[1]
            ℐlinminusSchw= conj(Integrate(IlinminusSchw, TheContourdown,pertparam=pert_ϵ,abstol=1e-6)[1])
            ω2s=Computeω2(∂ω𝒪plusSchw,∂ω𝒪minusSchw,ℋlinplusSchw,ℋlinminusSchw,ℐlinplusSchw,ℐlinminusSchw,ψ)
            @show ω2s

            push!(freqpertsmatrix, (l, m, n, a, ω2s[1], ω2s[3], ω2s[2], ω2s[4]))
        end
    end

    CSV.write("/home/dgw763/Documents/LinearizedSpin/SpectralitySplit/Eikonaln2.csv", freqpertsmatrix)
end

@testset "LinearizedJP221" begin

    dwOfile = "/home/dgw763/Documents/LinearizedSpin/Operators/dwOcoefficients.csv"
    Hlinfile = "/home/dgw763/Documents/LinearizedSpin/Operators/Hlincoefficients.csv"
    Ilinfile = "/home/dgw763/Documents/LinearizedSpin/Operators/Ilincoefficients.csv"

    ∂ωO = OperatorShift(dwOfile)
    Hlin = OperatorShift(Hlinfile)
    Ilin = OperatorShift(Ilinfile)

    println("Made operator shifts")

    l=2
    m=2
    n=1
    a=0.
    ψ = qnmfunctionnew(-2,l,m,n,a)
    ψconj = qnmfunctionnew(-2,l,m,n,a,is_conjugate=true)
    ψm = qnmfunctionnew(-2,l,m,n,a,is_minus=true)
    ψmconj = qnmfunctionnew(-2,l,m,n,a,is_minus=true,is_conjugate=true)
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
    ∂ωOplusSchw = OperatorSandwich(ψ,∂ωO,weight,ψ).Op
    ∂ωOminusSchw = OperatorSandwich(ψm,∂ωO,weight,ψm).Op
    println("Done ∂ωO's")
    HlinplusSchw = OperatorSandwich(ψ,Hlin,weight,ψ).Op
    HlinminusSchw = OperatorSandwich(ψm,Hlin,weight,ψm).Op
    println("Done Hlin's")
    IlinplusSchw = OperatorSandwich(ψ,Ilin,weight,ψmconj).Op
    IlinminusSchw = OperatorSandwich(ψm,Ilin,weight,ψconj).Op
    println("Done I's")
    println("Made Operators")
    
    pert_ϵ=1
    ∂ω𝒪plusSchw= Integrate(∂ωOplusSchw, TheContourup,pertparam=pert_ϵ,abstol=1e-7)[1]
    ∂ω𝒪minusSchw= conj(Integrate(∂ωOminusSchw, TheContourdown,pertparam=pert_ϵ,abstol=1e-7)[1])
    ℋlinplusSchw= Integrate(HlinplusSchw, TheContourup,pertparam=pert_ϵ,abstol=1e-7)[1]
    ℋlinminusSchw= conj(Integrate(HlinminusSchw, TheContourdown,pertparam=pert_ϵ,abstol=1e-7)[1])
    ℐlinplusSchw= Integrate(IlinplusSchw, TheContourup,pertparam=pert_ϵ,abstol=1e-7)[1]
    ℐlinminusSchw= conj(Integrate(IlinminusSchw, TheContourdown,pertparam=pert_ϵ,abstol=1e-7)[1])

    @show ∂ω𝒪plusSchw
    @show ∂ω𝒪minusSchw
    @show ℋlinplusSchw
    @show ℋlinminusSchw
    @show ℐlinplusSchw
    @show ℐlinminusSchw

    ω2s=Computeω2(∂ω𝒪plusSchw,∂ω𝒪minusSchw,ℋlinplusSchw,ℋlinminusSchw,ℐlinplusSchw,ℐlinminusSchw,ψ)
    @show ω2s
    
end

@testset "LinearizedJPMany" begin

    dwOfile = "/home/dgw763/Documents/LinearizedSpin/Operators/dwOcoefficients.csv"
    Hlinfile = "/home/dgw763/Documents/LinearizedSpin/Operators/Hlincoefficients.csv"
    Ilinfile = "/home/dgw763/Documents/LinearizedSpin/Operators/Ilincoefficients.csv"

    ∂ωO = OperatorShift(dwOfile)
    Hlin = OperatorShift(Hlinfile)
    Ilin = OperatorShift(Ilinfile)

    println("Made operator shifts")

    # freqpertsmatrix =DataFrame(l = Int[], m = Int[], n = Int[], a = Float64[] , Reωup = Float64[], Imωup = Float64[], Reωdown = Float64[], Imωdown = Float64[], CoarseReωup = Float64[], CoarseImωup = Float64[], CoarseReωdown = Float64[], CoarseImωdown = Float64[])
    # freqpertsmatrix =DataFrame(l = Int[], m = Int[], n = Int[], a = Float64[] , Reωup = Float64[], Imωup = Float64[], Reωdown = Float64[], Imωdown = Float64[], γup = ComplexF64[], γdown = ComplexF64[])
    freqpertsmatrix =DataFrame(l = Int[], m = Int[], n = Int[], a = Float64[] , ωup = ComplexF64[], γup = ComplexF64[], ωdown = ComplexF64[], γdown = ComplexF64[])
    lmax=3

    for l in 2:lmax
        for m in 2:min(l,3)
            for n in 0:1
                for a in [0.0, 0.01, 0.02]
                    ψ = qnmfunctionnew(-2,l,m,n,a)
                    ψconj = qnmfunctionnew(-2,l,m,n,a,is_conjugate=true)

                    ψm = qnmfunctionnew(-2,l,m,n,a,is_minus=true)
                    ψmconj = qnmfunctionnew(-2,l,m,n,a,is_minus=true,is_conjugate=true)

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

                    ∂ωOplusSchw = OperatorSandwich(ψ,∂ωO,weight,ψ).Op
                    ∂ωOminusSchw = OperatorSandwich(ψm,∂ωO,weight,ψm).Op
                    println("Done ∂ωO's")
                    HlinplusSchw = OperatorSandwich(ψ,Hlin,weight,ψ).Op
                    HlinminusSchw = OperatorSandwich(ψm,Hlin,weight,ψm).Op
                    println("Done Hlin's")
                    IlinplusSchw = OperatorSandwich(ψ,Ilin,weight,ψmconj).Op
                    IlinminusSchw = OperatorSandwich(ψm,Ilin,weight,ψconj).Op
                    println("Done I's")

                    println("Made Operators")
                    
                    pert_ϵ=1

                    ∂ω𝒪plusSchw= Integrate(∂ωOplusSchw, TheContourup,pertparam=pert_ϵ,abstol=1e-7)[1]
                    ∂ω𝒪minusSchw= conj(Integrate(∂ωOminusSchw, TheContourdown,pertparam=pert_ϵ,abstol=1e-7)[1])
                    ℋlinplusSchw= Integrate(HlinplusSchw, TheContourup,pertparam=pert_ϵ,abstol=1e-7)[1]
                    ℋlinminusSchw= conj(Integrate(HlinminusSchw, TheContourdown,pertparam=pert_ϵ,abstol=1e-7)[1])
                    ℐlinplusSchw= Integrate(IlinplusSchw, TheContourup,pertparam=pert_ϵ,abstol=1e-7)[1]
                    ℐlinminusSchw= conj(Integrate(IlinminusSchw, TheContourdown,pertparam=pert_ϵ,abstol=1e-7)[1])
                    ω2s=Computeω2(∂ω𝒪plusSchw,∂ω𝒪minusSchw,ℋlinplusSchw,ℋlinminusSchw,ℐlinplusSchw,ℐlinminusSchw,ψ)
                    @show ω2s

                    # Coarse∂ω𝒪plusSchw= Integrate(∂ωOplusSchw, TheContourup,pertparam=pert_ϵ,abstol=1e-6)[1]
                    # Coarse∂ω𝒪minusSchw= conj(Integrate(∂ωOminusSchw, TheContourdown,pertparam=pert_ϵ,abstol=1e-6)[1])
                    # CoarseℋlinplusSchw= Integrate(HlinplusSchw, TheContourup,pertparam=pert_ϵ,abstol=1e-6)[1]
                    # CoarseℋlinminusSchw= conj(Integrate(HlinminusSchw, TheContourdown,pertparam=pert_ϵ,abstol=1e-6)[1])
                    # CoarseℐlinplusSchw= Integrate(IlinplusSchw, TheContourup,pertparam=pert_ϵ,abstol=1e-6)[1]
                    # CoarseℐlinminusSchw= conj(Integrate(IlinminusSchw, TheContourdown,pertparam=pert_ϵ,abstol=1e-6)[1])

                    # Coarseω2s=Computeω2(Coarse∂ω𝒪plusSchw,Coarse∂ω𝒪minusSchw,CoarseℋlinplusSchw,CoarseℋlinminusSchw,CoarseℐlinplusSchw,CoarseℐlinminusSchw,ψ)
                    # @show Coarseω2s

                    # push!(freqpertsmatrix, (l, m, n, a, real(ω2s[1]),imag(ω2s[1]),real(ω2s[2]),imag(ω2s[2]),real(Coarseω2s[1]),imag(Coarseω2s[1]),real(Coarseω2s[2]),imag(Coarseω2s[2])))
                    push!(freqpertsmatrix, (l, m, n, a, ω2s[1], ω2s[3], ω2s[2], ω2s[4]))
                end
            end
        end
    end

    CSV.write("/home/dgw763/Documents/LinearizedSpin/SpectralitySplit/Final/freqpertsmatrix.csv", freqpertsmatrix)
end

@testset "LinearizedJP" begin     #CHECKED AND GOOD
    pert_ϵ=1
    
    ψ = qnmfunctionnew(-2,2,2,0,0.01)
    ψconj = qnmfunctionnew(-2,2,2,0,0.01,is_conjugate=true)

    ψm = qnmfunctionnew(-2,2,2,0,0.01,is_minus=true)
    ψmconj = qnmfunctionnew(-2,2,2,0,0.01,is_minus=true,is_conjugate=true)

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

    
    dwOfile = "./OperatorShifts/LinearizedJP/dwOcoefficients.csv"
    Hlinfile = "./OperatorShifts/LinearizedJP/Hlincoefficients.csv"
    Ilinfile = "./OperatorShifts/LinearizedJP/Ilincoefficients.csv"

    ∂ωO = OperatorShift(dwOfile)
    Hlin = OperatorShift(Hlinfile)
    Ilin = OperatorShift(Ilinfile)

    println("Made operator shifts")
    

    ∂ωOplus = OperatorSandwich(ψ,∂ωO,weight,ψ).Op
    ∂ωOminus = OperatorSandwich(ψm,∂ωO,weight,ψm).Op

    Hlinplus = OperatorSandwich(ψ,Hlin,weight,ψ).Op
    Hlinminus = OperatorSandwich(ψm,Hlin,weight,ψm).Op
    Ilinplus = OperatorSandwich(ψ,Ilin,weight,ψmconj).Op
    Ilinminus = OperatorSandwich(ψm,Ilin,weight,ψconj).Op

    println("Made Operators")
    
    ∂ω𝒪plus= Integrate(∂ωOplus, TheContourup,pertparam=pert_ϵ,abstol=1e-6)[1]
    ∂ω𝒪minus= conj(Integrate(∂ωOminus, TheContourdown,pertparam=pert_ϵ,abstol=1e-6)[1])
    ℋlinplus= Integrate(Hlinplus, TheContourup,pertparam=pert_ϵ,abstol=1e-6)[1]
    ℋlinminus= conj(Integrate(Hlinminus, TheContourdown,pertparam=pert_ϵ,abstol=1e-6)[1])
    ℐlinplus= Integrate(Ilinplus, TheContourup,pertparam=pert_ϵ,abstol=1e-6)[1]
    ℐlinminus= conj(Integrate(Ilinminus, TheContourdown,pertparam=pert_ϵ,abstol=1e-6)[1])
    
    ω2s=Computeω2(∂ω𝒪plus,∂ω𝒪minus,ℋlinplus,ℋlinminus,ℐlinplus,ℐlinminus,ψ)
    @show ω2s

end

@testset "KerrNewmanpertp2" begin   #CHECKED AND GOOD

    dqHfile = "./OperatorShifts/KerrNewman/sp2KNTeukolskyChargeDeriv.csv"
    dwHfile = "./OperatorShifts/KerrNewman/sp2KNTeukolskyFreqDeriv.csv"

    ∂qH = OperatorShift(dqHfile)
    ∂ωH = OperatorShift(dwHfile)

    println("Made operator shifts")
    
    ψ = qnmfunctionnew(2,2,2,0,0.1)

    # Compile ψ
    ψ(1,.5)
    println("Past ψ compile")

    ## Define the useful contours
    r₊ = ψ.R.r₊ ; r₋ = ψ.R.r₋ ; s = ψ.s ; Δr = 0.1*(r₊-r₋); ϵ = eps(0.1);

    # Define the Weight from Mark 2014
    weight = let a= ψ.a
        (r,z) ->(r-r₊)^s * (r-r₋)^s  
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
    ∂qHp2 = OperatorSandwich(ψ,∂qH,weight,ψ).Op
    ∂ωHp2 = OperatorSandwich(ψ,∂ωH,weight,ψ).Op
    println("Made Operators")

    ∂qℋp2 = Integrate(∂qHp2, TheContourup,abstol=1e-6)[1]
    ∂ωℋp2 = Integrate(∂ωHp2, TheContourup,abstol=1e-6)[1]

    @show δω = -(∂qℋp2/∂ωℋp2)
end

@testset "KerrNewmanpertm2" begin   #CHECKED AND GOOD

    dqHfile = "./OperatorShifts/KerrNewman/sm2KNTeukolskyChargeDeriv.csv"
    dwHfile = "./OperatorShifts/KerrNewman/sm2KNTeukolskyFreqDeriv.csv"

    ∂qH = OperatorShift(dqHfile)
    ∂ωH = OperatorShift(dwHfile)

    println("Made operator shifts")
    
    ψ = qnmfunctionnew(-2,2,2,0,0.1)

    # Compile ψ
    ψ(1,.5)
    println("Past ψ compile")

    ## Define the useful contours
    r₊ = ψ.R.r₊ ; r₋ = ψ.R.r₋ ; s = ψ.s ; Δr = 0.1*(r₊-r₋); ϵ = eps(0.1);

    # Define the Weight from Mark 2014
    weight = let a= ψ.a
        (r,z) ->(r-r₊)^s * (r-r₋)^s 
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
    ∂qHm2 = OperatorSandwich(ψ,∂qH,weight,ψ).Op
    ∂ωHm2 = OperatorSandwich(ψ,∂ωH,weight,ψ).Op
    println("Made Operators")

    ∂qℋm2 = Integrate(∂qHm2, TheContourup,abstol=1e-6)[1]
    ∂ωℋm2 = Integrate(∂ωHm2, TheContourup,abstol=1e-6)[1]

    @show δω = -(∂qℋm2/∂ωℋm2)
end

@testset "Kerr as Perturbation to Schwarzschild" begin #CHECKED AND GOOD
    pert_a=1
    
    ψ = qnmfunctionnew(-2,2,2,0,0.)
    ψconj = qnmfunctionnew(-2,2,2,0,0.,is_conjugate=true)

    ψm = qnmfunctionnew(-2,2,2,0,0.,is_minus=true)
    ψmconj = qnmfunctionnew(-2,2,2,0,0.,is_minus=true,is_conjugate=true)


    # Compile ψ
    ψ(1,.5)
    println("Past ψ compile")

    # Weight for integration
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

    dwOfile = "./OperatorShifts/Schwarzschild/dwOcoefficients.csv"
    Hfile = "./OperatorShifts/Schwarzschild/Hcoefficients.csv"
    Ifile = "./OperatorShifts/Schwarzschild/Icoefficients.csv"

    ∂ωO = OperatorShift(dwOfile)
    H = OperatorShift(Hfile)
    I = OperatorShift(Ifile)

    println("Made operator shifts")
    
    ∂ωOplusSchw = OperatorSandwich(ψ,∂ωO,weight,ψ).Op
    ∂ωOminusSchw = OperatorSandwich(ψm,∂ωO,weight,ψm).Op
    HplusSchw = OperatorSandwich(ψ,H,weight,ψ).Op
    HminusSchw = OperatorSandwich(ψm,H,weight,ψm).Op
    IplusSchw = OperatorSandwich(ψ,I,weight,ψmconj).Op
    IminusSchw = OperatorSandwich(ψm,I,weight,ψconj).Op

    println("Made Operators")

    𝒪plusSchw= Integrate(OSchw1, TheContourup,pertparam=pert_a,abstol=1e-6)[1]
    𝒪minusSchw= conj(Integrate(OSchw2, TheContourdown,pertparam=pert_a,abstol=1e-6)[1])
    ∂ω𝒪plusSchw= Integrate(∂ωOplusSchw, TheContourup,pertparam=pert_a,abstol=1e-6)[1]
    ∂ω𝒪minusSchw= conj(Integrate(∂ωOminusSchw, TheContourdown,pertparam=pert_a,abstol=1e-6)[1])
    ℋplusSchw= Integrate(HplusSchw, TheContourup,pertparam=pert_a,abstol=1e-6)[1]
    ℋminusSchw= conj(Integrate(HminusSchw, TheContourdown,pertparam=pert_a,abstol=1e-6)[1])
    ℐplusSchw= Integrate(IplusSchw, TheContourup,pertparam=pert_a,abstol=1e-6)[1]
    ℐminusSchw= conj(Integrate(IminusSchw, TheContourdown,pertparam=pert_a,abstol=1e-6)[1])
    
    ω2s=Computeω2(∂ω𝒪plusSchw,∂ω𝒪minusSchw,ℋplusSchw,ℋminusSchw,ℐplusSchw,ℐminusSchw,ψ)
    @show ω2s

end