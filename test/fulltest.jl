using ContourIntegrals
using KerrQNMShifts
using KerrQuasinormalModes
using Test
using CSV, DataFrames

println("Done usings")

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

@testset "LinearizedJPSingle" begin
    pert_ϵ=1
    
    ψ = qnmfunctionnew(-2,2,2,0,0.1)
    ψconj = qnmfunctionnew(-2,2,2,0,0.1,is_conjugate=true)

    ψm = qnmfunctionnew(-2,2,2,0,0.1,is_minus=true)
    ψmconj = qnmfunctionnew(-2,2,2,0,0.1,is_minus=true,is_conjugate=true)

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

    # Ofile = "/home/dgw763/Documents/LinearizedSpin/Operators/Ocoefficients.csv"
    # Olinfile = "/home/dgw763/Documents/LinearizedSpin/Operators/Olincoefficients.csv"
    dwOfile = "/home/dgw763/Documents/LinearizedSpin/Operators/dwOcoefficients.csv"
    # dwOlinfile = "/home/dgw763/Documents/LinearizedSpin/Operators/dwOlincoefficients.csv"
    Hlinfile = "/home/dgw763/Documents/LinearizedSpin/Operators/Hlincoefficients.csv"
    Ilinfile = "/home/dgw763/Documents/LinearizedSpin/Operators/Ilincoefficients.csv"

    # O = OperatorShift(Ofile)
    # Olin = OperatorShift(Olinfile)
    ∂ωO = OperatorShift(dwOfile)
    # ∂ωOlin = OperatorShift(dwOlinfile)
    Hlin = OperatorShift(Hlinfile)
    Ilin = OperatorShift(Ilinfile)

    println("Made operator shifts")
    
    # Oplus = OperatorSandwich(ψ,O,weight,ψ).Op
    # Ominus = OperatorSandwich(ψm,O,weight,ψm).Op

    # Olinplus = OperatorSandwich(ψ,Olin,weight,ψ).Op
    # Olinminus = OperatorSandwich(ψm,Olin,weight,ψm).Op

    ∂ωOplus = OperatorSandwich(ψ,∂ωO,weight,ψ).Op
    ∂ωOminus = OperatorSandwich(ψm,∂ωO,weight,ψm).Op
    
    # ∂ωOlinplus = OperatorSandwich(ψ,∂ωOlin,weight,ψ).Op
    # ∂ωOlinminus = OperatorSandwich(ψm,∂ωOlin,weight,ψm).Op

    Hlinplus = OperatorSandwich(ψ,Hlin,weight,ψ).Op
    Hlinminus = OperatorSandwich(ψm,Hlin,weight,ψm).Op
    Ilinplus = OperatorSandwich(ψ,Ilin,weight,ψmconj).Op
    Ilinminus = OperatorSandwich(ψm,Ilin,weight,ψconj).Op

    println("Made Operators")

    # @show Oplus(8+im,0.6,pertparam=pert_ϵ)
    # @show Ominus(8+im,0.6,pertparam=pert_ϵ)

    # @show Olinplus(8+im,0.6,pertparam=pert_ϵ)
    # @show Olinminus(8+im,0.6,pertparam=pert_ϵ)

    @show ∂ωOplus(8+im,0.6,pertparam=pert_ϵ)
    @show ∂ωOminus(8+im,0.6,pertparam=pert_ϵ)

    @show ∂ωOlinplus(8+im,0.6,pertparam=pert_ϵ)
    @show ∂ωOlinminus(8+im,0.6,pertparam=pert_ϵ)

    @show Hlinplus(8+im,0.6,pertparam=pert_ϵ)
    @show Hlinminus(8+im,0.6,pertparam=pert_ϵ)
    @show Ilinplus(8+im,0.6,pertparam=pert_ϵ)
    @show Ilinminus(8+im,0.6,pertparam=pert_ϵ)

    println("Complied Operators")

    # 𝒪plus= Integrate(Oplus, TheContourup,pertparam=pert_ϵ,abstol=1e-6)[1]
    # 𝒪minus= conj(Integrate(Ominus, TheContourdown,pertparam=pert_ϵ,abstol=1e-6)[1])

    # 𝒪linplus= Integrate(Olinplus, TheContourup,pertparam=pert_ϵ,abstol=1e-6)[1]
    # 𝒪linminus= conj(Integrate(Olinminus, TheContourdown,pertparam=pert_ϵ,abstol=1e-6)[1])
    
    ∂ω𝒪plus= Integrate(∂ωOplus, TheContourup,pertparam=pert_ϵ,abstol=1e-6)[1]
    ∂ω𝒪minus= conj(Integrate(∂ωOminus, TheContourdown,pertparam=pert_ϵ,abstol=1e-6)[1])
    ∂ω𝒪linplus= Integrate(∂ωOlinplus, TheContourup,pertparam=pert_ϵ,abstol=1e-6)[1]
    ∂ω𝒪linminus= conj(Integrate(∂ωOlinminus, TheContourdown,pertparam=pert_ϵ,abstol=1e-6)[1])
    ℋlinplus= Integrate(Hlinplus, TheContourup,pertparam=pert_ϵ,abstol=1e-6)[1]
    ℋlinminus= conj(Integrate(Hlinminus, TheContourdown,pertparam=pert_ϵ,abstol=1e-6)[1])
    ℐlinplus= Integrate(Ilinplus, TheContourup,pertparam=pert_ϵ,abstol=1e-6)[1]
    ℐlinminus= conj(Integrate(Ilinminus, TheContourdown,pertparam=pert_ϵ,abstol=1e-6)[1])

    # @show ∂ω𝒪plus
    # @show ∂ω𝒪minus
    # @show ℋplus
    # @show ℋminus
    # @show ℐplus
    # @show ℐminus
    
    ω2s=Computeω2(∂ω𝒪plus,∂ω𝒪minus,ℋlinplus,ℋlinminus,ℐlinplus,ℐlinminus,ψ)
    @show ω2s

end

@testset "KerrNewmanpertp2" begin

    dqHfile = "C:/Users/dwuuu/Documents/UT Academics/Research/RingdownDataAnalysis/ChargedQNMs/dqHp2coefficients.csv"
    dwHfile = "C:/Users/dwuuu/Documents/UT Academics/Research/RingdownDataAnalysis/ChargedQNMs/dwHp2coefficients.csv"

    ∂qH = OperatorShift(dqHfile)
    ∂ωH = OperatorShift(dwHfile)

    println("Made operator shifts")
    
    ψ = qnmfunctionnew(2,2,2,0,0.)

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

    @show ∂qHp2(8+im,0.6)
    @show ∂ωHp2(8+im,0.6)

    println("Complied Operators")

    ∂qℋp2 = Integrate(∂qHp2, TheContourup,abstol=1e-6)[1]
    ∂ωℋp2 = Integrate(∂ωHp2, TheContourup,abstol=1e-6)[1]

    @show δω = -(∂qℋp2/∂ωℋp2)
end

@testset "KerrNewmanpertm2" begin

    dqHfile = "C:/Users/dwuuu/Documents/UT Academics/Research/RingdownDataAnalysis/ChargedQNMs/dqHm2coefficients.csv"
    dwHfile = "C:/Users/dwuuu/Documents/UT Academics/Research/RingdownDataAnalysis/ChargedQNMs/dwHm2coefficients.csv"

    ∂qH = OperatorShift(dqHfile)
    ∂ωH = OperatorShift(dwHfile)

    println("Made operator shifts")
    
    ψ = qnmfunctionnew(-2,2,2,0,0.)

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

    @show ∂qHm2(8+im,0.6)
    @show ∂ωHm2(8+im,0.6)

    println("Complied Operators")

    ∂qℋm2 = Integrate(∂qHm2, TheContourup,abstol=1e-6)[1]
    ∂ωℋm2 = Integrate(∂ωHm2, TheContourup,abstol=1e-6)[1]

    @show δω = -(∂qℋm2/∂ωℋm2)
end

@testset "FinalTestMany" begin
    
    dwOfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/Final/dwOcoefficients.csv"
    Hfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/Final/Hcoefficients.csv"
    Ifile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/Final/Icoefficients.csv"
    Ofile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/Final/Ocoefficients.csv"

    ∂ωO = OperatorShift(dwOfile)
    H = OperatorShift(Hfile)
    I = OperatorShift(Ifile)
    O = OperatorShift(Ofile)

    println("Made operator shifts")


    freqpertsmatrix =DataFrame(l = Int[], m = Int[], n = Int[], a = Float64[], value = ComplexF64[])
    lmax=5

    for l in 2:lmax
        for m in 2:min(l,3)
            for n in 0:2
                ψ = qnmfunctionnew(-2,l,m,n,0.)
                ψconj = qnmfunctionnew(-2,l,m,n,0.,is_conjugate=true)

                ψm = qnmfunctionnew(-2,l,m,n,0.,is_minus=true)
                ψmconj = qnmfunctionnew(-2,l,m,n,0.,is_minus=true,is_conjugate=true)

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

                OplusSchw = OperatorSandwich(ψ,O,weight,ψ).Op

                ∂ωOplusSchw = OperatorSandwich(ψ,∂ωO,weight,ψ).Op
                ∂ωOminusSchw = OperatorSandwich(ψm,∂ωO,weight,ψm).Op
                println("Done ∂ωO's")
                HplusSchw = OperatorSandwich(ψ,H,weight,ψ).Op
                HminusSchw = OperatorSandwich(ψm,H,weight,ψm).Op
                println("Done H's")
                IplusSchw = OperatorSandwich(ψ,I,weight,ψmconj).Op
                IminusSchw = OperatorSandwich(ψm,I,weight,ψconj).Op
                println("Done I's")

                println("Made Operators")

                for pert_a_iter in 1:3

                    pert_a=.02*pert_a_iter

                    OplusSchw(8+im,0.6,pertparam=.1)

                    ∂ω𝒪plusSchw= Integrate(∂ωOplusSchw, TheContourup,pertparam=pert_a,abstol=1e-6)[1]
                    ∂ω𝒪minusSchw= conj(Integrate(∂ωOminusSchw, TheContourdown,pertparam=pert_a,abstol=1e-6)[1])
                    ℋplusSchw= Integrate(HplusSchw, TheContourup,pertparam=pert_a,abstol=1e-6)[1]
                    ℋminusSchw= conj(Integrate(HminusSchw, TheContourdown,pertparam=pert_a,abstol=1e-6)[1])
                    ℐplusSchw= Integrate(IplusSchw, TheContourup,pertparam=pert_a,abstol=1e-6)[1]
                    ℐminusSchw= conj(Integrate(IminusSchw, TheContourdown,pertparam=pert_a,abstol=1e-6)[1])

                    ω2s=Computeω2(∂ω𝒪plusSchw,∂ω𝒪minusSchw,ℋplusSchw,ℋminusSchw,ℐplusSchw,ℐminusSchw,ψ)
                    @show ω2s[1],pert_a

                    push!(freqpertsmatrix, (l, m, n, pert_a, ω2s[1]))
                end
            end
        end
    end

    CSV.write("C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/Final/juliafreqpertsmatrix.csv", freqpertsmatrix)

end

@testset "FinalTestSingle" begin
    
    dwOfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/Final/dwOcoefficients.csv"
    Hfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/Final/Hcoefficients.csv"
    Ifile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/Final/Icoefficients.csv"
    Ofile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/Final/Ocoefficients.csv"

    ∂ωO = OperatorShift(dwOfile)
    H = OperatorShift(Hfile)
    I = OperatorShift(Ifile)
    O = OperatorShift(Ofile)

    println("Made operator shifts")

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

    OplusSchw = OperatorSandwich(ψ,O,weight,ψ).Op
    OminusSchw = OperatorSandwich(ψm,O,weight,ψm).Op

    ∂ωOplusSchw = OperatorSandwich(ψ,∂ωO,weight,ψ).Op
    ∂ωOminusSchw = OperatorSandwich(ψm,∂ωO,weight,ψm).Op
    println("Done ∂ωO's")
    HplusSchw = OperatorSandwich(ψ,H,weight,ψ).Op
    HminusSchw = OperatorSandwich(ψm,H,weight,ψm).Op
    println("Done H's")
    IplusSchw = OperatorSandwich(ψ,I,weight,ψmconj).Op
    IminusSchw = OperatorSandwich(ψm,I,weight,ψconj).Op
    println("Done I's")

    println("Made Operators")

    pert_a=.1

    OplusSchw(8+im,0.6,pertparam=pert_a)
    OminusSchw(8+im,0.6,pertparam=pert_a)

    𝒪plusSchw= Integrate(OplusSchw, TheContourup,pertparam=pert_a,abstol=1e-6)[1]
    𝒪minusSchw= conj(Integrate(OminusSchw, TheContourdown,pertparam=pert_a,abstol=1e-6)[1])
    ∂ω𝒪plusSchw= Integrate(∂ωOplusSchw, TheContourup,pertparam=pert_a,abstol=1e-6)[1]
    ∂ω𝒪minusSchw= conj(Integrate(∂ωOminusSchw, TheContourdown,pertparam=pert_a,abstol=1e-6)[1])
    ℋplusSchw= Integrate(HplusSchw, TheContourup,pertparam=pert_a,abstol=1e-6)[1]
    ℋminusSchw= conj(Integrate(HminusSchw, TheContourdown,pertparam=pert_a,abstol=1e-6)[1])
    ℐplusSchw= Integrate(IplusSchw, TheContourup,pertparam=pert_a,abstol=1e-6)[1]
    ℐminusSchw= conj(Integrate(IminusSchw, TheContourdown,pertparam=pert_a,abstol=1e-6)[1])

    ω2s=Computeω2(∂ω𝒪plusSchw,∂ω𝒪minusSchw,ℋplusSchw,ℋminusSchw,ℐplusSchw,ℐminusSchw,ψ)
    @show ω2s[1],pert_a
end

@testset "Weyl a Linearity" begin
    r=2.1
    z=.5
    lmax=3
    for l in 2:lmax
        for m in 2:min(l,3)
            for n in 0:1
                ψ0 = qnmfunctionnew(-2,l,m,n,0.)
                ψsmall = qnmfunctionnew(-2,l,m,n,.00001)
                # slopefine=(ψ0(r,z)-ψsmall(r,z))/.00001
                slopefine=(ψ0.R(r)-ψsmall.R(r))/.00001
                @show slopefine
                for a in [0.01, 0.02]
                     println(l,", ",m,", ",n,", ",a)
                     ψ = qnmfunctionnew(-2,l,m,n,a)
                    #  slopecoarse = (ψ0(r,z)-ψ(r,z))/a
                     slopecoarse = (ψ0.R(r)-ψ.R(r))/a
                     @show slopecoarse
                     # ψconj = qnmfunctionnew(-2,l,m,n,a,is_conjugate=true)
                     # ψm = qnmfunctionnew(-2,l,m,n,a,is_minus=true)
                     # ψmconj = qnmfunctionnew(-2,l,m,n,a,is_minus=true,is_conjugate=true)
                     # Compile ψ
                end

                # slopefine=(ψ0(r,z)-ψsmall(r,z))/.00001
                slopefine=(ψ0.S(z)-ψsmall.S(z))/.00001
                @show slopefine
                for a in [0.05, 0.1]
                     println(l,", ",m,", ",n,", ",a)
                     ψ = qnmfunctionnew(-2,l,m,n,a)
                    #  slopecoarse = (ψ0(r,z)-ψ(r,z))/a
                     slopecoarse = (ψ0.S(z)-ψ.S(z))/a
                     @show slopecoarse
                end
            end
        end
    end
end

@testset "integration" begin
    abstract type Callable end
    abstract type CallableAtom <: Callable end
    abstract type CallableCombination <: Callable end

    ### Radial Functions
    struct Testfunc{T} <: CallableAtom
        is_conjugate::Bool;is_minus::Bool;s::Int64; m::Int64; a::Float64; ω::Complex{Float64}
    end

    function (f::Testfunc)(r,z)
        log(im*r)*exp(Complex(r)^2)
    end

    f=Testfunc{Complex{Float64}}(false,false,1,1,1.,1.)

    import KerrQuasinormalModes: ∂r
    function ∂r(f::Testfunc)
        f
    end

    import KerrQuasinormalModes: ∂θ
    function ∂θ(f::Testfunc)
        f
    end

    weight = let
        (r,z) ->1
    end

     ## Define the useful contours
    ϵ = eps(0.1);

     #The upwards pointing contour
     point1up = .1-.1*im
     point2up = -.1-.1*
     
     angular = LineSegment(-1.0+100*ϵ , 1.0-100*ϵ , true) #to avoid the NaNs at the edges
 
     radial1up = SemiInfiniteLine(point1up , point1up + .1*im , false)
     C1up = radial1up ⊗ angular
 
     radial2up = LineSegment(point1up,point2up,true)
     C2up = radial2up ⊗ angular

     radial3up = SemiInfiniteLine(point2up , point2up + .1*im , true)
     C3up = radial3up ⊗ angular
 
     TheContourup = C1up⊕C2up⊕C3up

    unitfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/FullTest/UnitOperator.csv"
 
    unit = OperatorShift(unitfile)

    unittest = OperatorSandwich(weight,unit,weight,f).Op

    @show unittest(8+im,0.6,pertparam=0)

    println("Complied Operators")

    unitintegral= Integrate(unittest, TheContourup,pertparam=0,abstol=1e-6)[1]
    # Should be 11.1367 + 4.44089*10^-16 I which it looks like it is

end

@testset "Testfreqpert" begin
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

    dwOfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/FullTest/dwOcoefficients.csv"
    KerrδOfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/FullTest/KerrδOcoefficients.csv"

    ∂ωO = OperatorShift(dwOfile)
    KerrδO = OperatorShift(KerrδOfile)

    println("Made operator shifts")
    
    ∂ωOplusSchw = OperatorSandwich(ψ,∂ωO,weight,ψ).Op
    ∂ωOminusSchw = OperatorSandwich(ψm,∂ωO,weight,ψm).Op
    KerrδOplus = OperatorSandwich(ψ,KerrδO,weight,ψ).Op
    KerrδOminus = OperatorSandwich(ψm,KerrδO,weight,ψm).Op
    
    println("Made Operators")

    @show ∂ωOplusSchw(8+im,0.6,pertparam=pert_a)
    @show ∂ωOminusSchw(8+im,0.6,pertparam=pert_a)
    @show KerrδOplus(8+im,0.6,pertparam=pert_a)
    @show KerrδOminus(8+im,0.6,pertparam=pert_a)

    println("Complied Operators")

    ∂ω𝒪plusSchw= Integrate(∂ωOplusSchw, TheContourup,pertparam=pert_a,abstol=1e-6)[1]
    ∂ω𝒪minusSchw= conj(Integrate(∂ωOminusSchw, TheContourdown,pertparam=pert_a,abstol=1e-6)[1])
    δ𝒪plusKerr= Integrate(KerrδOplus, TheContourup,pertparam=pert_a,abstol=1e-6)[1]
    δ𝒪minusKerr= conj(Integrate(KerrδOminus, TheContourdown,pertparam=pert_a,abstol=1e-6)[1])

    @show ∂ω𝒪plusSchw
    @show ∂ω𝒪minusSchw
    @show δ𝒪plusKerr
    @show δ𝒪minusKerr

    # You don't need to multiply by pert_a here because it's already built into
    # the operator when you mutliplied by 𝒶 in Mathematica
    @show -(δ𝒪plusKerr/∂ω𝒪plusSchw)
    @show δ𝒪minusKerr/∂ω𝒪minusSchw
end

@testset "FullTest" begin
    pert_a=.01
    
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


    # Ofile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/FullTest/1_14_Old/Ocoefficients.csv"
    # dwOfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/FullTest/1_14_Old/dwOcoefficients.csv"
    # Hfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/FullTest/1_14_Old/Hcoefficients.csv"
    # Ifile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/FullTest/1_14_Old/Icoefficients.csv"
    Ofile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/Final/Ocoefficients.csv"
    dwOfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/Final/dwOcoefficients.csv"
    Hfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/Final/Hcoefficients.csv"
    Ifile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/Final/Icoefficients.csv"

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
    @show OSchw2(8+im,0.6,pertparam=pert_a)
    @show ∂ωOplusSchw(8+im,0.6,pertparam=pert_a)
    @show ∂ωOminusSchw(8+im,0.6,pertparam=pert_a)
    @show HplusSchw(8+im,0.6,pertparam=pert_a)
    @show HminusSchw(8+im,0.6,pertparam=pert_a)
    @show IplusSchw(8+im,0.6,pertparam=pert_a)
    @show IminusSchw(8+im,0.6,pertparam=pert_a)

    println("Complied Operators")
    # pert_a=0.05

    𝒪plusSchw= Integrate(OSchw1, TheContourup,pertparam=pert_a,abstol=1e-6)[1]
    𝒪minusSchw= conj(Integrate(OSchw2, TheContourdown,pertparam=pert_a,abstol=1e-6)[1])
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
    Computeω2(∂ω𝒪plusSchw,-∂ω𝒪plusSchw,ℋplusSchw,ℋplusSchw,ℐplusSchw,ℐplusSchw,ψ)
    @show ω2s

    # test for general γ
    ComputeDplus(ψ)*ℋplusSchw
    (12*im*ψ.ω)*ℐplusSchw

    ℐplustest=ComputeDplus(ψ)*ℋplusSchw/(12*im*ψ.ω)
    Computeω2(∂ω𝒪plusSchw,-∂ω𝒪plusSchw,ℋplusSchw,ℋplusSchw,ℐplustest,ℐplustest,ψ)

    ℋplustest=(12*im*ψ.ω)*ℐplusSchw/ComputeDplus(ψ)
    Computeω2(∂ω𝒪plusSchw,-∂ω𝒪plusSchw,ℋplustest,ℋplustest,ℐplusSchw,ℐplusSchw,ψ)


end


@testset "derivative test" begin
    pert_a=.1
    
    ψ = qnmfunctionnew(-2,2,2,0,0.1)
    ψs = qnmfunctionnew(2,2,2,0,0.1)
    ψm = qnmfunctionnew(-2,2,2,0,0.1,is_minus=true)
    ψsm = qnmfunctionnew(2,2,2,0,0.1,is_minus=true)
    ψconj = qnmfunctionnew(-2,2,2,0,0.1,is_conjugate=true)
    ψsconj = qnmfunctionnew(2,2,2,0,0.1,is_conjugate=true)
    ψmconj = qnmfunctionnew(-2,2,2,0,0.1,is_minus=true,is_conjugate=true)
    ψsmconj = qnmfunctionnew(2,2,2,0,0.1,is_minus=true,is_conjugate=true)

    ψ(3), ψ.R(3)
    ψm(3), ψm.R(3)
    ψconj(3), ψconj.R(3)
    ψmconj(3), ψmconj.R(3)

    ψs(3), ψs.R(3)
    ψsm(3),     ψsm.R(3)
    ψsconj(3),     ψsconj.R(3)
    ψsmconj(3),     ψsmconj.R(3)
    
    ψ(3+im), ψ.R(3+im)
    ψm(3+im),     ψm.R(3+im)
    ψconj(3+im),     ψconj.R(3+im)
    ψmconj(3+im),     ψmconj.R(3+im)
 
    ψs(3+im),     ψs.R(3+im)
    ψsm(3+im),     ψsm.R(3+im)
    ψsconj(3+im),     ψsconj.R(3+im)
    ψsmconj(3+im),     ψsmconj.R(3+im)


    ψ.S(.4)
    ψs.S(.4)
    ψm.S(.4)
    ψsm.S(.4)
    ψconj.S(.4)
    ψsconj.S(.4)
    ψmconj.S(.4)
    ψsmconj.S(.4)



    ∂r(ψ)(3+im,.4)-(ψ(3.001+1.001im,.4)-ψ(3+im,.4))/(.001+.001im)
    ∂r(ψ)(3,.4)-(ψ(3.001,.4)-ψ(3,.4))/(.001)

    ∂r(ψm)(3+im,.4)-(ψm(3.001+1.001im,.4)-ψm(3+im,.4))/(.001+.001im)
    ∂r(ψm)(3,.4)-(ψm(3.001,.4)-ψm(3,.4))/(.001)

    ∂r(ψconj)(3+im,.4)-(ψconj(3.001+1.001im,.4)-ψconj(3+im,.4))/(.001+.001im)
    ∂r(ψconj)(3,.4)-(ψconj(3.001,.4)-ψconj(3,.4))/(.001)

    ∂r(ψmconj)(3+im,.4)-(ψmconj(3.001+1.001im,.4)-ψmconj(3+im,.4))/(.001+.001im)
    ∂r(ψmconj)(3,.4)-(ψmconj(3.001)-ψmconj(3,.4))/(.001)

    ∂r(ψs)(3+im,.4)-(ψs(3.001+1.001im,.4)-ψs(3+im,.4))/(.001+.001im)
    ∂r(ψs)(3,.4)-(ψs(3.001,.4)-ψs(3,.4))/(.001)

    ∂r(ψsm)(3+im,.4)-(ψsm(3.001+1.001im,.4)-ψsm(3+im,.4))/(.001+.001im)
    ∂r(ψsm)(3,.4)-(ψsm(3.001,.4)-ψsm(3,.4))/(.001)

    ∂r(ψsconj)(3+im,.4)-(ψsconj(3.001+1.001im,.4)-ψsconj(3+im,.4))/(.001+.001im)
    ∂r(ψsconj)(3,.4)-(ψsconj(3.001,.4)-ψsconj(3,.4))/(.001)

    ∂r(ψsmconj)(3+im,.4)-(ψsmconj(3.001+1.001im,.4)-ψsmconj(3+im,.4))/(.001+.001im)
    ∂r(ψsmconj)(3,.4)-(ψsmconj(3.001,.4)-ψsmconj(3,.4))/(.001)


    ∂θ(ψ)(3+im,.4)-(ψ(3+im,.4001)-ψ(3+im,.4))/.0001*(-sqrt(1-.4^2))

    ∂θ(ψm)(3+im,.4)-(ψm(3+im,.4001)-ψm(3+im,.4))/.0001*(-sqrt(1-.4^2))

    ∂θ(ψconj)(3+im,.4)-(ψconj(3+im,.4001)-ψconj(3+im,.4))/.0001*(-sqrt(1-.4^2))

    ∂θ(ψmconj)(3+im,.4)-(ψmconj(3+im,.4001)-ψmconj(3+im,.4))/.0001*(-sqrt(1-.4^2))

    ∂θ(ψs)(3+im,.4)-(ψs(3+im,.4001)-ψs(3+im,.4))/.0001*(-sqrt(1-.4^2))

    ∂θ(ψsm)(3+im,.4)-(ψsm(3+im,.4001)-ψsm(3+im,.4))/.0001*(-sqrt(1-.4^2))

    ∂θ(ψsconj)(3+im,.4)-(ψsconj(3+im,.4001)-ψsconj(3+im,.4))/.0001*(-sqrt(1-.4^2))

    ∂θ(ψsmconj)(3+im,.4)-(ψsmconj(3+im,.4001)-ψsmconj(3+im,.4))/.0001*(-sqrt(1-.4^2))
end