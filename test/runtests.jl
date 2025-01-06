# push!(LOAD_PATH, "C:/Users/dwuuu/Documents/GitHub/ContourIntegrals.jl")
# push!(LOAD_PATH, "C:/Users/dwuuu/Documents/GitHub/SchwarzschildQNMs")
# push!(LOAD_PATH, "C:/Users/dwuuu/Documents/GitHub/KerrQuasinormalModes.jl")

using ContourIntegrals
using KerrQNMShifts
using KerrQuasinormalModes
using Test
using Zygote

println("Done usings")

using Plots


function make_plots(func, zbound,zstep; filename="savedfile",pertparam=.1)
    # Define the complex plane grid
    xmin, xmax = 2.15, 2.25
    ymin, ymax = -.3, .3
    n_points = 100
    x = range(xmin, xmax, length=n_points)
    y = range(ymin, ymax, length=n_points)
    grid = [x[i] + im * y[j] for i in 1:n_points, j in 1:n_points]
    # Get the global color limits for consistent scaling
    clim = (-3, 3)
    zstepper=-zbound
    while zstepper < zbound
        # Compute function values
        values = [func(r,zstepper,pertparam=pertparam) for r in grid]
        # Exclude the region 1.99 < x < 2.01 and -0.01 < y < 0.01
        for i in 1:n_points, j in 1:n_points
            if 1.9 < real(grid[i, j]) < 2.1 && -0.1 < imag(grid[i, j]) < 0.1
                values[i, j] = NaN
            end
        end

        # Extract real and imaginary parts
        real_part = real.(values)
        imag_part = imag.(values)

        # Plot the real part
        plot_real = heatmap(x, y, real_part',
            title = "Real Part of operator(r)",
            xlabel = "Re(r)",
            ylabel = "Im(r)",
            clim=clim,
            colorbar_title = "Re(operator)")

        # Plot the imaginary part
        plot_imag = heatmap(x, y, imag_part',
            title = "Imaginary Part of operator(r)",
            xlabel = "Re(r)",
            ylabel = "Im(r)",
            clim=clim,
            colorbar_title = "Im(operator)")

        # Display the plots
        plot(plot_real, plot_imag, layout = (1, 2), size = (800, 400))
        if zstepper < 0
            savefig(filename*"m"*string(abs(round(zstepper,digits=1)))*".png")
        elseif zstepper >=0
            savefig(filename*string(abs(round(zstepper,digits=1)))*".png")
        end

        zstepper+=zstep
    end
end

function evaltime(ψ,start_r, num_steps,step_size)
    value=ψ(start_r,.5, pertparam=.1)
    for i in 0:num_steps
        value+=ψ(start_r+i*step_size,.5,pertparam=.1)
    end
    return value
end

function finitedifference(ψ,r,θ,dr;isconjugate=false)
    return (ψ(r+dr,θ,isconjugate=isconjugate)-ψ(r,θ,isconjugate=isconjugate))/dr
end

function compute_derivative_matrix(ψ, r, θ)
    # Initialize a 5x5 matrix to hold the derivatives
    derivatives = Matrix{ComplexF64}(undef, 7, 7)
    for i in 0:6
        for j in 0:6
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

println("Done deriv funcs")

# @testset "KerrAsPertub" begin
    
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

    # testradial = LineSegment(r₊ + Δr+1im,r₊ + Δr+.5im,true)
    # TestC = testradial ⊗ angular
    
    # @time Integrate(HplusSchw, TestC,pertparam=.1)[1]
    # @time Integrate(HminusSchw, TestC,pertparam=.1)[1]

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

    @show ∂ωOplusSchw(3,0.4,pertparam=.1)
    @show ∂ωOminusSchw(3,0.4,pertparam=.1)
    @show HplusSchw(3,0.4,pertparam=.1)
    @show HminusSchw(3,0.4,pertparam=.1)
    @show IplusSchw(3,0.4,pertparam=.1)
    @show IminusSchw(3,0.4,pertparam=.1)

    @show ∂ωOplusSchw(8+im,0.6,pertparam=.1)
    @show ∂ωOminusSchw(8+im,0.6,pertparam=.1)
    @show HplusSchw(8+im,0.6,pertparam=.1)
    @show HminusSchw(8+im,0.6,pertparam=.1)
    @show IplusSchw(8+im,0.6,pertparam=.1)
    @show IminusSchw(8+im,0.6,pertparam=.1)

    println("Complied Operators")

    ∂ω𝒪plusSchw= Integrate(∂ωOplusSchw, TheContour,pertparam=.1,abstol=1e-6)[1]
    ∂ω𝒪minusSchw= conj(Integrate(∂ωOminusSchw, TheContour,pertparam=.1,abstol=1e-6)[1])
    ℋplusSchw= Integrate(HplusSchw, TheContour,pertparam=.1,abstol=1e-6)[1]
    ℋminusSchw= conj(Integrate(HminusSchw, TheContour,pertparam=.1,abstol=1e-6)[1])
    ℐplusSchw= Integrate(IplusSchw, TheContour,pertparam=.1,abstol=1e-6)[1]
    ℐminusSchw= conj(Integrate(IminusSchw, TheContour,pertparam=.1,abstol=1e-6)[1])

    @show ∂ω𝒪plusSchw
    @show ∂ω𝒪minusSchw
    @show ℋplusSchw
    @show ℋminusSchw
    @show ℐplusSchw
    @show ℐminusSchw 
    
    ω2s=Computeω2(∂ω𝒪plusSchw,∂ω𝒪minusSchw,ℋplusSchw,ℋminusSchw,ℐplusSchw,ℐminusSchw,ψ)
    @show ω2s

# end

# @testset "TeukolskyOperator Hermiticity" begin
#     ψref = qnmfunctionnew(-2,2,2,0,0.0)
#     ψrefm = qnmfunctionnew(-2,2,2,0,0.0,modesign="minus")

#     ψ = qnmfunctionnew(-2,4,3,0,0.0)
#     ψ.ω=ψref.ω
#     ψ.a=ψref.a
#     ψ.m=ψref.m
#     ψ.s=ψref.s

#     ψm = qnmfunctionnew(-2,4,3,0,0.0,modesign="minus")
#     ψm.ω=ψrefm.ω
#     ψm.a=ψrefm.a
#     ψm.m=ψrefm.m
#     ψm.s=ψrefm.s

#     ψtest = qnmfunctionnew(2,5,3,0,0.0)
#     ψtest.ω=ψref.ω
#     ψtest.a=ψref.a
#     ψtest.m=ψref.m
#     ψtest.s=ψref.s

#     ψtestm = qnmfunctionnew(2,5,3,0,0.0,modesign="minus")
#     ψtestm.ω=ψrefm.ω
#     ψtestm.a=ψrefm.a
#     ψtestm.m=ψrefm.m
#     ψtestm.s=ψrefm.s

#     Σ = let a= ψref.a
#         (r,z) -> Complex(r)^2+a^2*z^2
#     end
#     Δ = let a= ψref.a
#         (r,z) -> Complex(r)^2+a^2-2*r
#     end
#     ζ = let a= ψref.a
#         (r,z) -> r-im*a*z
#     end

#     weightplus = let a= ψref.a
#         (r,z) ->8*ζ(r,z)^(4)*Σ(r,z)/((Δ(r,z))^2)
#     end
#     weightminus = let a= ψref.a
#         (r,z) ->8*ζ(conj(r),z)^(4)*Σ(conj(r),z)/((Δ(conj(r),z))^2)
#     end

#      ## Define the useful contours
#      r₊ = ψref.R.r₊ ; r₋ = ψref.R.r₋ ; s = ψref.s ; Δr = 0.1*(r₊-r₋); ϵ = eps(0.1);

#      point1 = r₊ + Δr - Δr*im
#      point2 = r₊ - Δr - Δr*im
 
#      radial1 = SemiInfiniteLine(point1 , point1 + Δr*im , false)
#      angular = LineSegment(-1.0+100*ϵ , 1.0-100*ϵ , true) #to avoid the NaNs at the edges
#      C1 = radial1 ⊗ angular
 
#      radial2 = LineSegment(point1,point2,true)
#      C2 = radial2 ⊗ angular
 
#      radial3 = SemiInfiniteLine(point2 , point2 + Δr*im , true)
#      C3 = radial3 ⊗ angular
 
#      TheContour = C1⊕C2⊕C3

#     Oplusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/Opluscoefficients.csv"
#     Ominusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/Ominuscoefficients.csv"
#     Oplus = OperatorShift(Oplusfile)
#     Ominus = OperatorShift(Ominusfile)

#     OplusHermiticity1 = OperatorSandwich(ψ,Oplus,weightplus,ψtest).Op
#     OplusHermiticity2 = OperatorSandwich(ψtest,Oplus,weightplus,ψ).Op
#     OminusHermiticity1 = OperatorSandwich(ψm,Ominus,weightminus,ψtestm).Op
#     OminusHermiticity2= OperatorSandwich(ψtestm,Ominus,weightminus,ψm).Op

#     @show OplusHermiticity1(3+im,0.5)
#     @show OplusHermiticity2(3+im,0.5) 
#     @show OminusHermiticity1(3+im,0.5)
#     @show OminusHermiticity2(3+im,0.5)

#     @show 𝒪plusHermiticity1= Integrate(OplusHermiticity1, TheContour)[1]
#     @show 𝒪plusHermiticity2= Integrate(OplusHermiticity2, TheContour)[1]
#     @show 𝒪minusHermiticity1= Integrate(OplusHermiticity1, TheContour)[1]
#     @show 𝒪minusHermiticity2= Integrate(OplusHermiticity2, TheContour)[1]

# end

# @testset "TeukolskyOperator" begin
    
    # ψ = qnmfunctionnew(-2,2,2,0,0.1)
    # ψ0 = qnmfunctionnew(-2,2,2,0,0.)

    # ψm = qnmfunctionnew(-2,2,2,0,0.1,modesign="minus")
    # ψ0m = qnmfunctionnew(-2,2,2,0,0.,modesign="minus")

    # ω2=ψ.ω-ψ0.ω

    # ψtest = qnmfunctionnew(2,5,3,0,0.)
    # ψtest.ω=ψ0.ω
    # ψtest.a=ψ0.a
    # ψtest.m=ψ0.m
    # ψtest.s=ψ0.s

    # # Compile ψ
    # ψ(1,.5)
    # println("Past ψ compile")

    # Σ = let a= ψ.a
    #     (r,z) -> Complex(r)^2+a^2*z^2
    # end
    # Δ = let a= ψ.a
    #     (r,z) -> Complex(r)^2+a^2-2*r
    # end
    # ζ = let a= ψ.a
    #     (r,z) -> r-im*a*z
    # end

    # weightplus = let a= ψ.a
    #     (r,z) ->8*ζ(r,z)^(4)*Σ(r,z)/((Δ(r,z))^2)
    # end

    # weightminus = let a= ψ.a
    #     (r,z) ->8*ζ(conj(r),z)^(4)*Σ(conj(r),z)/((Δ(conj(r),z))^2)
    # end


    # Σ0 = let a= ψ0.a
    #     (r,z) -> Complex(r)^2+a^2*z^2
    # end
    # Δ0 = let a= ψ0.a
    #     (r,z) -> Complex(r)^2+a^2-2*r
    # end
    # ζ0 = let a= ψ0.a
    #     (r,z) -> r-im*a*z
    # end
    # weight0plus = let a= ψ0.a
    #     (r,z) ->8*ζ0(r,z)^(4)*Σ0(r,z)/((Δ0(r,z))^2)
    # end
    # weight0minus = let a= ψ0.a
    #     (r,z) ->8*ζ0(conj(r),z)^(4)*Σ0(conj(r),z)/((Δ0(conj(r),z))^2)
    # end
    # println("Past weight")

    # ## Define the useful contours
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

    # r₊0 = ψ0.R.r₊ ; r₋0 = ψ0.R.r₋ ; s0 = ψ0.s ; Δr0 = 0.1*(r₊0-r₋0); ϵ = eps(0.1);

    # point10 = r₊0 + Δr0 - Δr0*im
    # point20 = r₊0 - Δr0 - Δr0*im

    # radial10 = SemiInfiniteLine(point10 , point10 + Δr0*im , false)
    # angular0 = LineSegment(-1.0+100*ϵ , 1.0-100*ϵ , true) #to avoid the NaNs at the edges
    # C10 = radial10 ⊗ angular0

    # radial20 = LineSegment(point10,point20,true)
    # C20 = radial20 ⊗ angular0

    # radial30 = SemiInfiniteLine(point20 , point20 + Δr0*im , true)
    # C30 = radial30 ⊗ angular0

    # TheContour0 = C10⊕C20⊕C30

    # println("Done Contours")

    # KerrOplusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/KerrOpluscoefficients.csv"
    # Oplusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/Opluscoefficients.csv"
    # dwOplusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/dwOpluscoefficients.csv"
    # KerrOminusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/KerrOminuscoefficients.csv"
    # Ominusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/Ominuscoefficients.csv"
    # dwOminusfile = "C:/Users/dwuuu/Documents/UT Academics/Research/Ringdown/Mathematica/SavedFiles/dwOminuscoefficients.csv"

    # KerrOplus = OperatorShift(KerrOplusfile)
    # KerrOminus = OperatorShift(KerrOminusfile)
    # Oplus = OperatorShift(Oplusfile)
    # ∂ωOplus = OperatorShift(dwOplusfile)
    # Ominus = OperatorShift(Ominusfile)
    # ∂ωOminus = OperatorShift(dwOminusfile)

    # println("Made operator shifts")

    # KerrOplusKerr = OperatorSandwich(ψ,KerrOplus,weightplus,ψ).Op
    # KerrOminusKerr = OperatorSandwich(ψm,KerrOminus,weightminus,ψm).Op
    # # HermiticityTest1 = OperatorSandwich(ψ0,Oplus,weight0,ψtest).Op
    # # HermiticityTest2 = OperatorSandwich(ψtest,Oplus,weight0,ψ0).Op

    # println("Made the first set of operators")
    
    # OplusKerr = OperatorSandwich(ψ0,Oplus,weight0plus,ψ).Op
    # OplusSchw = OperatorSandwich(ψ0,Oplus,weight0plus,ψ0).Op
    # ∂ωOplusSchw = OperatorSandwich(ψ0,∂ωOplus,weight0plus,ψ0).Op
    # OminusKerr = OperatorSandwich(ψ0m,Ominus,weight0minus,ψm).Op
    # OminusSchw = OperatorSandwich(ψ0m,Ominus,weight0minus,ψ0m).Op
    # ∂ωOminusSchw = OperatorSandwich(ψ0m,∂ωOminus,weight0minus,ψ0m).Op

    # println("Made Operators")

    # # @show HermiticityTest1(3,0.5)
    # # @show HermiticityTest2(3,0.5)

    # @show KerrOplusKerr(3+im,0.5)
    # @show KerrOminusKerr(3+im,0.5)

    # @show OplusKerr(3+im,0.5)
    # @show OplusSchw(3+im,0.5)
    # @show ∂ωOplusSchw(3+im,0.5)

    # @show OminusKerr(3+im,0.5)
    # @show OminusSchw(3+im,0.5)
    # @show ∂ωOminusSchw(3+im,0.5)

    # println("Complied Operators")

    # # HermiticityTestResult1 = Integrate(HermiticityTest1,TheContour0)[1]
    # # @show HermiticityTestResult1
    # # HermiticityTestResult2 = Integrate(HermiticityTest2,TheContour0)[1]
    # # @show HermiticityTestResult2
   
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
