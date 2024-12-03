push!(LOAD_PATH, "C:/Users/dwuuu/Documents/GitHub/ContourIntegrals.jl")
push!(LOAD_PATH, "C:/Users/dwuuu/Documents/GitHub/SchwarzschildQNMs")
push!(LOAD_PATH, "C:/Users/dwuuu/Documents/GitHub/KerrQuasinormalModes.jl")

using ContourIntegrals
using KerrQNMShifts
using KerrQuasinormalModes
using Test
using Zygote

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

function take_derivative(f,r,z,wrt)
    if wrt=="r"
        if isreal(r)
            return gradient(x -> real(f(x,z)), r)[1]+im*gradient(x -> imag(f(x,z)), r)[1]
        elseif imag(r)!=0
            return gradient(x -> real(f(x,z)), r)[1] |> conj
        end
    end
    if wrt=="z"
        return gradient(x -> real(f(r,x)), z)[1]+im*gradient(x -> imag(f(r,x)), z)[1]
    end
end

function take_second_derivative(f,r,z,wrt)
    take_derivative(f,r,z,wrt)
    if wrt=="r"
        if isreal(r)
            return gradient(x -> real(take_derivative(f,x,z,wrt)), r)[1]+im*gradient(x -> imag(take_derivative(f,x,z,wrt)), r)[1]
        elseif imag(r)!=0
            return gradient(x -> real(take_derivative(f,x,z,wrt)), r)[1] |> conj
        end
    end
    if wrt=="z"
        return gradient(x -> real(take_derivative(f,r,x,wrt)), z)[1]+im*gradient(x -> imag(take_derivative(f,r,x,wrt)), z)[1]
    end
end

function teukop(f,r,x)
    let omega=1,M=1,m=2
        (12 * im * M .* omega .* r + 2 * M - 4 * m .* x .* (2 * M - r) ./ (1 - x .^ 2) + omega .^ 2 .* r .^ 3 - 4 * im * omega .* r .^ 2 - r + x .^ 2 .* (6 * M - 3 * r) ./ (1 - x .^ 2) + (2 * M - r) .* (m .^ 2 + 1) ./ (1 - x .^ 2)) ./ (2 * r .^ 6 .* (2 * M - r))*f(r,x)-x ./ (2 * r .^ 6 .* sqrt(1 - x .^ 2))*(-take_derivative(f,r,x,"z")*sqrt(1-x^2))-1 ./ (2 * r .^ 6)*(take_second_derivative(f,r,x,"z")*(1-x^2)-take_derivative(f,r,x,"z")*x)+(-M + r) ./ r .^ 6*take_derivative(f,r,x,"r")-(-2 * M + r) ./ (2 * r .^ 5)*take_second_derivative(f,r,x,"r")
    end
end

Σ = let a=0
    (r,z) -> Complex(r)^2+a^2*z^2
end
Δ = let a=0
    (r,z) -> Complex(r)^2+a^2-2*r
end
ζ = let a=0
    (r,z) -> r-im*a*z
end

weight = let s = -2 , a=0
    (r,z) ->8*ζ(r,z)^(4)*Σ(r,z)/((Δ(r,z))^2)
end


# @testset "Integration" begin

    function ψ(r,z; isconjugate=false)
        Complex(im*(r))^(1/2)*exp(im*Complex(r))*sqrt(1-z^2)
    end
    function ϕ(r,z; isconjugate=false)
        Complex(im*(r))^(3/2)*exp(im*Complex(r))*(1-z^2)
    end

    r₊ = 0; Δr = 0.01; ϵ = eps(0.1);

    point1 = r₊ + Δr - Δr*im
    point2 = r₊ - Δr - Δr*im

    AngularContour = LineSegment(-1+100*ϵ , 1.0-100*ϵ , true)

    radial1 = SemiInfiniteLine(point1 , point1 + Δr*im , false)
    C1 = radial1 ⊗ AngularContour

    radial2 = LineSegment(point1,point2,true)
    C2 = radial2 ⊗ AngularContour

    radial3 = SemiInfiniteLine(point2 , point2 + Δr*im , true)
    C3 = radial3 ⊗ AngularContour

    TestContour = C1⊕C2⊕C3

    AngularContour = LineSegment(100*ϵ , 1.0-100*ϵ , true)

    function integrand1(r,z; isconjugate=false)
        ϕ(r,z)*weight(r,z)*teukop(ψ,r,z)
    end

    function integrand2(r,z; isconjugate=false)
        ψ(r,z)*weight(r,z)*teukop(ϕ,r,z)
    end

    function testintegrand(r,z; isconjugate=false)
        r*z
    end

    TestrContour1 = LineSegment(-1 , 2.0, true)
    TestrContour2 = LineSegment(2.0 , 2.0+im, true)
    TestzContour = LineSegment(0 , 3.0, true)
    TestrzContour1= TestrContour1 ⊗ TestzContour
    TestrzContour2= TestrContour2 ⊗ TestzContour
    TestrzContour=TestrzContour1⊕TestrzContour2



    @show Integrate(integrand1,TestContour)[1]
    @show Integrate(integrand2,TestContour)[1]

#end