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

function take_derivative(f,r)
    if isreal(r)
        return gradient(x -> real(f(x)), r)[1]+im*gradient(x -> imag(f(x)), r)[1]
    elseif imag(r)!=0
        return gradient(x -> real(f(x)), r)[1] |> conj
    end
end

function teukop(f,r)
    (12 * im .* r + 2 - 4 * 2 .* 0.995004 .* (2  - r) ./ (1 - 0.995004 .^ 2) + 1^ 2 .* r .^ 3 -4 * im * 1 .* r .^ 2 - r + 0.995004 .^ 2 .* (6 - 3 * r) ./ (1 - 0.995004 .^ 2) + (2 - r) .* (2 .^ 2 + 1) ./ (1 - 0.995004 .^ 2)) ./ (2 * r .^ 6 .* (2  - r))*f(r)+(-1 + r) ./ r .^ 6*take_derivative(f,r)-(-2 + r) ./ (2 * r .^ 5)*take_derivative(rₚ -> take_derivative(f, rₚ), r)
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

weight = let s = -2 , a=0, z= 0.995004
    (r) ->2^(5+s)*ζ(r,z)^(8+2*s)*Σ(r,z)*sqrt(1-z^2)/((Δ(r,z))^2)
end


# function teukop(f,r,x)
#     (12 * im  .* omega .* r + 2 - 4 * m .* x .* (2  - r) ./ (1 - x .^ 2) + omega .^ 2 .* r .^ 3 -
#     4 * im * omega .* r .^ 2 - r + x .^ 2 .* (6 - 3 * r) ./ (1 - x .^ 2) + (2 * M - r) .* 
#     (m .^ 2 + 1) ./ (1 - x .^ 2)) ./ (2 * r .^ 6 .* (2  - r))*f(r)+
#     (-1 + r) ./ r .^ 6*ForwardDiff.derivative(f,r)+
#     -(-2 + r) ./ (2 * r .^ 5)*ForwardDiff.derivative(rₚ -> ForwardDiff.derivative(f, rₚ), r)
# end



# @testset "Integration" begin

    function ψ(r; isconjugate=false)
        Complex(im*(r-1))^(1+im)*exp(im*Complex(r)^2)
    end
    function ϕ(r; isconjugate=false)
        Complex(im*(r-1))^(3+2*im)*exp(2*im*Complex(r)^2)
    end

    r₊ = 1; Δr = 0.01; ϵ = eps(0.1);

    point1 = r₊ + Δr - Δr*im
    point2 = r₊ - Δr - Δr*im

    radial1 = SemiInfiniteLine(point1 , point1 + Δr*im , false)
    C1 = radial1

    radial2 = LineSegment(point1,point2,true)
    C2 = radial2

    radial3 = SemiInfiniteLine(point2 , point2 + Δr*im , true)
    C3 = radial3

    TestContour = C1⊕C2⊕C3

    AngularContour = LineSegment(100*ϵ , 1.0-100*ϵ , true)

    function integrand1(r; isconjugate=false)
        ϕ(r)*weight(r)*teukop(ψ,r)
    end

    function integrand2(r; isconjugate=false)
        ψ(r)*weight(r)*teukop(ϕ,r)
    end

    function angularintegrand(z; isconjugate=false)
        sqrt(1-z^2)
    end

    @show Integrate(integrand1,TestContour)[1]
    @show Integrate(integrand2,TestContour)[1]

#end