# function NP_quantities(r, θ, a)
#     # Assumes M=1
#     Δ = r^2 - 2*r + a^2

#     ρ = -1/(r-im*a*cos(θ))
#     fd = -1/(ρ*(r+im*a*cos(θ)))
#     β = -fd*ρ*cot(θ)/(2*sqrt(2))
#     τ = -im*a*ρ^2*fd*sin(θ)/sqrt(2)
#     μ = ρ^3*fd*Δ/2
#     ξ = ρ^2*fd*(r-1)/2

#     return (; ρ, fd, β, τ, μ, ξ)
# end

function rho(r,z,a)
    -1/(r-im*a*z)
end

function fdagger(r,z,a,ρ)
    -1/(ρ*(r+im*a*z))
end

function beta(z,ρ,fd)
    -fd*ρ*(z/sqrt(1-z^2))/(2*sqrt(2))
end

function tau(z,a,ρ,fd)
    -im*a*ρ^2*fd*sqrt(1-z^2)/sqrt(2)
end

function mu(r,a,ρ,fd)
    Δ = r^2 - 2*r + a^2
    ρ^3*fd*Δ/2
end

function xi(r,z,a,ρ,fd)
    ρ^2*fd*(r-1)/2
end