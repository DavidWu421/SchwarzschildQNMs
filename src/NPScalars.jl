function NP_vars(r,z,a)

    ρ  = -1/(r-im*a*z)
    fd = -1/(ρ*(r+im*a*z))
    β = -fd*ρ*(z/sqrt(1-z^2))/(2*sqrt(2))
    τ = -im*a*ρ^2*fd*sqrt(1-z^2)/sqrt(2)
    μ =ρ^3*fd*( r^2 - 2*r + a^2)/2
    ξ = ρ^2*fd*(r-1)/2

    return (; ρ, fd, β, τ, μ, ξ)
end