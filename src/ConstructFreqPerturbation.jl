struct FreqPertub
    ω2::ComplexF64 
    Error::Float64
end

function ComputeDplus(ψ::QuasinormalModeFunction)
    Alm=ψ.Alm; a=ψ.a; m=ψ.m; ω=ψ.ω;
    λlm=Alm+a^2*ω^2-2*a*m*ω
    D²=λlm^2*(λlm+2)^2-8*λlm*(5*λlm+6)*(a^2*ω^2-a*m*ω)+96*λlm*a^2*ω^2
        +144*(a^2*ω^2-a*m*ω)^2
    Dplus=sqrt(D²)
    println("Dplus: ", Dplus)
    Dplus
end

function Compute𝒞plus(ψ::QuasinormalModeFunction,Dplus)
    ω=ψ.ω;
    𝒞plus=Dplus+144*ω^2
    println("𝒞plus: ", 𝒞plus)
    𝒞plus
end

function Computeγs(∂ωOplusInt,∂ωOminusInt,HplusInt,HminusInt,IplusInt,IminusInt,Dplus,m,ω)
    # denom=(2*((-1)^m * Dplus*HplusInt*∂ωOminusInt- 12*im*IplusInt*∂ωOminusInt*ω))
    # γ1num=(-(-1)^m*Dplus*IminusInt*∂ωOplusInt + 12*im*HplusInt*∂ωOminusInt*ω - 
    #     IplusInt*∂ωOminusInt*(-1)^m*Dplus + 12*im*HminusInt*∂ωOplusInt*ω+
    #     sqrt(-4*(-12*im*IminusInt*∂ωOplusInt*ω+ HminusInt*∂ωOplusInt*(-1)^m
    #     *Dplus)*((-1)^m*Dplus* HplusInt*∂ωOminusInt- 12*im*IplusInt*∂ωOminusInt
    #     *ω) + ((-1)^m * Dplus*IminusInt*∂ωOplusInt - 12*im*HplusInt*∂ωOminusInt
    #     *ω + IplusInt*∂ωOminusInt*(-1)^m*Dplus-12*im*HminusInt*∂ωOplusInt*ω)^2))
    # γ1=γ1num/denom
    # γ2num=(-(-1)^m*Dplus*IminusInt*∂ωOplusInt + 12*im*HplusInt*∂ωOminusInt*ω - 
    #     IplusInt*∂ωOminusInt*(-1)^m*Dplus + 12*im*HminusInt*∂ωOplusInt*ω-
    #     sqrt(-4*(-12*im*IminusInt*∂ωOplusInt*ω+ HminusInt*∂ωOplusInt*(-1)^m
    #     *Dplus)*((-1)^m*Dplus* HplusInt*∂ωOminusInt-12*im*IplusInt*∂ωOminusInt
    #     *ω) + ((-1)^m * Dplus*IminusInt*∂ωOplusInt - 12*im*HplusInt*∂ωOminusInt
    #     *ω + IplusInt*∂ωOminusInt*(-1)^m*Dplus-12*im*HminusInt*∂ωOplusInt*ω)^2))
    denom=(2*(Dplus*HplusInt*∂ωOminusInt- 12*im*IplusInt*∂ωOminusInt*ω))
    γ1num=(-Dplus*IminusInt*∂ωOplusInt + 12*im*HplusInt*∂ωOminusInt*ω - 
        IplusInt*∂ωOminusInt*Dplus + 12*im*HminusInt*∂ωOplusInt*ω+
        sqrt(-4*(-12*im*IminusInt*∂ωOplusInt*ω+ HminusInt*∂ωOplusInt
        *Dplus)*(Dplus* HplusInt*∂ωOminusInt- 12*im*IplusInt*∂ωOminusInt
        *ω) + ( Dplus*IminusInt*∂ωOplusInt - 12*im*HplusInt*∂ωOminusInt
        *ω + IplusInt*∂ωOminusInt*Dplus-12*im*HminusInt*∂ωOplusInt*ω)^2))
    γ1=γ1num/denom
    γ2num=(-Dplus*IminusInt*∂ωOplusInt + 12*im*HplusInt*∂ωOminusInt*ω - 
        IplusInt*∂ωOminusInt*Dplus + 12*im*HminusInt*∂ωOplusInt*ω-
        sqrt(-4*(-12*im*IminusInt*∂ωOplusInt*ω+ HminusInt*∂ωOplusInt
        *Dplus)*(Dplus* HplusInt*∂ωOminusInt-12*im*IplusInt*∂ωOminusInt
        *ω) + (Dplus*IminusInt*∂ωOplusInt - 12*im*HplusInt*∂ωOminusInt
        *ω + IplusInt*∂ωOminusInt*Dplus-12*im*HminusInt*∂ωOplusInt*ω)^2))
    γ2=γ2num/denom
    γs=(γ1,γ2)
    println("γs: ", γs)
    γs
end

function ComputeAs(𝒞plus,Dplus,ψ::QuasinormalModeFunction,γs)
    m=ψ.m;ω=ψ.ω;
    # A1=8((-1)^m*Dplus*γs[1]-12*im*ω)/𝒞plus
    # A2=8((-1)^m*Dplus*γs[2]-12*im*ω)/𝒞plus
    A1=8(Dplus*γs[1]-12*im*ω)/𝒞plus
    A2=8(Dplus*γs[2]-12*im*ω)/𝒞plus
    As=(A1,A2)
    As
end

function ComputeBs(𝒞plus,Dplus,ψ::QuasinormalModeFunction,γs)
    m=ψ.m;ω=ψ.ω;
    # B1=8((-1)^m*conj(Dplus)+12*im*conj(ω)*conj(γs[1]))/conj(𝒞plus)
    # B2=8((-1)^m*conj(Dplus)+12*im*conj(ω)*conj(γs[2]))/conj(𝒞plus)
    B1=8(conj(Dplus)+12*im*conj(ω)*conj(γs[1]))/conj(𝒞plus)
    B2=8(conj(Dplus)+12*im*conj(ω)*conj(γs[2]))/conj(𝒞plus)
    Bs=(B1,B2)
    Bs
end

function Computeω2(∂ωOplusInt,∂ωOminusInt,HplusInt,HminusInt,IplusInt,IminusInt,ψ)
    Dplus=ComputeDplus(ψ)
    # println(Dplus)
    𝒞plus= Compute𝒞plus(ψ,Dplus)
    γs=Computeγs(∂ωOplusInt,∂ωOminusInt,HplusInt,HminusInt,IplusInt,IminusInt,Dplus,ψ.m,ψ.ω)
    As=ComputeAs(𝒞plus,Dplus,ψ,γs)
    Bs= ComputeBs(𝒞plus,Dplus,ψ,γs)

    println("Bs: ",Bs)
    println("As: ",As)
    println("numerator: ",(As[1]*HplusInt+conj(Bs[1])IplusInt))
    println("∂ωOplusInt: ",∂ωOplusInt)

    ω2up=-((As[1]*HplusInt+conj(Bs[1])IplusInt)/∂ωOplusInt)
    ω2down=-((As[2]*HplusInt+conj(Bs[2])IplusInt)/∂ωOplusInt)

    ωs=(ω2up,ω2down)
    ωs
end