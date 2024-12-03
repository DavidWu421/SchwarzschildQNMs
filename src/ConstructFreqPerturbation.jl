struct FreqPertub
    ω2::ComplexF64 
    Error::Float64
end

function ComputeDplus(ψ::QuasinormalModeFunction)
    Alm=ψ.Alm; a=ψ.a; m=ψ.m; ω=ψ.ω;
    λlm=Alm+a^2*ω^2-2*a*m*ω
    D²=λlm^2*(λlm+2)^2-8*λlm*(5*λlm+6)*(a^2*ω^2-a*m*ω)+96*λlm*a^2*ω^2+
        144*(a^2*ω^2-a*m*ω)^2
    Dplus=sqrt(D²)
    Dplus
end

function Compute𝒞plus(ψ::QuasinormalModeFunction,Dplus)
    ω=ψ.ω;
    𝒞plus=Dplus^2+144*ω^2
    println("𝒞plus: ", 𝒞plus)
    𝒞plus
end

function Computeγs(∂ωOplusInt,∂ωOminusInt,HplusInt,HminusInt,IplusInt,IminusInt,Dplus,m,ω)
    a=Dplus*HplusInt*∂ωOminusInt-12*im*ω*IplusInt*∂ωOminusInt
    b=Dplus*(IplusInt*∂ωOminusInt+IminusInt*∂ωOplusInt)-12*im*ω*(HplusInt*∂ωOminusInt+HminusInt*∂ωOplusInt)
    c=Dplus*HminusInt*∂ωOplusInt-12*im*ω*IminusInt*∂ωOplusInt
    γ1=(-b+sqrt(b^2-4*a*c))/(2*a)
    γ2=(-b-sqrt(b^2-4*a*c))/(2*a)
    γs=(γ1,γ2)
    println("γs: ", γs)
    γs
end


function ComputeAs(𝒞plus,Dplus,ψ::QuasinormalModeFunction,γs)
    ω=ψ.ω;
    A1=16(Dplus*γs[1]-12*im*ω)/𝒞plus
    A2=16(Dplus*γs[2]-12*im*ω)/𝒞plus
    As=(A1,A2)
    As
end

function ComputeBs(𝒞plus,Dplus,ψ::QuasinormalModeFunction,γs)
    ω=ψ.ω;
    B1=16(conj(Dplus)+12*im*conj(ω)*conj(γs[1]))/conj(𝒞plus)
    B2=16(conj(Dplus)+12*im*conj(ω)*conj(γs[2]))/conj(𝒞plus)
    Bs=(B1,B2)
    Bs
end

function Computeω2(∂ωOplusInt,∂ωOminusInt,HplusInt,HminusInt,IplusInt,IminusInt,ψ)
    Dplus=ComputeDplus(ψ)
    𝒞plus= Compute𝒞plus(ψ,Dplus)
    γs=Computeγs(∂ωOplusInt,∂ωOminusInt,HplusInt,HminusInt,IplusInt,IminusInt,Dplus,ψ.m,ψ.ω)
    As=ComputeAs(𝒞plus,Dplus,ψ,γs)
    Bs= ComputeBs(𝒞plus,Dplus,ψ,γs)

    println("As: ",As)
    println("Bs: ",Bs)
    println("numerator: ",(As[1]*HplusInt+conj(Bs[1])IplusInt))
    println("∂ωOplusInt: ",∂ωOplusInt)

    ω2up=-((As[1]*HplusInt+conj(Bs[1])IplusInt)/∂ωOplusInt)
    ω2down=-((As[2]*HplusInt+conj(Bs[2])IplusInt)/∂ωOplusInt)

    ωs=(ω2up,ω2down)
    ωs
end