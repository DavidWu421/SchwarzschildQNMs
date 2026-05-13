# abstract type Callable end
# abstract type CallableAtom <: Callable end
# abstract type CallableCombination <: Callable end

# ### Radial Functions
# struct Testfunc{T} <: CallableAtom
#     is_conjugate::Bool;is_minus::Bool;s::Int64; m::Int64; a::Float64; ω::Complex{Float64}
# end

struct OperatorSandwich{T1,T2}
    ψ::T1
    ψL::T2
    Op::Function
    weight::Function
end

function OperatorSandwich(ψL,OpShift::OperatorShift,weight,ψ)
    Ops = OpShift.Op

        ψ00 = ψ;
        ψ01 = ∂θ(ψ00);
        
        ψ10 = ∂r(ψ00); 
        ψ11 = ∂r(ψ01);
        
        # Determine ω and m based on is_conjugate
        ω = ψ.is_conjugate ? -conj(ψ.ω) : ψ.ω
        m = ψ.is_conjugate ? -ψ.m : ψ.m
        
        FF = let a = ψ.a, m = m ,ω = ω ,s = ψ.s, ψ=ψ, ψL=ψL, Ops=Ops,ψ01=ψ01,ψ10=ψ10,ψ11=ψ11, weight=weight
        function F(r,z; pertparam=0)
            NPvars = NP_vars(r,z,a)
            ψL(r,z)*Ops(r,z,a,m,ω,s,NPvars,ψ,ψ01,ψ10,ψ11; pertparam=pertparam)*weight(r,z)
        end
    end
    
    OperatorSandwich(ψ,ψL,FF,weight)
end
