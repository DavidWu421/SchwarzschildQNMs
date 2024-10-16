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
    ψ02 = ∂θ(ψ01);
    ψ03 = ∂θ(ψ02);
    ψ04 = ∂θ(ψ03);
    ψ05 = ∂θ(ψ04);
    ψ06 = ∂θ(ψ05);

    ψ10 = ∂r(ψ00); 
    ψ11 = ∂r(ψ01);
    ψ12 = ∂r(ψ02);
    ψ13 = ∂r(ψ03);
    ψ14 = ∂r(ψ04);
    ψ15 = ∂r(ψ05);
    ψ16 = ∂r(ψ06);

    ψ20 = ∂r(ψ10); 
    ψ21 = ∂r(ψ11);
    ψ22 = ∂r(ψ12);
    ψ23 = ∂r(ψ13);
    ψ24 = ∂r(ψ14);
    ψ25 = ∂r(ψ15);
    ψ26 = ∂r(ψ16);

    ψ30 = ∂r(ψ20); 
    ψ31 = ∂r(ψ21);
    ψ32 = ∂r(ψ22);
    ψ33 = ∂r(ψ23);
    ψ34 = ∂r(ψ24);
    ψ35 = ∂r(ψ25);
    ψ36 = ∂r(ψ26);

    ψ40 = ∂r(ψ30); 
    ψ41 = ∂r(ψ31);
    ψ42 = ∂r(ψ32);
    ψ43 = ∂r(ψ33);
    ψ44 = ∂r(ψ34);
    ψ45 = ∂r(ψ35);
    ψ46 = ∂r(ψ36);

    ψ50 = ∂r(ψ40); 
    ψ51 = ∂r(ψ41);
    ψ52 = ∂r(ψ42);
    ψ53 = ∂r(ψ43);
    ψ54 = ∂r(ψ44);
    ψ55 = ∂r(ψ45);
    ψ56 = ∂r(ψ46);

    ψ60 = ∂r(ψ50); 
    ψ61 = ∂r(ψ51);
    ψ62 = ∂r(ψ52);
    ψ63 = ∂r(ψ53);
    ψ64 = ∂r(ψ54);
    ψ65 = ∂r(ψ55);
    ψ66 = ∂r(ψ56);
    
    
    FF = let a = ψ.a, m = ψ.m ,ω = ψ.ω ,s = ψ.s, ψ=ψ, ψL=ψL, Ops=Ops,ψ01=ψ01,ψ02=ψ02,ψ03=ψ03,ψ04=ψ04,ψ05=ψ05,ψ06=ψ06,ψ10=ψ10,ψ11=ψ11,ψ12=ψ12,ψ13=ψ13,ψ14=ψ14,ψ15=ψ15,ψ16=ψ16,ψ20=ψ20,ψ21=ψ21,ψ22=ψ22,ψ23=ψ23,ψ24=ψ24,ψ25=ψ25,ψ26=ψ26,ψ30=ψ30,ψ31=ψ31,ψ32=ψ32,ψ33=ψ33,ψ34=ψ34,ψ35=ψ35,ψ36=ψ36,ψ40=ψ40,ψ41=ψ41,ψ42=ψ42,ψ43=ψ43,ψ44=ψ44,ψ45=ψ45,ψ46=ψ46,ψ50=ψ50,ψ51=ψ51,ψ52=ψ52,ψ53=ψ53,ψ54=ψ54,ψ55=ψ55,ψ56=ψ56,ψ60=ψ60,ψ61=ψ61,ψ62=ψ62,ψ63=ψ63,ψ64=ψ64,ψ65=ψ65,ψ66=ψ66, weight=weight
    function F(r,z;isconjugate=false,isminus=false,LHSisminus=false)
        ψL(r,z,isminus=LHSisminus)*Ops(r,z,a,m,ω,s,ψ,ψ01,ψ02,ψ03,ψ04,ψ05,ψ06,ψ10,ψ11,ψ12,ψ13,ψ14,ψ15,ψ16,ψ20,ψ21,ψ22,ψ23,ψ24,ψ25,ψ26,ψ30,ψ31,ψ32,ψ33,ψ34,ψ35,ψ36,ψ40,ψ41,ψ42,ψ43,ψ44,ψ45,ψ46,ψ50,ψ51,ψ52,ψ53,ψ54,ψ55,ψ56,ψ60,ψ61,ψ62,ψ63,ψ64,ψ65,ψ66;isconjugate=isconjugate,isminus=isminus)*weight(r,z)
    end
    end
    
    OperatorSandwich(ψ,ψL,FF,weight)
end
