# I think this file should be good too. This just calls in Leo Stein's QNM code.

const qnm = PyNULL()

function __init__()
    copy!(qnm, pyimport("qnm"))
    qnm.download_data()
end

function qnmfunctionnew(s,l,m,n,a; qnm=qnm,is_conjugate=false,is_minus=false)
    # qnmfunction defined on line 202 of KerrQuasinormalModes.jl/src/ModeFunctionInterface/Interface.jl
    # m should be the m of the plus mode, regardless of if you want the plus or minus mode
    grav_freq = qnm.modes_cache(s=s,l=l,m=m,n=n)
    ω, Alm, Cllʼ = grav_freq(a=a)
    if is_minus==true
        m=-m
        ω=-conj(ω)
        Alm=conj(Alm)
        Cllʼ=[i % 2 == 0 ? -x : x for (i, x) in enumerate(Cllʼ)]
        Cllʼ=(-1)^l*conj.(Cllʼ)
        # The following if statement might be buggy for s!=-2,0. It's needed because of how the
        # Cllʼ are indexed in the SpinWeightedSpheroidal of the KerrQuasinormalModes package. In particular
        # how the max(abs(s),abs(m)) condition triggers.
        if (abs(m)<abs(s)) & (isodd(abs(m)))
            Cllʼ=-Cllʼ
        end
    end
    if is_conjugate==true
        Alm=conj(Alm)
    end
    qnmfunction(Custom; s=s,l=l,m=m,n=n,a=a,ω=ω,Alm=Alm,Cllʼ=Cllʼ,is_conjugate=is_conjugate,is_minus=is_minus)
end
