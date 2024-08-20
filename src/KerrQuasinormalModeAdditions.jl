# I think this file should be good too. This just calls in Leo Stein's QNM code.

const qnm = PyNULL()

function __init__()
    copy!(qnm, pyimport("qnm"))
    qnm.download_data()
end

function qnmfunctionnew(s,l,m,n,a; qnm=qnm,modesign="plus")
    # qnmfunction defined on line 202 of KerrQuasinormalModes.jl/src/ModeFunctionInterface/Interface.jl
    # m should be the m of the plus mode, regardless of if you want the plus or minus mode
    grav_freq = qnm.modes_cache(s=s,l=l,m=m,n=n)
    ω, Alm, Cllʼ = grav_freq(a=a)
    if modesign=="minus"
        m=-m
        ω=-conj(ω)
        Alm=conj(Alm)
        Cllʼ=[i % 2 == 0 ? -x : x for (i, x) in enumerate(Cllʼ)]
        Cllʼ=(-1)^l*conj.(Cllʼ)
    end
    qnmfunction(Custom; s=s,l=l,m=m,n=n,a=a,ω=ω,Alm=Alm,Cllʼ=Cllʼ)
end
