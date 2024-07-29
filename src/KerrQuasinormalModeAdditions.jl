# I think this file should be good too. This just calls in Leo Stein's QNM code.

const qnm = PyNULL()

function __init__()
    copy!(qnm, pyimport("qnm"))
    qnm.download_data()
end

function qnmfunctionnew(s,l,m,n,a,plusminus::AbstractString; qnm=qnm)
    # qnmfunction defined on line 202 of KerrQuasinormalModes.jl/src/ModeFunctionInterface/Interface.jl
    # m should be the m of the plus mode, regardless of if you want the plus or minus mode
    grav_freq = qnm.modes_cache(s=s,l=l,m=m,n=n)
    ω, Alm, Cllʼm = grav_freq(a=a)
    if plusminus == "plus"
        qnmfunction(Custom; s=s,l=l,m=m,n=n,a=a,ω=ω,Alm=Alm,Cllʼ=Cllʼm)
    elseif plusminus == "minus"
        grav_freq_minus = qnm.modes_cache(s=s,l=l,m=-m,n=n)
        ω, Almm, Cllʼmm = grav_freq_minus(a=a)
        qnmfunction(Custom; s=s,l=l,m=m,n=n,a=a,ω=-conj(ω),Alm=Almm,Cllʼ=Cllʼmm)
    else
        println("plusminus argument is not valid. Must be either 'plus' or 'minus'")
    end
end
