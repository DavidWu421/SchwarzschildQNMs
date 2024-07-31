# I think the first two functions in this file are fine. Last two may need changing

function GetExprMetric(file)
    df = DataFrame(CSV.File(file,types=Dict(1=>String)))
    exprs = []
    for (d,f) ∈ eachrow(df[df.Coefficient .!= "0",[:Index,:Coefficient]])
        # println((d,f))
        d = d == "-" ? "" : d
        d = replace(d,"th"=>"θ")
        f = replace(f,"omega"=>"ω")
        f = replace(f,"spina"=>"a")
        ps = "ψ"*d
        push!(exprs,"($f)*$(ps)(r,x;isconjugate=isconjugate)")
    end
    func_def = "$(join(exprs," + "))"
    # println(func_def)
    Meta.parse(func_def)
end

function MakeOp(file)
    thisexpr = GetExprMetric(file)
    Op = eval(Meta.parse("((r,x,a,m,ω,s,ψ00,ψ01,ψ02,ψ03,ψ04,ψ05,ψ06,ψ10,ψ11,ψ12,ψ13,ψ14,ψ15,ψ16,ψ20,ψ21,ψ22,ψ23,ψ24,ψ25,ψ26,ψ30,ψ31,ψ32,ψ33,ψ34,ψ35,ψ36,ψ40,ψ41,ψ42,ψ43,ψ44,ψ45,ψ46,ψ50,ψ51,ψ52,ψ53,ψ54,ψ55,ψ56,ψ60,ψ61,ψ62,ψ63,ψ64,ψ65,ψ66; isconjugate=false, M=1) -> $thisexpr)"))
    println(thisexpr)
    Op
end
