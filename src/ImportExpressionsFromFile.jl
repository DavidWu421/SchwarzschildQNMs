# I think the first two functions in this file are fine. Last two may need changing

function GetExprMetric(file)
    df = DataFrame(CSV.File(file,types=Dict(1=>String)))
    exprs = []
    for (d,f) ∈ eachrow(df[df.Coefficient .!= "0",[:Index,:Coefficient]])
        # println((d,f))
        # d = d == "-" ? "" : d
        # d = replace(d,"th"=>"θ")
        f = replace(f,
            "omega" => "ω",
            "rho"   => "NPvars.ρ",
            "fdagger" => "NPvars.fd",
            "beta"  => "NPvars.β",
            "tau"   => "NPvars.τ",
            "mu"    => "NPvars.μ",
            "xi"    => "NPvars.ξ"
        )
        ps = "ψ"*d
        push!(exprs,"($f)*$(ps)(r,x)")
    end
    func_def = "$(join(exprs," + "))"
    println(func_def)
    Meta.parse(func_def)
end

function MakeOp(file)
    thisexpr = GetExprMetric(file)
    Op = eval(Meta.parse("((r,x,a,m,ω,s,NPvars,ψ00,ψ01,ψ10,ψ11; pertparam=0, M=1) -> $thisexpr)"))
    # println(thisexpr)
    Op
end
