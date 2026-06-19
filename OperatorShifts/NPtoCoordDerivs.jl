using CSV
using DataFrames

df = CSV.read("operators.csv", DataFrame)

# Convert "01" -> (0,1), etc.
function parse_index(s)
    (parse(Int, s[1]), parse(Int, s[2]))
end

# Build dictionary
coeffs = Dict(
    parse_index(row.Index) => row.Coefficient
    for row in eachrow(df)
)

println(coeffs)

coeffs[(0,0)] = coeffs[(0,0)] +
    coeffs[(1,0)]* (-((im*(-a + (a^2 + r^2)*omega))/(a^2 + r*(-2*M + r))))+ 
    coeffs[(0,1)]*(im*(1 + a*(-1 + x^2)*omega))/((-im*r + a*x)*sqrt(2 - 2*x^2)) +
    coeffs[(1,1)]*((1 + a*(-1 + x^2)*omega)*(a^3*x*omega - 
    a^2*(1 + x + im*r*omega) + a*r*(im + r*x*omega) + 
    r*(2*M + r*(-1 - im*r*omega))))/((a^2 + r*(-2*M + r))*(-im*r + 
    a*x)^2*sqrt(2 - 2*x^2))

coeffs[(1,0)] = coeffs[(1,0)] + coeffs[(1,1)]*(im*(1 + a*(-1 + x^2)*omega))/((-im*r + a*x)*sqrt(2 - 2*x^2))

coeffs[(0,1)] = coeffs[(0,1)]*1/(sqrt(2)*(r + im*a*x)) + coeffs[(1,1)]* (-a^3*x*omega + a^2*(1 + x + im*r*omega) + 
    r*(-2*M + r + im*r^2*omega) - 
    a*r*(im + r*x*omega))/(sqrt(2)*(a^2 + r*(-2*M + r))*(-im*r + a*x)^2)