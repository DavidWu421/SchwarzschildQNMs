# This defines a structure where the first entry is some filename for the file corresponding to operator, and the second
# entry is the filename of the file with the derivative of that operator. The third and fourth entries are the operator
# and the derivative of that operator respectively
struct OperatorShift
    filename::AbstractString
    Op::Function
end

# This function takes in two file names. Then it reads this file using Make... defined in ImportExpressionsFromFile.jl
# to construct the actual operators from the csv file. Then it makes an OperatorShift construct with the names and
# operators
function OperatorShift(file::AbstractString)
    Op = MakeOp(file)
    OperatorShift(file,Op)
end
