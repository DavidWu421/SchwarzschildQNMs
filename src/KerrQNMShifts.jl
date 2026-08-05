__precompile__()
module KerrQNMShifts

using KerrQuasinormalModes
using PyCall
using CSV
using DataFrames
using Plots
using Cubature
using ContourIntegrals
# Write your package code here.

include("KerrQuasinormalModeAdditions.jl")
include("ImportExpressionsFromFile.jl")
include("OperatorShiftType.jl")
include("OperatorSandwichType.jl")
include("ConstructFreqPerturbation.jl")
include("NPScalars.jl")


export qnm, qnmfunctionnew, importqnm
export OperatorShift, OperatorSandwich, Computeω2, ComputeDplus,Compute𝒞plus, ComputeDplusm, ComputeAs, ComputeBs, CheckComputeω2Error

end
