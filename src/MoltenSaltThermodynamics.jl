module MoltenSaltThermodynamics
using LsqFit
using Statistics: mean
using Plots
using LinearAlgebra
using LaTeXStrings
using Optim 
using ForwardDiff 
using Debugger
using NLsolve

const R_GAS = 1.380649*6.02214076 # Boltzmann*Avogadro, values from CODATA18
const FARADAY = 96485.3365

include("extrapolation.jl")
include("modifiedquasichemical.jl")
include("numericalmethods.jl")
include("test_mqcm.jl")

export infinitedilution_extrapolate, plot_model, compare_actcoeffs, mytest1

end # module
