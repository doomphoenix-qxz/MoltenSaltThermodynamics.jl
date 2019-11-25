module MoltenSaltThermodynamics
using LsqFit
using Statistics: mean

const R_GAS = 1.380649*6.02214076 # Boltzmann*Avogadro, values from CODATA18
const FARADAY = 96485.3365

include("extrapolation.jl")

export infinitedilution_extrapolate, plot_model  

end # module
