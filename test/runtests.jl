using MoltenSaltThermodynamics
using Test
using DataFrames, CSV

mydata1 = CSV.File("cocl2_naclkcl.csv") |> DataFrame 
n = 2.0
T = 1073.0
edata = mydata1[!, Symbol("800C")].* -1  
xdata = mydata1[!, Symbol("CoCl2_NaClKCl")]
model = infinitedilution_extrapolate(xdata, edata, n, T)

@testset "MoltenSaltThermodynamics.jl" begin
  @test isapprox(model.standardpotential, -1.1259; atol=1e-4)
  @test isapprox(model.slope, 0.04623; atol=1e-5)
end
