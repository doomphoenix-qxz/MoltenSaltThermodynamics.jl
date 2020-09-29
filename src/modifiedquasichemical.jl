
"""
Holds the parameters unique to a particular binary subsystem modeled by
the MQCM; namely, the parameters in Eqns 17, 19, and 20 of Pelton et al 
paper I. 
"""
struct BinaryMQCModel{T}
    gab0::T 
  	gaparams::Vector{T}
    gbparams::Vector{T}
    coordparams::Vector{T}
end

"""
Calculate both coordination numbers Za and Zb given a particular 
Modified Quasichemical Model and its configuration in pair fractions.
Eqns. 19-20 of Pelton et al. 'The Modified Quasichemical 
Model I--Binary Solutions'
"""
function coordination(data::BinaryMQCModel, pairfracs)
  Za_aa, Zb_bb, Za_ab, Zb_ba = data.coordparams 
  Xaa, Xbb, Xab = pairfracs 
  inv_Za = 2Xaa/(Za_aa*(2Xaa + Xab)) + Xab/(Za_ab*(2Xaa + Xab))
  inv_Zb = 2Xbb/(Zb_bb*(2Xbb + Xab)) + Xab/(Zb_ba*(2Xbb + Xab))
  return [1/inv_Za, 1/inv_Zb]
end 

"""
Calculate Δgab given a particular Modified Quasichemical Model and 
its current configuration. Equivalent to Eqn 17 in Pelton et al. I
"""
function gab(data::BinaryMQCModel, pairfracs)
	aterm = 0
  Xaa, Xbb, Xab = pairfracs 
  for i in 1:length(data.gaparams)
    aterm += data.gaparams[i] * Xaa^i 
  end 
  for j in 1:length(data.gbparams)
    bterm += data.gbparams[j] * Xbb^j
  end 
  return data.gab0 + aterm + bterm 
end 

"""
Calculates configurational entropy of a Modified Quasichemical System
given mole factions and configuration. Eqn 10 of Pelton et al I.  
"""
function Sconfig(molefracs, pairfracs)
  Xa, Xb = molefracs 
  Xaa, Xbb, Xab = pairfracs 
  Ya = Xaa + Xab/2 
  Yb = Xbb + Xab/2 
  term1 = 0
  if Xa != 0 
    term1 += Xa*log(Xa) 
  end
  if Xb!= 0
    term1 += Xb*log(Xb)
  end 
  term2 = 0 
  if Xaa != 0
    term2 += Xaa*log(Xaa/Ya^2) 
  end 
  if Xbb != 0
    term2 += Xbb*log(Xbb/Yb^2)
  end 
  if Xab != 0 
    term2 += Xab*log(Xab/(2*Ya*Yb))
  end 
  return -R_GAS*(term1 + term2)
end 

"""
Calculates overall molar Gibbs energy of solution for a Quasichemical 
system. Eqn 9 of Pelton et al. I 
"""
function gibbs_binary_solution(data::BinaryMQCModel,molefracs, 
                               pairfracs, purecompdata,temp)
  if Xb == 0
    return gibbsenergy(purecompdata[1], temp)
  elseif Xa ==0
    return gibbsenergy(purecompdata[2], temp)
  end 
  Xa, Xb = molefracs 
  Xaa, Xbb, Xab = pairfracs 
  #Za, Zb = coordination(data, pairfracs)
  ΔSconfig = Sconfig(molefracs, pairfracs)
  Δgab = gab(data, pairfracs)
  ga = gibbsenergy(purecompdata[1], temp)
  gb = gibbsenergy(purecompdata[2], temp)
  return Xa*ga + Xb*gb - temp*ΔSconfig + (Xab/2)*Δgab 
end 


struct PureComp{T}
  params::Vector{T}
end



function heatcapacity(pc::PureComp, T)
	A,B,a,b,c,d,e = pc.params
  return a + b*T + c/(T^2) + d/sqrt(T) + e*T^2
end 

function enthalpy(pc::PureComp, T)
  A,B,a,b,c,d,e = pc.params
  return A + a*T + b*T^2 / 2 - c/T + 2d*sqrt(T) + e*T^3 /3
end 

function entropy(pc::PureComp, T)
  A,B,a,b,c,d,e = pc.params
  return B + a*log(T) + b*T - (c/2)/(T^2) - 2*d/sqrt(T) + (e/2)*T^2
end

function gibbsenergy(pc::PureComp, T)
  return enthalpy(pc, T) - T*entropy(pc, T)
end 

# notation: Xa is mole fraction of a, Ya is coordiation-equivalent
# fraction of a (see Pelton and Blander's paper)
#
function find_configuration(data::BinaryMQCModel,molefracs, 
                                purecompdata,temp,
                               x0 = [0.5,0.5,0.5], algo = LGBFS())
  f(x) = gibbsbinarysolution(data, molefracs, x, purecompdata, temp)
  myans1 = optimize(f, x0, algo; autodiff=:forward)
  println(myans)
end 
