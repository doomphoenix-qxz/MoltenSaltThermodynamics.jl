struct MQCModel{T}
    CNumbers::Vector{T}
      
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


