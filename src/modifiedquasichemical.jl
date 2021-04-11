"""
Holds the parameters unique to a particular binary subsystem modeled by
the MQCM; namely, the parameters in Eqns 17, 19, and 20 of Pelton et al 
paper I. 
"""
struct BinaryMQCModel{T}
    gab0::T 
  	gaparams::Vector{T}
    gbparams::Vector{T}
    # Order: Za_aa, Zb_bb, Za_ab, Zb_ba 
    coordparams::Vector{T}

end

"""
Calculate both coordination numbers Za and Zb given a particular 
Modified Quasichemical Model and its configuration in pair fractions.
Eqns. 19-20 of Pelton et al. 'The Modified Quasichemical 
Model I--Binary Solutions'
"""
function coordination(data::BinaryMQCModel, npairs)
  Za_aa, Zb_bb, Za_ab, Zb_ba = data.coordparams 
  naa, nbb, nab = npairs 
  inv_Za = 2naa/(Za_aa*(2naa + nab)) + nab/(Za_ab*(2naa + nab))
  inv_Zb = 2nbb/(Zb_bb*(2nbb + nab)) + nab/(Zb_ba*(2nbb + nab))
  return [1/inv_Za, 1/inv_Zb]
end 

"""
Calculate Δgab given a particular Modified Quasichemical Model and 
its current configuration. Equivalent to Eqn 17 in Pelton et al. I
"""
function gab(data::BinaryMQCModel, pairfracs)
	aterm = 0.0
  bterm = 0.0
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
function Sconfig(data::BinaryMQCModel, molefracs, pairfracs)
  Xa, Xb = molefracs 
  Xaa, Xbb, Xab = pairfracs 
  Za, Zb = coordination(data, pairfracs)
  totpairs = (Xa*Za+Xb*Zb)/2
  #println(totpairs)
  naa, nbb, nab = pairfracs .* totpairs 
  Ya = Xaa + Xab/2 
  Yb = Xbb + Xab/2 
  term1 = zero(eltype(molefracs))
  if Xa != zero(eltype(molefracs))
    term1 += Xa*log(Xa) 
  end
  if Xb!= zero(eltype(molefracs))
    term1 += Xb*log(Xb)
  end 
  term2 = zero(eltype(molefracs))
  if Xaa != zero(eltype(pairfracs))
    term2 += naa*log(Xaa/Ya^2) 
  end 
  if Xbb != zero(eltype(pairfracs))
    term2 += nbb*log(Xbb/Yb^2)
  end 
  if Xab != zero(eltype(pairfracs))
    term2 += nab*log(Xab/(2*Ya*Yb))
  end 
  return -R_GAS*(term1 + term2)
end 

"""
Calculates overall molar Gibbs energy of solution for a Quasichemical 
system. Eqn 9 of Pelton et al. I 
"""
function gibbs_binary_solution(data::BinaryMQCModel,molefracs,
                               pairfracs, purecompdata,temp)
  Xa, Xb = molefracs 
  if Xb == zero(eltype(molefracs))
    return gibbsenergy(purecompdata[1], temp)
  elseif Xa == zero(eltype(molefracs))
    return gibbsenergy(purecompdata[2], temp)
  end 
  Xaa, Xbb, Xab = pairfracs 
  Za, Zb = coordination(data, pairfracs)
  totpairs = (Xa*Za+Xb*Zb)/2
  #println(totpairs)
  naa, nbb, nab = pairfracs .* totpairs 
  ΔSconfig = Sconfig(data, molefracs, pairfracs)
  Δgab = gab(data, pairfracs)
  ga = gibbsenergy(purecompdata[1], temp)
  gb = gibbsenergy(purecompdata[2], temp)
  return Xa*ga + Xb*gb - temp*ΔSconfig + (nab/2)*Δgab 
end 

"""
Calculate excess Gibbs energy using Eqn 35 from Pelton and Blander.
"""
function gibbs_excess_eq35(data::BinaryMQCModel, npairs, ndiff=missing)
  if !ismissing(ndiff)
    myzero = zero(typeof(npairs[ndiff]))
  else
    myzero = zero(Float64)
  end
  aterms = myzero
  bterms = myzero
  naa, nbb, nab = npairs 
  totpairs = sum(npairs)
  Xaa, Xbb, Xab = npairs ./ totpairs  
  for i in 1:length(data.gaparams)
    aterms += data.gaparams[i]*Xaa^(i-1)
  end 
  for j in 1:length(data.gbparams)
    bterms += data.gbparams[j]*Xbb^(j-1)
  end 
  return (Xaa*Xab*aterms + Xbb*Xab*bterms) * totpairs/2 
end 

function gibbs_excess2(data::BinaryMQCModel, molfracs, pairfracs,
                      purecompdata, temp)
  lgibbs = gibbs_binary_solution(data, molfracs, pairfracs, purecompdata, temp)
  puregibbs = sum([molfracs[i] * gibbsenergy(purecompdata[i], temp) for i in 1:length(molfracs)])
  idealmixpiece = R_GAS*temp*sum([molfracs[i]*log(molfracs[i]) for i in 1:length(molfracs)])
  return lgibbs - (puregibbs + idealmixpiece)
end 

function entropy_excess(data::BinaryMQCModel, molfracs, pairfracs,
                       purecompdata, temp)
  return -FiniteDifferences.central_fdm(5,1)(T -> gibbs_excess2(data, molfracs,
                                pairfracs, purecompdata, T), temp)
end 


function enthalpy_excess(data::BinaryMQCModel, molfracs, pairfracs,
                       purecompdata, temp)
  return ForwardDiff.derivative(inv_T -> gibbs_excess2(data, molfracs,
            pairfracs, purecompdata, 1/inv_T)*inv_T, 1/temp)
  #return gibbs_excess2(data, molfracs, pairfracs, purecompdata, temp) + 
  #entropy_excess(data, molfracs, pairfracs, purecompdata, temp)*temp 
end 

function chemical_potential(data::BinaryMQCModel, purecompdata, molfracs,
                            npairs, temp, compnum=1)
  totpairs = sum(npairs)
  xpairs = npairs ./ totpairs  
  if compnum == 1 
    myfunction(Xa) = gibbs_binary_solution(data, [Xa, molfracs[2]],
                                          xpairs, purecompdata, temp)
    mycomp = molfracs[1]
  else 
    myfunction(Xb) = gibbs_binary_solution(data, [molfracs[1], Xb],
                                          xpairs, purecompdata, temp)
    mycomp = molfracs[2]
  end 
  g0 = gibbsenergy(purecompdata[compnum], temp)
  mu = ForwardDiff.derivative(myfunction, mycomp)
  gex = mu - g0
  #println("Pure Gibbs: "*string(g0)*".  Chemical potenital: " *string(mu) * ".  
  #Gibbs Ex: " * string(gex))
  return mu
end 

function gibbs_excess_comp(data::BinaryMQCModel, purecompdata, molfracs,
  npairs, temp, compnum=1)
  return chemical_potential(data, purecompdata, molfracs, npairs, temp, 
  compnum) - gibbsenergy(purecompdata[compnum], temp)
end 

# function gibbs_excess_comp2(data::BinaryMQCModel, purecompdata, molfracs,
#   npairs, temp, compnum=1)
#   pairfracs = npairs ./ sum(npairs)
#   Gextot = gibbs_excess2(data, molfracs, pairfracs, purecompdata, temp)
#   myfunc = if compnum == 1 
#     myfunca(xa) = gibbs_excess2(data, [xa, 1-xa], pairfracs, purecompdata, temp)
#     myfunca
#   elseif compnum == 2 
#     #println("Ought to get here")
#     myfuncb(xb) = gibbs_excess2(data, [1-xb, xb], pairfracs, purecompdata, temp)
#     myfuncb
#   else
#     println("SAD!")
#   end
#   #println("It's the right version")
#   mycomp = molfracs[compnum]
#   return Gextot - ForwardDiff.derivative(myfunc, mycomp) * mycomp
# end 
function activity_coeff(data::BinaryMQCModel, purecompdata, molfracs,
  npairs, temp, compnum=1)
  return exp(gibbs_excess_comp(data, purecompdata, molfracs, npairs, temp,
  compnum) / (R_GAS*temp))
end

function activity(data::BinaryMQCModel, purecompdata, molfracs,
  npairs, temp, compnum=1)
  return activity_coeff(data, purecompdata, molfracs,
  npairs, temp, compnum) * molfracs[compnum]
end

function chemical_potential2(data::BinaryMQCModel, purecompdata, molfracs,
  npairs, temp, compnum=1)
  totpairs = sum(npairs)
  xpairs = npairs ./ totpairs
  myfunction = x-> x  
  if compnum == 1 
    myfunction(naa) = gibbs_excess_eq35(data,[naa,npairs[2],npairs[3]])
  elseif compnum == 2
    myfunction(nbb) = gibbs_excess_eq35(data,[npairs[1],nbb,npairs[3]])
  else 
    throw(DomainError("compnum must be 1 or 2"))
  end 
  mycomp = purecompdata[compnum]
  mymfrac = molfracs[compnum]
  mynp = npairs[compnum]
  mypcomp = purecompdata[compnum]
  mycoorddata = data.coordparams[compnum]
  myxp = xpairs[compnum]
  myY = xpairs[compnum] + xpairs[3]/2
  #@infiltrate
  part1 = gibbsenergy(mycomp,temp) + R_GAS*temp*log(mymfrac)
  #part2_der = ForwardDiff.derivative(myfunction, mynp)
  part2_der = FiniteDifferences.central_fdm(5,1,max_range=9e-4)(myfunction, mynp)
  part2a = R_GAS*temp*log(myxp/myY^2)
  total = part1 + (mycoorddata/2)*(part2a + part2_der)
  #print("Part 1: " * string(part1) * "  ")
  #print("Part 2a: " * string(part2a) * "  ")
  #print("Derivative: "*string(part2_der)*"  ")
  #println("Total: " * string(total))
  return part1 + (mycoorddata/2)*(part2a + part2_der)
end 

function gibbs_excess_comp3(data::BinaryMQCModel, purecompdata, molfracs,
npairs, temp, compnum=1)
return chemical_potential2(data, purecompdata, molfracs, npairs, temp, 
compnum) - gibbsenergy(purecompdata[compnum], temp)
end 

function gibbs_excess_comp2(data::BinaryMQCModel, purecompdata, molfracs,
  npairs, temp, compnum=1)
  return molfracs[compnum]*(chemical_potential2(data, purecompdata, molfracs, npairs, temp, 
  compnum) - gibbsenergy(purecompdata[compnum], temp))
  end 

function activity_coeff3(data::BinaryMQCModel, purecompdata, molfracs,
npairs, temp, compnum=1)
return exp(gibbs_excess_comp2(data, purecompdata, molfracs, npairs, temp,
compnum) / (R_GAS*temp))
end

function activity_coeff2(data::BinaryMQCModel, purecompdata, molfracs,
  npairs, temp, compnum=1)
  return exp(-gibbs_excess_comp2(data, purecompdata, molfracs, npairs, temp,
  compnum) / (R_GAS*temp))
  end

function activity2(data::BinaryMQCModel, purecompdata, molfracs,
npairs, temp, compnum=1)
return activity_coeff2(data, purecompdata, molfracs,
npairs, temp, compnum) * molfracs[compnum]
end

struct PureComp{T}
  params::Vector{T}
end


function heatcapacity(pc::PureComp, T)
	A,B,a,b,c,d,e = pc.params
  return a + b*1e-3*T + c*1e5/(T^2) + d/sqrt(T) + e*1e-6*T^2
end 

function enthalpy(pc::PureComp, T)
  A,B,a,b,c,d,e = pc.params 
  T0 = 298.15 
  return A + a*T + b*1e-3*T^2 / 2 - c*1e5/T + 2d*sqrt(T) + e*1e-6*T^3 /3 - 
  (a*T0 + b*1e-3*T0^2 / 2 - c*1e5/T0 + 2d*sqrt(T0) + e*1e-6*T0^3 /3)
end 

function entropy(pc::PureComp, T)
  A,B,a,b,c,d,e = pc.params
  T0 = 298.15 
  return B + a*log(T) + b*1e-3*T - (c*1e5/2)/(T^2) - 2*d/sqrt(T) + (e*1e-6/2)*T^2 - 
  (a*log(T0) + b*1e-3*T0 - (c*1e5/2)/(T0^2) - 2*d/sqrt(T0) + (e*1e-6/2)*T0^2)
end

# TODO: There's a problem here but I'm not sure what it is. Figure out how to fix!!!
function gibbsenergy(pc::PureComp, T)
  return enthalpy(pc, T) - T*entropy(pc, T)
end 

# notation: Xa is mole fraction of a, Ya is coordiation-equivalent
# fraction of a (see Pelton and Blander's paper)
#
function find_configuration(data::BinaryMQCModel,molefracs, 
                                purecompdata,temp,
                                x0 = [0.5,0.5,0.0])
  if molefracs[1] == 1.0
    return [1.0, 0.0, 0.0]
  elseif molefracs[2] == 1.0 
    return [0.0,1.0, 0.0]
  end
  Za_aa, Zb_bb, Za_ab, Zb_ba = data.coordparams
  lower = [0.0,0.0,0.0]
  upper = [1.0, 1.0, 1.0]
  f(x) = gibbs_binary_solution(data, molefracs, x, purecompdata, temp)
  function con_c!(c, x)
    c[1] = x[1] + x[2] + x[3]
    yaa = x[1]*2/Za_aa + x[3]/Za_ab 
    ybb = x[2]*2/Zb_bb + x[3]/Zb_ba 
    c[2] = yaa /(yaa + ybb)
    c[3] = ybb /(yaa + ybb)
    c
  end 
  function con_jacobian!(J, x)
    c = Vector{Real}([0.0,0.0,0.0])
    J .= ForwardDiff.jacobian(con_c!, c, x)
    J
  end
  function con_h!(h, x, λ)

    c = Vector{Real}([0.0,0.0,0.0])
    h .= ForwardDiff.hessian(y -> sum(λ .* con_c!(c, y)), x)   
    h
  end 
  lx = [0.0,0.0,0.0]
  ux = [1.0,1.0,1.0]
  lc = [1.0, molefracs[1], molefracs[2]]
  uc = [1.0, molefracs[1], molefracs[2]]
  mygrad!(g, x) = ForwardDiff.gradient!(g, f, x)
  myhess!(h, x) = ForwardDiff.hessian!(h,f, x) 
  mydfc = TwiceDifferentiableConstraints(con_c!, con_jacobian!, con_h!,
                                        lx,ux,lc,uc)
  algo = IPNewton()
  #println("Function: ", f(x0))
  #println("Jacobian: ", ForwardDiff.jacobian(f, x0))
  #println("Hessian: ", ForwardDiff.hessian(f, x0))
  myans1 = optimize(f, mydfc, x0, algo; autodiff=:forward)
  #println(myans1)
  #println(Optim.minimizer(myans1))
  return myans1 
end

function  makeguess(Xa, Xa_maxorder, Xab_max=0.66)

    Xaa_infl = (1-Xab_max)/2 
    if Xa == Xa_maxorder 
      x0 = [Xaa_infl, Xaa_infl, Xab_max]
    else
      if Xa < Xa_maxorder
        x1aa, x2aa, y1aa, y2aa = 0, Xa_maxorder, 1, (1-Xab_max)/2
        x1bb, x2bb, y1bb, y2bb = 0, Xa_maxorder, 0, (1-Xab_max)/2
        x1ab, x2ab, y1ab, y2ab = 0, Xa_maxorder, 0, Xab_max 
      elseif Xa > Xa_maxorder 
        x1aa, x2aa, y1aa, y2aa = Xa_maxorder, 1, (1-Xab_max)/2, 0
        x1bb, x2bb, y1bb, y2bb = Xa_maxorder, 1, (1-Xab_max)/2, 1
        x1ab, x2ab, y1ab, y2ab = Xa_maxorder, 1, Xab_max, 0
      end 
      x0 = [(Xa-x1aa)*(y2aa-y1aa)/(x2aa-x1aa) + y1aa,
            (Xa-x1bb)*(y2bb-y1bb)/(x2bb-x1bb) + y1bb,
            (Xa-x1ab)*(y2ab-y1ab)/(x2ab-x1ab) + y1ab] 
    end
    return x0 
  end 
function find_configuration2(data::BinaryMQCModel,molefracs, 
                                purecompdata,temp,
                                xg = missing)
  Xa, Xb = molefracs
  Za_aa, Zb_bb, Za_ab,  Zb_ba= data.coordparams  
  function objectives!(err, pfracs)
    Xaa, Xbb, Xab =pfracs 
    yaa = Xaa*2/Za_aa + Xab/Za_ab 
    ybb = Xbb*2/Zb_bb + Xab/Zb_ba 
    err[1] = Xab^2/(Xaa*Xbb) - 4*exp(-gab(data, pfracs)/(R_GAS*temp))
    err[2] = yaa/(yaa+ybb) - Xa 
    err[3] = Xaa + Xbb+ Xab - 1.0
  end 
  if ismissing(xg) || length(xg) != 3
    # set an initital guess 
    # At max short-range ordering, assume Xab = 2/3, Xaa = Xbb = 1/6 
    # This is a significant assumtion, come back to it 
    # Else linearly interpolate. 
    Xa_maxorder = Za_ab /(Za_ab + Zb_ba)
    for Xab_max in [0.66, 0.75, 0.50, 0.25, 0.90, 0.10, 0.01, 0.99]
      x0 = makeguess(Xa, Xa_maxorder, Xab_max)
      answer = nlsolve(objectives!, x0, autodiff=:forward)
      xaag, xbbg, xabg = answer.zero
      if 0.0 <= xaag <= 1.0 && 0.0 <= xbbg <= 1.0 && 0.0 <= xabg <= 1.0 
        return answer 
      end  
    end   
    println("Warning, no good solution was found")
  else 
    x0 = xg 
  end 
  answer = nlsolve(objectives!, x0, autodiff=:forward)
  #println(answer)
  answer 
end 
function find_configuration3(data::BinaryMQCModel,molefracs, 
                                purecompdata,temp,
                                xg = missing)
  Xa, Xb = molefracs
  Za_aa, Zb_bb, Za_ab,  Zb_ba= data.coordparams 
  function solvethis(nab)
    naa = (Xa - nab/Za_ab)*Za_aa/2 
    nbb = (Xb - nab/Zb_ba)*Zb_bb/2 
    ptot = naa+nbb+nab 
    pairfracs = [naa, nbb, nab] ./ ptot 
    Xaa, Xbb, Xab = pairfracs
    return Xab^2/(Xaa*Xbb) - 4*exp(-gab(data, pairfracs)/(R_GAS*temp))
  end 
  if ismissing(xg) || length(xg) != 3
    # set an initital guess 
    x0 = (0.0, max(Za_aa/2, Zb_bb/2))
  else 
    x0 = xg 
  end 
  xl, xu = x0 
  answer = optimize(solvethis, xl, xu)
  #println(answer)
  nab_final = Optim.minimizer(answer)
  naa_final = (Xa - nab_final/Za_ab)*Za_aa/2 
  nbb_final = (Xb - nab_final/Zb_ba)*Zb_bb/2 
  npairs = [naa_final, nbb_final, nab_final]
  xpairs = [i/sum(npairs) for i in npairs]
  return npairs, xpairs 
end

function bruteforce_configuration(data::BinaryMQCModel,molefracs, 
                                purecompdata,temp,
                                x0 = [0.5,0.5,0.0])
  XaaXbbspace = Iterators.product(0.00:0.01:1.00, 0.00:0.01:1.00, 0.00:0.01:1.00) |> collect
  paramspace = filter(x -> x[1] + x[2] + x[3] <= 1.0, XaaXbbspace) 
end
