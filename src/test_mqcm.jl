function mytest1()
x_data_a2 = [0.1, 0.102, 0.153, 0.2, 0.202, 
    0.248, 0.297, 0.3, 0.348, 0.399, 0.4, 
    0.5, 0.5, 0.6, 0.601, 0.699, 0.7, 0.794, 
    0.8, 0.9]

a2_data_a2 = [0.9, 0.866488436, 0.76805011, 
    0.63, 0.649506857, 0.524300393, 0.38397441, 
    0.347, 0.256229362, 0.162509557, 0.149, 
    0.073868064, 0.0665, 0.0281, 0.030561125, 
    0.014647489, 0.0116, 0.006033884, 
    0.0042, 0.0012]

T_data_a2 = [998., 1073., 1073., 998., 1073., 
    1073., 1073., 998., 1073., 1073., 998., 
    1073., 998., 998., 1073., 1073., 998., 
    1073., 998., 998.]

x_data_Tf = Float64[0, 0, 0, 0.049533151, 0.097097104, 0.1,
    0.149600304, 0.192848845, 0.2, 0.200826782, 0.25, 
    0.250912166, 0.280680928, 0.295130904, 0.3, 0.304846515, 
    0.306117172, 0.31, 0.312350413, 0.31546046, 0.32, 
    0.324720543, 0.33, 0.331772771, 0.333, 0.34, 0.346155446, 
    0.35, 0.350583994, 0.354739648, 0.369980613, 0.37, 
    0.375348261, 0.391228458, 0.395210589, 0.4, 0.418663535, 
    0.42, 0.432446439, 0.446446984, 0.45, 0.451241823, 
    0.465687834, 0.47, 0.480656179, 0.485668569, 0.49637962, 
    0.5, 0.5, 0.501330506, 0.516401957, 0.52, 0.526433927, 
    0.530199064, 0.535751614, 0.546510668, 0.55, 0.55092091, 
    0.555777151, 0.57, 0.577064846, 0.6, 0.611741422, 0.631046965, 
    0.648471584, 0.65, 0.68, 0.7, 0.751358143, 0.8, 
    0.848845459, 0.9, 0.947644142, 1, 1, 1]

Tf_data_Tf = Float64[1041, 1043.45, 1044.626691, 1023.27028, 
    997.9162418, 983, 962.5594888, 913.7211707, 908, 
    891.9362264, 833.5, 803.3174683, 755.19996, 738.1305908, 
    716, 700.78878, 707.6685787, 700, 703.9170541, 
    701.2203147, 699, 704.0910172, 701, 706.0509132, 
    705.35, 702, 703.3818541, 713, 705.6608689, 
    701.9512391, 711.9000928, 710, 713.7865827, 722.8688513, 
    729.2058882, 740, 742.475943, 746, 748.6008856, 
    750.7296288, 760, 750.7372196, 755.7544464, 761, 
    759.1474261, 757.9121455, 758.2744076, 761.5, 
    757.95, 761.5921327, 755.7726235, 761.5, 752.1919898, 
    757.0571252, 748.6111, 743.5970352, 761.5, 751.9226969, 
    745.7508941, 752, 761.5291335, 748, 751.051037, 
    771.8169362, 792.59026, 747, 746, 868, 882.1841289, 
    908.5, 932.9572053, 956.5, 974.5672105, 985, 987.05, 987.00111]

γ1_data_γ1 = [0.0048, 0.0042, 0.0046, 3.90E-03, 0.0056, 
    0.009, 0.011, 0.029, 0.057, 0.073, 0.09, 1.13E-01, 
    0.13, 3.00E-01, 0.29, 0.45, 6.14E-01, 0.66, 0.79, 
    0.7875]

x_data_γ1 = [0.015, 0.02, 0.088, 0.1, 0.168, 0.2, 
    0.252, 0.3, 0.325, 0.343, 0.398, 0.4, 
    0.469, 0.5, 0.563, 0.6, 0.7, 0.71, 0.792, 
    0.8]

T_data_γ1 = [1073., 1073., 1073., 1073., 1073., 
    1073., 1073., 1073., 1073., 1073., 1073., 
    1073., 1073., 1073., 1073., 1073., 1073., 
    1073., 1073., 1073.]

x_data_Hex = [0.0608, 0.1059, 0.2047, 
    0.2602, 0.3524, 0.4495, 0.53, 
    0.623, 0.6577, 0.8147, 0.8816, 
    0.9415]

Hex_data_Hex = [-3313.728, -5614.928, 
    -10447.448, -12991.32, -15627.24, 
    -15627.24, -15008.008, -13167.048, 
    -12811.408, -8393.104, -5640.032, 
    -3016.664]

T_data_Hex = [1073., 1073., 1073., 1073., 
    1073., 1073., 1073., 1073., 1073., 
    1073., 1073., 1073.]

  # Pure Comp data from Table I of Pelton and Chartrand, Modified Quasichemical I
  KCls = MoltenSaltThermodynamics.PureComp([-436684.1, 82.55023, 40.01578, 25.46801, 3.64845, 0.,0.])
  KCll = MoltenSaltThermodynamics.PureComp([-421824.9, 86.5525, 73.59656, 0.,0.,0.,0.])
  MgCl2s = MoltenSaltThermodynamics.PureComp([-641616.0, 89.629, 54.58434, 21.42127, -11.12119, 399.1767, -2.356672])
                     MgCl2l = MoltenSaltThermodynamics.PureComp([-606887.4, 117.29708, 92.048, 0.,0.,0.,0.])
                     K2MgCl4s = MoltenSaltThermodynamics.PureComp([-1550013.0, 216.8, 263.05,0.,0.,0.,0.])
                     KMgCl3s = MoltenSaltThermodynamics.PureComp([-1100924.0, 162.3,157.65,0.,0.,0.,0.0])
  Temp = 1073.0 
  mygab0 = -17947.0 
  mygai = [-1026.0]
  mygbj = [-14801.0]
  mycoordparams = [6.0,6.0,3.0,6.0]
  KClMgCl2_model = BinaryMQCModel(mygab0, mygai, mygbj, mycoordparams) 
  initguess = [0.990, 0.005, 0.005]
  my_mfracs = [[x, 1-x] for x in 0.999:-0.001:0.001] 
  myfunction(mf, pfg=missing) = find_configuration2(KClMgCl2_model, mf,[KCll, MgCl2l], Temp, pfg)
  manyans = Vector{Any}(undef, length(my_mfracs))
  manyans[1] = myfunction(my_mfracs[1], initguess)
  max_nab = mycoordparams[1]/2  
  for i in 2:length(manyans)
    manyans[i] = myfunction(my_mfracs[i], manyans[i-1].zero)
  end 
  coordnums = [coordination(KClMgCl2_model, manyans[i].zero) for i in 1:length(my_mfracs)]
  totpairs = [sum(coordnums[i] .* my_mfracs[i])/2 for i in 1:length(my_mfracs)]
  Xminims = [x.zero for x in manyans]
  npairspace = Xminims .* totpairs 
  xaa = [thing[1] for thing in Xminims]
  xbb = [thing[2] for thing in Xminims]
  xab = [thing[3] for thing in Xminims]
  naa = [thing[1] for thing in npairspace]
  nbb = [thing[2] for thing in npairspace]
  nab = [thing[3] for thing in npairspace]
  Xaspace = [my_mfracs[i][1] for i in 1:length(my_mfracs)]
  # Yes, I know that I don't have to use a different internal variable
  # for the following functions. I just want to. 
  mygibbs(i) = gibbs_binary_solution(KClMgCl2_model, my_mfracs[i], Xminims[i], [KCll, MgCl2l], Temp)
  myent(m) = Sconfig(KClMgCl2_model, my_mfracs[m], Xminims[m])
  mygibbsex(j) = gibbs_excess_eq35(KClMgCl2_model, [naa[j], nbb[j], nab[j]])
  mygibbsex2(k) = gibbs_excess2(KClMgCl2_model, my_mfracs[k], Xminims[k], [KCll, MgCl2l], Temp)
  mygibbsexc1(k) = gibbs_excess_comp2(KClMgCl2_model,[KCll, MgCl2l], my_mfracs[k], npairspace[k],  Temp)
  mygibbsexc2(k) = gibbs_excess_comp2(KClMgCl2_model,[KCll, MgCl2l], my_mfracs[k], npairspace[k], Temp,2)
  myenthex2(l) = enthalpy_excess(KClMgCl2_model, my_mfracs[l], Xminims[l], [KCll, MgCl2l], Temp)
  myact(m) = activity(KClMgCl2_model, [KCll, MgCl2l], my_mfracs[m], npairspace[m], Temp)
  gibbsltot = [mygibbs(i) for i in 1:length(manyans)]
  entconfig = [myent(m) for m in 1:length(manyans)]
  gibbsl = [mygibbsex(j) for j in 1:length(manyans)]
  gibbsx2 = [mygibbsex2(k) for k in 1:length(manyans)]
  gibbsxc1 = [mygibbsexc1(k) for k in 1:length(manyans)]
  gibbsxc2 = [mygibbsexc2(k) for k in 1:length(manyans)]
  enthex = [myenthex2(l) for l in 1:length(manyans)]
  activkcl = [myact(m) for m in 1:length(manyans)]
  #println(activkcl)
  #println(length(gibbsl))
  sc1 = plot(xab, label="Xab")
  plot!(sc1, xaa, label="Xaa")
  plot!(sc1, xbb, label="Xbb")
  plot!(sc1, naa, label="naa")
  plot!(sc1, nab, label="nab")
  plot!(sc1, nbb, label="nbb")
  savefig(sc1, "/home/doomphoenix/ENTER/Research/MQCM/testfig1.png")
  display(sc1)
  #println(my_mfracs)
  #print(Xaspace)
  sc2 = plot(Xaspace ,gibbsl, label="GEx by Eqn 35")
  plot!(sc2, Xaspace, gibbsx2, label="GEx by definition")
  plot!(sc2, Xaspace, enthex, label="Hex")
  scatter!(sc2, x_data_Hex, Hex_data_Hex, label="Hex data")
  savefig(sc2, "/home/doomphoenix/ENTER/Research/MQCM/testfig2.png")
  ##display(sc2) 
  #println("Gibbs of solid KCl at 1400: ", gibbsenergy(KCls, 1400.0))
  #println("Gibbs of solid MgCl2 at 1400: ", gibbsenergy(MgCl2s, 1400.0))
  sc3 = plot(Xaspace, gibbsltot, label="Total Gibbs liquid")
  savefig(sc3, "/home/doomphoenix/ENTER/Research/MQCM/testfig3.png")
  ##display(sc3)
  sc4 = plot(Xaspace, entconfig, label="Config. entropy")
  savefig(sc4, "/home/doomphoenix/ENTER/Research/MQCM/testfig4.png")
  ##display(sc4)
  sc5 = plot(Xaspace ,activkcl, label="Act. KCl calc")
  plot!(sc5, [0.0,1.0],[1.0,0.0])
  scatter!(sc5, x_data_a2, a2_data_a2, label="Activ. data KCl")
  savefig(sc5, "/home/doomphoenix/ENTER/Research/MQCM/testfig5.png")
  sc6 = plot(Xaspace ,gibbsxc1, label="GEx KCl")
  plot!(sc6, Xaspace, gibbsxc2, label="GEx MgCl2")
  savefig(sc6, "/home/doomphoenix/ENTER/Research/MQCM/testfig6.png")
  println("Gibbs of total liquid at mole frac ", Xaspace[500], " KCl: ", gibbsltot[500])
  println("Gibbs excess of total liquid at mole frac ", Xaspace[500], " KCl: ", gibbsx2[500])
  println("Entropy of total liquid at mole frac ", Xaspace[500], " KCl: ", gibbsltot[500])
  println("Config Entropy of total liquid at mole frac ", Xaspace[500], " KCl: ", entconfig[500])
end
