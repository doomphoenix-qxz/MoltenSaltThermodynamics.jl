const R = R_GAS 
const F = FARADAY 

struct SimpleLinearModel{T<:Real}
    standardpotential::T
    slope::T
    rsquared::T
    SimpleLinearModel{T}(a,b,c) where {T<:Real} = new(a,b,c)
end 

SimpleLinearModel(a::Float64, b::Float64, c::Float64) = SimpleLinearModel{Float64}(a,b,c)
#SimpleLinearModel(a::Real, b::Real)  = SimpleLinearModel{Real}(promote(a, b)...)
function (s::SimpleLinearModel)(x)
    return s.slope .* x .+ s.standardpotential
end 

struct ModelResult{S, T}
    themodel::S
    theplot::T
end 

function gen_model(n, T=723.0)
    coeff = R*T/(n*F)
    model(x, E0) = E0[1] .+ coeff .* x
    return model
end 

function get_rsquared(model, xdata, ydata)
    ȳ = mean(ydata)
    sstot = sum((ydata .- ȳ).^2)
    ssres = sum((ydata .- model(xdata)).^2)
    R² = 1 - ssres/sstot 
    return R²
end 

function fitamodel(tofit, xdata, ydata)
    p0 = [1.0]
    fit = curve_fit(tofit, xdata, ydata, p0)
    e0 = fit.param[1]
    slope = tofit(1.0, e0)  - tofit(0.0, e0)
    Rsq = get_rsquared(x->tofit(x,e0), xdata, ydata)
    return SimpleLinearModel(e0, slope, Rsq)
end 

function infinitedilution_extrapolate(xdata, vdata, n, T)
    logxs = log.(xdata)
    fitmodel = gen_model(n, T)
    initialmodel = fitamodel(fitmodel, logxs, vdata) 
    Rsq = initialmodel.rsquared
    finalmodel = initialmodel 
    modeltwopoints = fitamodel(fitmodel, logxs[1:2], vdata[1:2])
    oldRsq = modeltwopoints.rsquared 
    # this is to initialize finalmodel with the right type; we're going
    # to figure out what the REAL final model is in a bit
    for i in 3:length(xdata)
        if Rsq > 0.99
            break # assume if Rsq is already high, model is already good
        elseif i == 3
            if oldRsq > Rsq finalmodel = modeltwopoints end  
        end 
        xn = logxs[1:i]
        vn = vdata[1:i]
        newmodel = fitamodel(fitmodel, xn, vn)
        Rsq = get_rsquared(newmodel, xn, vn)
        if Rsq < oldRsq 
            if Rsq - oldRsq > 0.02 || Rsq < 0.95 || Rsq /oldRsq < 0.95
                break 
            end 
        end 
        oldRsq = Rsq 
        finalmodel = newmodel 
    end 
    return finalmodel 
end 

function plot_model(model, xdata, vdata, mixname="Mixture", addname="Additive")
    logxs = log.(xdata)
    myxrange = minimum(logxs) .- 0.1:0.01:0.0
    mynewxs = collect(myxrange)
    mynewys = broadcast(model, mynewxs)
    scene = scatter(logxs, vdata, label=:Experimental, legend=:topleft)
    plot!(scene, mynewxs, mynewys, label=:Extrapolated)
    title!(scene, "Standard Potential Extrapolation")
    xlabel!(scene, "Natural log of Concetration of $addname in $mixname")
    ylabel!(scene, "Measured Potential (V)")
    return ModelResult(model, scene) 
end 

function compare_actcoeffs(extrap_model, purecomp_enot, xdata, edata,mixname="Mixture", addname="Additive")
    logxs = log.(xdata)
    rightslope = extrap_model(1.0) - extrap_model(0.0)
    puremodel(x) = x .* rightslope .+ purecomp_enot 
    unrefined_pure = edata .- puremodel(logxs)
    unrefined_infd = edata .- extrap_model(logxs)
    newcoeff = 1.0/rightslope 
    pureγ = exp.(unrefined_pure .* newcoeff)
    infdγ = exp.(unrefined_infd .* newcoeff)
    print("Pure Component Activities for $addname in $mixname: ", pureγ)
    print("Infinite Dilution Activities for $addname in $mixname: ", infdγ)
    scene = scatter(log10.(xdata), log10.(infdγ), label=:Extrapolated, legend=:topleft)
    scatter!(scene, log10.(xdata), log10.(pureγ), label=Symbol("Pure Component"))
    title!(scene, latexstring("\\mathrm{Activity \\: Coefficient \\: Comparison \\: for \\: $addname \\: in \\: $mixname}"))
    xlabel!(scene, latexstring("\\mathrm{Log_{10} \\: of \\: Concentration}"))
    ylabel!(scene, latexstring("\\mathrm{Log}_{10} \\mathrm{ \\: of \\: Activity \\: Coefficient}"))
    return scene
end 
