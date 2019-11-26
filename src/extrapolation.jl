const R = R_GAS 
const F = FARADAY 

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
    return x -> tofit(x, fit.param[1])
end 

function infinitedilution_extrapolate(xdata, vdata, n, T)
    logxs = log.(xdata)
    fitmodel = gen_model(n, T)
    initialmodel = fitamodel(fitmodel, logxs, vdata) 
    Rsq = get_rsquared(initialmodel, logxs, vdata)
    finalmodel = initialmodel 
    modeltwopoints = fitamodel(fitmodel, logxs[1:2], vdata[1:2])
    oldRsq = get_rsquared(modeltwopoints, logxs[1:2], vdata[1:2])
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
    scene = scatter(logxs, edata, label=:Experimental, legend=:topleft)
    plot!(scene, mynewxs, mynewys, label=:Extrapolated)
    title!(scene, "Standard Potential Extrapolation")
    xlabel!(scene, "Natural log of Concetration of $addname in $mixname")
    ylabel!(scene, "Measured Potential (V)")
end 

        
