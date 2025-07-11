using CSV
using DataFrames
using Turing
using DifferentialEquations
using StatsPlots
using Random

Random.seed!(544);


# Define logistic model.
function logistic_growth(dx, x, p, t)
    # Model parameters.
    r, K = p

    # Evaluate differential equations.
    dx[1] = r*x[1]*( 1 - x[1]/K)

    return nothing
end


@model function fit_logistic_all(obsdata)
    # Prior distributions.
    σ ~ InverseGamma(2, 2) #Cauchy(0, 2)#
    
    # Random effects for strain, temperature, batch and their combinations on growth rate.
    sigma_r0 ~ InverseGamma(2, 2) #Cauchy(0, 2)#
    sigma_rt ~ InverseGamma(2, 2) #Cauchy(0, 2)#
    sigma_rs ~ InverseGamma(2, 2) #Cauchy(0, 2)#
    sigma_rst ~ InverseGamma(2, 2) #Cauchy(0, 2)#
    sigma_rb ~ InverseGamma(2, 2) #Cauchy(0, 2)#

    # Random effects for strain, temperature, batch and their combinations on carrying capacity.
    sigma_kt ~ InverseGamma(2, 2) #Cauchy(0, 2)#
    sigma_ks ~ InverseGamma(2, 2) #Cauchy(0, 2)#
    sigma_kst ~ InverseGamma(2, 2) #Cauchy(0, 2)#
    sigma_kb ~ InverseGamma(2, 2) #Cauchy(0, 2)#

    r_0 ~ Uniform(0, 3)
    K_0 ~ Uniform(0, 2)#LogNormal(log(150), 0.1)
    
    b_rt = Dict()
    b_rs = Dict()
    b_rst = Dict()
    b_rb = Dict()
    b_kt = Dict()
    b_ks = Dict()
    b_kst = Dict()
    # b_kb = Dict()
    
    for strain in obsdata[4]
        b_rs[string(strain)] ~ Normal(0, sigma_rs^2)
        b_ks[string(strain)] ~ Normal(0, sigma_ks^2)

        for temp in obsdata[3]
            b_rst[string(strain, "_", temp)] ~ Normal(0, sigma_rst^2)
            b_kst[string(strain, "_", temp)] ~ Normal(0, sigma_kst^2)
            if !haskey(b_rt, temp)
                b_rt[temp] ~ Normal(0, sigma_rt^2)
            end
            if !haskey(b_kt, temp)
                b_kt[temp] ~ Normal(0, sigma_kt^2)
            end
        end 
    end
    
    for batch in obsdata[2]
        b_rb[batch] ~ Normal(0, sigma_rb^2)
        # b_kb[batch] ~ Normal(0, sigma_kb^2)
    end

    # i is for the combination of batch, strain and temperature
    for i in eachindex(obsdata[1])
        x0i = [obsdata[1][i][1][1]] # initial condition for the i-th experiment
        tspan = extrema(obsdata[1][i][2])

        # r is the growth rate for the combination of batch, strain and temperature
        r = r_0 + b_rs[obsdata[1][i][5]] + b_rt[obsdata[1][i][3]] + b_rst[string(obsdata[1][i][5], "_", obsdata[1][i][3])] + b_rb[obsdata[1][i][4]]
        # K is the growth rate for the combination of batch, strain and temperature
        K = K_0 + b_ks[obsdata[1][i][5]] + b_kt[obsdata[1][i][3]] + b_kst[string(obsdata[1][i][5], "_", obsdata[1][i][3])] 

        p = [r, K] # parameters for the logistic growth model

        probh = ODEProblem(logistic_growth, x0i, tspan, p) #remake(prob; x0 = x0i, p = [r, K])
        predicted  = solve(probh, Tsit5(); saveat = obsdata[1][i][2])
    
        # k is time (observation)
        for k in eachindex(predicted)
            obsdata[1][i][1][k] ~  Normal(predicted[k][1], σ^2)
        end 
    end

    return nothing
end


Dat = CSV.read("/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/single_strains/pilot_strain_growth_curves_filtered.tsv", 
    DataFrame, delim='\t')
outdir = "/Users/sur/lab/exp/2025/today3/all_strains_logistic/"

Strains = unique(Dat.strain)
batches = unique(Dat.batch)
temps = unique(Dat.temp)

obsdata = Vector{Any}(undef, 4)
obsdata[1] = Array{Tuple{Vector, Vector, Float64, String7, String7}}(undef, length(unique(string.(Dat.strain, "_", Dat.batch, "_", Dat.temp))))
obsdata[2] = batches
obsdata[3] = float(temps)
obsdata[4] = Strains
 
# Filter data for the specific strain
# Loop through each batch and temperature and make an array of arrays
# with each batch and temperature combination having two arrays: one for
# OD600 and one for total_time_h. Also include in the array the number
# of batches and temperatures.
i = 1
for strain in Strains
    for temp in temps
        for batch in batches
            ii = (Dat.strain .== strain) .& (Dat.batch .== batch) .& (Dat.temp .== temp)
            if(!any(ii))
                continue
            end
            
            obsdata[1][i] = (Dat.OD600[ii], Dat.total_time_h[ii], temp, batch, strain)
            global i += 1
        end
    end

end
obsdata

# Run model parameter infeference
n_samples = 500;
n_warmup = 1000;
model = fit_logistic_all(obsdata);
# map_estimate = maximum_a_posteriori(model)
# map_estimate.values
chain = sample(model, NUTS(),  MCMCThreads(),  n_samples, 4; num_warmup=n_warmup)

# describe(chain)
# hpd(chain; alpha=0.2)
# plot(chain)

# Save posterior
Post = DataFrame(chain)
outfile = joinpath(outdir, "all_logistic_fit.tsv")
CSV.write(outfile, Post; delim='\t')