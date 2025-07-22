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


@model function fit_logistic_multidata_all(obsdata)
    # Prior distributions.
    sigma ~ InverseGamma(2, 2) #Cauchy(0, 2)#
    # sigma_r0 ~ InverseGamma(2, 2) #Cauchy(0, 2)#
    sigma_r0 = 0.5
    sigma_rt ~ InverseGamma(2, 2) #Cauchy(0, 2)#
    sigma_rb ~ InverseGamma(2, 2) #Cauchy(0, 2)#
    # sigma_k0 ~ InverseGamma(2, 2) #Cauchy(0, 2)#
    sigma_k0 = 0.5
    sigma_kt ~ InverseGamma(2, 2) #Cauchy(0, 2)#
    sigma_kb ~ InverseGamma(2, 2) #Cauchy(0, 2)#

    b_rt = Dict()
    b_kt = Dict()
    for temp in obsdata[3]
       b_rt[temp] ~ Normal(0, sigma_rt^2)
       b_kt[temp] ~ Normal(0, sigma_kt^2)
    end

    b_rb = Dict()
    b_kb = Dict()
    for batch in obsdata[2]
       b_rb[batch] ~ Normal(0, sigma_rb^2)
       b_kb[batch] ~ Normal(0, sigma_kb^2)
    end

    K_0 ~ Uniform(0, 2)#LogNormal(log(150), 0.1)
    r_0 ~ Uniform(0, 3)

    # p = [r, K]

    # i is for number of experiment
    for i in eachindex(obsdata[1])
        x0i = [obsdata[1][i][1][1]] # initial condition for the i-th experiment
        tspan = extrema(obsdata[1][i][2])

        r_m = r_0 + b_rt[obsdata[1][i][3]] + b_rb[obsdata[1][i][4]]# r is a function of temperature
        K_m = K_0 + b_kt[obsdata[1][i][3]] + b_kb[obsdata[1][i][4]]# r is a function of temperature

        r ~ LogNormal(r_m, sigma_r0) # r is a function of temperature and batch
        K ~ LogNormal(K_m, sigma_k0) # K is a function of temperature and batch

        p = [r, K] # parameters for the logistic growth model

        probh = ODEProblem(logistic_growth, x0i, tspan, p) #remake(prob; x0 = x0i, p = [r, K])
        predicted  = solve(probh, Tsit5(); saveat = obsdata[1][i][2])
    
        # k is time (observation)
        for k in eachindex(predicted)
            obsdata[1][i][1][k] ~  Normal(predicted[k][1], sigma^2)
        end 
    end

    return nothing
end

Dat = CSV.read("/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/single_strains/pilot_strain_growth_curves_filtered.tsv", 
    DataFrame, delim='\t')
outdir = "/Users/sur/lab/exp/2025/today/single_strain_all_logistic/"

Strains = unique(Dat.strain)
# strain = Strains[2]

for strain in Strains
    # Filter data for the specific strain
    ii = Dat.strain .== strain
    dat = Dat[ii, :]
    println("Fitting strain: $strain with $(sum(ii)) data points")

    # Loop through each batch and temperature and make an array of arrays
    # with each batch and temperature combination having two arrays: one for
    # OD600 and one for total_time_h. Also include in the array the number
    # of batches and temperatures.
    batches = unique(dat.batch)
    temps = unique(dat.temp)
    obsdata = Vector{Any}(undef, 3)
    obsdata[1] = Array{Tuple{Vector, Vector, Float64, String7}}(undef, length(batches) * length(temps))
    obsdata[2] = batches
    obsdata[3] = temps

    i = 1
    for temp in temps
        for batch in batches
            ii = (dat.batch .== batch) .& (dat.temp .== temp)
            obsdata[1][i] = (dat.OD600[ii], dat.total_time_h[ii], temp, batch)
            i += 1
        end
    end

    # Run model parameter infeference
    n_warmup = 1000;
    n_samples = 1500;
    model = fit_logistic_multidata_all(obsdata);
    # map_estimate = maximum_a_posteriori(model)
    # map_estimate.values
    chain = sample(model, NUTS(),  MCMCThreads(),  n_samples, 4; num_warmup=n_warmup)

    # describe(chain)
 
    # hpd(chain; alpha=0.2)

    # plot(chain)

    Post = DataFrame(chain)
    outfile = joinpath(outdir, "$(strain)_all_logistic_fit.tsv")
    CSV.write(outfile, Post; delim='\t')
end
