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
    σ ~ InverseGamma(2, 2) #Cauchy(0, 2)#
    sigma_r0 ~ InverseGamma(2, 2) #Cauchy(0, 2)#
    sigma_rt ~ InverseGamma(2, 2) #Cauchy(0, 2)#

    b_rt = Dict()
    for temp in obsdata[3]
       b_rt[temp] ~ Normal(0, sigma_rt)
    end

    K ~ Uniform(0, 2)#LogNormal(log(150), 0.1)
    r_0 ~ LogNormal(log(1), 1)


    # p = [r, K]

    # i is for number of experiment
    for i in eachindex(obsdata[1])
        x0i = [obsdata[1][i][1][1]] # initial condition for the i-th experiment
        tspan = extrema(obsdata[1][i][2])

        r = r_0 + b_rt[obsdata[1][i][3]] # r is a function of temperature
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
outdir = "/Users/sur/lab/exp/2025/today3/single_strain_all_logisitc/"

Strains = unique(Dat.strain)
strain = Strains[1]

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
    obsdata[1] = Array{Tuple{Vector, Vector, Float64, String7}}(undef, 6)
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
    n_samples = 1000;
    model = fit_logistic_multidata_all(obsdata);
    chain = sample(model, NUTS(),  MCMCThreads(),  n_samples, 4)
    describe(chain)
    # maximum_a_posteriori(model)
    # hpd(chain; alpha=0.2)

    # plot(chain)

    Post = DataFrame(chain)
    # combine(Post, [:r, :K] .=> [mean, median], renamecols = false)
    # combine(Post, [:r] .=> (x -> [quantile(x, (0.1,0.2,0.8,0.9))]) => [:q10, :q20, :q80, :q90], renamecols = false)
    # res = summarize(chain)
    outfile = joinpath(outdir, "$(strain)_logistic_fit.tsv")
    CSV.write(outfile, Post; delim='\t')
end







for i in eachindex(obsdata[1])
    b_
    println(obsdata[1][i][3])
end