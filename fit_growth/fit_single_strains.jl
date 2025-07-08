using CSV
using DataFrames
using Turing
using DifferentialEquations
using StatsPlots
using Random
using MCM


Random.seed!(544);


# Define logistic model.
function logistic_growth(dx, x, p, t)
    # Model parameters.
    r, K = p

    # Evaluate differential equations.
    dx[1] = r*x[1]*( 1 - x[1]/K)

    return nothing
end


@model function fit_logistic_multidata(obsdata)
    # Prior distributions.
    σ ~ InverseGamma(2, 2) #Cauchy(0, 2)# 

    K ~ Uniform(0, 2)#LogNormal(log(150), 0.1)
    r ~ Uniform(0, 3) #LogNormal(log(1), 0.1)


    p = [r, K]

    # i is for number of experiment
    for i in eachindex(obsdata)
        x0i = [obsdata[i][1][1]] # initial condition for the i-th experiment
        tspan = extrema(obsdata[i][2])
        probh = ODEProblem(logistic_growth, x0i, tspan, p) #remake(prob; x0 = x0i, p = [r, K])
        predicted  = solve(probh, Tsit5(); saveat = obsdata[i][2])
    
        # k is time (observation)
        for k in eachindex(predicted)
            obsdata[i][1][k] ~  Normal(predicted[k][1], σ^2)
        end 

    end

    return nothing
end


Dat = CSV.read("/Users/sur/lab/exp/2025/today3/pilot_strain_growth_curves_filtered.tsv", 
    DataFrame, delim='\t')
outdir = "/Users/sur/lab/exp/2025/today3/single_strain_logistic/"

Strains = unique(Dat.strain)

strain = Strains[1]


# Filter data for the specific strain
ii = Dat.strain .== strain
dat = Dat[ii, :]

# Loopm through each batch and temperature and make an array of arrays
# with each batch and temperature combination having two arrays: one for
# OD600 and one for total_time_h
batches = unique(dat.batch)
temps = unique(dat.temp)
obsdata = Array{Tuple{Vector, Vector}}(undef, length(batches) * length(temps))
i = 1
for temp in temps
    for batch in batches
        ii = (dat.batch .== batch) .& (dat.temp .== temp)
        # println(sum(ii))
        obsdata[i] = (dat.OD600[ii], dat.total_time_h[ii]) 
        i += 1
    end
end

# Run model parameter infeference

n_samples = 1000;
model = fit_logistic_multidata(obsdata);
chain = sample(model, NUTS(),  MCMCThreads(),  n_samples, 4)
# describe(chain)
# maximum_a_posteriori(model)
# hpd(chain; alpha=0.2)

# plot(chain)

Post  = DataFrame(chain)
# combine(Post, [:r, :K] .=> [mean, median], renamecols = false)
# combine(Post, [:r] .=> (x -> [quantile(x, (0.1,0.2,0.8,0.9))]) => [:q10, :q20, :q80, :q90], renamecols = false)
# res = summarize(chain)
outfile = joinpath(outdir, "$(strain)_logistic_fit.tsv")
CSV.write(outfile, Post; delim='\t')
