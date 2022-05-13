#==============================================================================
file: ps1-codeonly.jl
description: estimating model for Problem Set #1 of Econ 220A (UC Berkeley, 2022),
             a model from Adverse Selection and Inertia in Health Insurance 
             Markets: When Nudging Hurts by Ben Handel (2013). There are two
             companion PDFs describing the problem set and data:
             - ProblemSet_DemandEst_Handel2013-5.pdf
             - Data_Description_Handel_ASIN-1.pdf
author: Aaron C Watt (UCB Grad Student, Ag & Resource Econ)
notes:
    - Some test code and many bulky comments omitted from this version.
==============================================================================#

#==============================================================================
                                PACKAGES
==============================================================================#
# For importing MATLAB file (.mat)
using MAT
# For mean(), std() of matricies
using Statistics
# For dataframes
using DataFrames
# Writing dataframes to CSV
using CSV
# Summary Statistics
using StatsBase
# Random Variable distributions
using Random, Distributions
Random.seed!(0);
# Optimizing functions
using Optim, NLSolversBase
using LinearAlgebra: diag
using JuMP, KNITRO



#==============================================================================
                                 SETUP
==============================================================================#
dir = "/media/a/E/Programming/github/econ220a/problem set 1/"
file = "ASIN-ChoiceModelData-FINAL.mat"
vars = matread(string(dir, file))
# Fix some of the data types
# Integers
for varname in ["Sim", "K", "nIs", "nPlans"]
    vars[varname] = Int(vars[varname])
end


# Read dictionary items into memory
function dict_to_mem(dict)
    println("Reading variables into memory:")
    println("VARNAME \tSIZE")
    for (key, value) in dict
        tabs = 2 - (length(key) ÷ 8)  # ÷ is floor division
        size_ = size(value)
        if size_ == ()
            size_ = (length(value),)
        end
        println(key, "\t"^tabs, size_)
        @eval ($(Symbol(key)) = ($value))
    end
end


# Read all dictionary items into memory (so each variable can be used directly)
println("Adding the following variables to memory from the matlab file:")
dict_to_mem(vars)

macro size(x)  # Easily print size of object with @size on front
    return :( println(size($x)) )
end



#==============================================================================
        QUESTION 1: SUMMARY STATISTICS & TRANSFORM DATA AND VARIABLE CREATION
==============================================================================#
include(string(dir,"ps1_summarystats.jl"))
include(string(dir,"ps1_transformdata.jl"))



#==============================================================================
        QUESTION 2 & 3: LOG LIKELIHOOD FUNCTION
==============================================================================#
# Simulation functions

function simulate_risk_dist(; α_mean, α_sd, α_income, α_age, rand_coef)
    sd = α_sd/100 * rand_coef  # dim = 1 x Sim
    sd = repeat(sd, nIs)       # dim = nIs x Sim
    means = α_mean/100 .+ α_income/100*Inc_perm .+ α_age/100*Ages_max  # dim = nIs x 1
    means = repeat(means, outer=(1, Sim))                              # dim = nIs x Sim
    dist = means .+ sd
    dist = repeat(dist, outer=(1, 1, K))  # dim = nIs x Sim x K
    return max.(dist, 0.000001)
end


function simulate_ppo1200_dist(; αm_sgl, αsd_sgl, αm_fam, αsd_fam, rand_coef)
    dist = IND' .* (αm_sgl .+ αsd_sgl*rand_coef) .+ (1 .- IND') .* (αm_fam .+ αsd_fam*rand_coef)  # dim = nIs x Sim
    dist = repeat(dist, outer=(1, 1, K))  # dim = nIs x Sim x K
    return max.(dist, 0.000001)
end


"""Create ​ ϵₖⱼₜ with sd = σ​²_ϵⱼ(​Yₖ) (vector of family-plan-time specific shocks)
"""
function simulate_ϵ(α3)
    α500 = 2(1 .- IND')  .+ 1  # vector of indices to get from α3
    α500 = getindex(α3, Int.(α500))  # α(3,1) or α(3,3)
    α1200 = 2(1 .- IND')  .+ 2  # vector of indices to get from α3
    α1200 = getindex(α3, Int.(α1200))  # α(3,2) or α(3,4)
    ϵₖⱼₜ = Dict("500" => Dict(  # plan = 500
            "1" => repeat(α500 * ϵ[1, :]', outer=(1, 1, K)),  # year 1
            "2" => repeat(α500 * ϵ[2, :]', outer=(1, 1, K)),
            "3" => repeat(α500 * ϵ[3, :]', outer=(1, 1, K))),
        "1200" => Dict(  # plan = 1200
            "1" => repeat(α1200 * ϵ[4, :]', outer=(1, 1, K)),  # year 1
            "2" => repeat(α1200 * ϵ[5, :]', outer=(1, 1, K)),
            "3" => repeat(α1200 * ϵ[6, :]', outer=(1, 1, K)))
    )
    return ϵₖⱼₜ
end


function simulate_inertia(α)
    # time-constant (for years 2,3) family vectors (nIs x 1)
    η₀ = α[4,2]*IND' .+ α[4,3]*(1 .- IND')
    INCOME = α[5,1]*Inc2  # No Inc for year 3, assume Inc3 = Inc2
    SOPHIST = α[5,2]*QS
    MANAGER = α[5,3]*managerX
    CHRONIC = α[5,4]*CC2  # No CC for year 3, assume CC3 = CC2

    η = Dict()
    for year ∈ 2:3, plan ∈ plans
        FSA = α[4,4] * vars["FSAY$year"]
        ΔEXPEND = α[6,1] * vars["CSAL$year"]
        CHOICE = mat_dict["choice$plan$year"]
        η_ = η₀ + INCOME + SOPHIST + MANAGER + CHRONIC + FSA + ΔEXPEND
        η_ = repeat(η_, outer=(1,Sim)) .* CHOICE
        η["$plan$year"] = repeat(η_, outer=(1,1,K))
    end
    return η
end


# PART A: 
"""Return the negative log likelihood from the input parameter matrix α.
"""
function nLL(α)
    # PART B: 
    γₖ = simulate_risk_dist(α_mean = α[1,1],
                            α_income = α[1,2], 
                            α_age = α[1,3], 
                            α_sd = α[1,4], 
                            rand_coef = randcoef_risk) # this is defined in ps1_transformdata.jl
    
    # PART C:
    δₖ₁₂₀₀ = simulate_ppo1200_dist(αm_sgl = α[2,1],
                                      αsd_sgl = α[2,2], 
                                      αm_fam = α[2,3], 
                                      αsd_fam = α[2,4], 
                                      rand_coef = randcoef_1200)

    # PART D: 
    ϵₖⱼₜ = simulate_ϵ(α[3, :])

    # PART E:
    H = α[4,1] .* HTCi

    # PART F:
    η = simulate_inertia(α)


    # QUESTION 3: Utility-choice estimation ======================================
    # PART A
    Usim = Dict()
    for plan ∈ plans, year ∈ 1:3
        if year > 1 && plan != "250"
            x = mat_dict["EUvector$plan$year"] + η["$plan$year"] + δₖ₁₂₀₀ + H + ϵₖⱼₜ[plan]["$year"]
        elseif year > 1 && plan == "250"  # omit ϵ
            x = mat_dict["EUvector$plan$year"] + η["$plan$year"] + δₖ₁₂₀₀ + H
        elseif year == 1 && plan != "250"  # omit η
            x = mat_dict["EUvector$plan$year"] + δₖ₁₂₀₀ + H + ϵₖⱼₜ[plan]["$year"]
        else  # 250-1: omit ϵ & η
            x = mat_dict["EUvector$plan$year"] + δₖ₁₂₀₀ + H
        end
        Usim["$plan$year"] = -exp.(-γₖ .* x) ./ γₖ
    end

    # PART B
    Umean = Dict()
    for year ∈ 1:3
        Umean["$year"] = zeros(nIs, Sim, nPlans)
        for (i, plan) ∈ enumerate(plans)
            Umean["$year"][:, :, i] = mean(Usim["$plan$year"], dims=3)
        end
    end

    # Normalize by EU of PPO1200 for each year
    for year ∈ 1:3, i ∈ 1:3
        Umean["$year"][:, :, i] = Umean["$year"][:, :, 3] ./ Umean["$year"][:, :, i]
    end

    # PART C
    Usum = Dict()
    for year ∈ 1:3
        s = sum(Umean["$year"].^6, dims=3)
        Usum["$year"] = Umean["$year"].^6 ./ s
    end

    # PART D
    probs = ones(nIs, Sim)
    for year ∈ 1:3
        yearprobs = zeros(nIs, Sim)
        for plan ∈ 1:3
            yearplanindicator = repeat(Int.(choice[:, year] .== plan), outer=(1, Sim))
            yearprobs .+= Usum["$year"][:,:,plan] .* yearplanindicator
        end
        probs .*= yearprobs
    end

    # Take the average over all simulations (along the Sim dimension)
    meanprobs = mean(probs, dims=2)
    
    # return the log-likelihood
    value = -sum(log.(meanprobs))
    # println("Function Value: $(round(value))")
    return value
end



#==============================================================================
        QUESTION 4: OPTIMIZE OVER THE LIKELIHOOD
==============================================================================#
"""Wrapper to pass individual parameter elements to likelihood function"""
function nLL_wrapper(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21)
    α = [
        a1   a7  a12  a17
        a2   a8  a13  a18
        a3   a9  a14  a19
        a4  a10  a15  a20
        a5  a11  a16  a21
        a6    0.    0.    0.
    ]
    return nLL(copy(α))
end

"""Wrapper to pass 21x1 parameter vector to likelihood function"""
function nLL_vec(a::Vector)
    return nLL_wrapper(a...)
end

"""Transform 21 input arguments to a 6x4 matrix"""
function transform_to_mat(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21)
    α = [
        a1   a7  a12  a17
        a2   a8  a13  a18
        a3   a9  a14  a19
        a4  a10  a15  a20
        a5  a11  a16  a21
        a6    0.    0.    0.
    ]
    return copy(α)
end

"""Transform 6x4 matrix to a 21x1 vector"""
function transform_to_21vec(x::AbstractMatrix)
    y = reshape(copy(x), 24, 1)
    # remove the (6,2), (6,3), and (6,4) elements
    y = deleteat!(vec(y), [12, 18, 24])
    return y
end



# PART A
lb = [
          0.0     -0.5     -0.5      0.0
        -Inf       0.1    -Inf       0.1
          1.0      1.0      1.0      1.0
        -Inf       0.0      0.0  -5000.0
      -5000.0  -5000.0  -5000.0  -5000.0
      -5000.0      0.0      0.0      0.0
]

α₀ = [
        0.06     0.003      0.0     0.04
    -2500.0    700.0    -2200.0   800.0
      300.0    800.0      300.0   800.0
     -500.0   1250.0     1750.0  -500.0
        0.0      0.0        0.0     0.0
        0.0      0.0        0.0     0.0
]


ub = [
        2.0      1.0      1.0      1.0
        Inf      Inf      Inf      Inf
    20000.0  20000.0  20000.0  20000.0
        Inf   20000.0  20000.0   5000.0
    5000.0   5000.0   5000.0   5000.0
    5000.0      0.0      0.0      0.0
]



# PART B 
nLL(α₀)  # 5525.114607540947


alpha0 = α₀

# Other starting values (modified from replication code to fit in bounds) 

alpha1 = [
    0.1    0.04     0.0     0.03
    900.0    0.1   1500.0   500.0
    800.0  300.0    800.0   300.0
  -2500.0  700.0   1500.0  -200.0
     20.0  100.0    100.0  -100.0
   -100.0    0.0      0.0     0.0
]
alpha2 = [
    0.3    0.1    -0.0003     0.08
    600.0    0.1  2500.0      300.0
    500.0  200.0   500.0      200.0
  -1500.0  700.0  1200.0     -800.0
    200.0    0.0     0.0     -200.0
   -100.0    0.0     0.0        0.0
]
alpha3 = [
    0.5    0.02    -0.0003     0.03
    1500.0    0.1   1700.0        1.0
     400.0  150.0    400.0      150.0
   -1000.0  500.0   1500.0     -500.0
       0.0   50.0    200.0     -100.0
     100.0    0.0      0.0        0.0
]
alpha4 = [
    0.9    0.05     0.003     0.06
    1000.0    0.1   1500.0     100.0
     150.0  150.0    150.0     150.0
   -2000.0  800.0   1050.0    -700.0
     -50.0  -50.0    100.0     100.0
     -50.0    0.0      0.0       0.0
]


# PART C
# Time limit for each of 5 iterations to finish in 8 hours (in seconds)
time_limit5 = round(8*3600/5)
max_iter = 40_000

#! Save results from 5 different starting parameter vectors using SAMIN algorithm
results = Dict(); i=0;
for a ∈ [alpha0, alpha1, alpha2, alpha3, alpha4]
    resSM = Optim.optimize(nLL_vec, transform_to_21vec(lb), transform_to_21vec(ub), 
        transform_to_21vec(a), SAMIN(), 
        Optim.Options(time_limit=time_limit5, iterations=max_iter))
    results[i] = resSM
    open("SAMIN-log.txt", "a") do io
        write(io, "\n"^6 * "="^30 * " alpha$i" * "="^30)
        write(io, "\nStarting Value:\n $(a)")
        write(io, "\nFinal Value:\n $(transform_to_mat(Optim.minimizer(resSM)...))")
        write(io, "\n\nresult:\n $(resSM)")
    end
    i += 1
end



#=========
KNITRO
=========#
using KNITRO
function callbackEvalF(kc, cb, evalRequest, evalResult, userParams)
    α = evalRequest.x
    # Evaluate nonlinear objective
    evalResult.obj[1] = nLL_vec(α)

    return 0
end

# function callbackEvalG!(kc, cb, evalRequest, evalResult, userParams)
#     x = evalRequest.x

#     # Evaluate gradient of nonlinear objective
#     dTmp = x[2] - x[1] * x[1]
#     evalResult.objGrad[1] = (-400.0 * dTmp * x[1]) - (2.0 * (1.0 - x[1]))
#     evalResult.objGrad[2] = 200.0 * dTmp

#     return 0
# end

# Create a new Knitro solver instance.
kc = KNITRO.KN_new()

n = 21
KNITRO.KN_add_vars(kc, n)


_inf = KNITRO.KN_INFINITY
lb2 = replace(lb, -Inf => -_inf, Inf => _inf)
α₀2 = replace(α₀, -Inf => -_inf, Inf => _inf)
ub2 = replace(ub, -Inf => -_inf, Inf => _inf)

KNITRO.KN_set_var_lobnds(kc,  transform_to_21vec(lb2)) # not necessary since infinite
KNITRO.KN_set_var_upbnds(kc,  transform_to_21vec(ub2))
# Define an initial point.  If not set, Knitro will generate one.
KNITRO.KN_set_var_primal_init_values(kc, transform_to_21vec(α₀2))

# Set verbose printing level
KNITRO.KN_set_param(kc, "outlev", 6)

# Enable multi-start (See notes below)
# KNITRO.KN_set_param(kc, "ms_enable", 1)

# Set multi-threading parallelism (<0 will automagically find a good number of threads to use)
KNITRO.KN_set_param(kc, "numthreads", 16)
# KNITRO.KN_set_param(kc, "ms_numthreads", 8)

cb = KNITRO.KN_add_eval_callback_all(kc, callbackEvalF)
# KNITRO.KN_set_cb_grad(kc, cb, callbackEvalG!)
nStatus = KNITRO.KN_solve(kc)

# Delete the Knitro solver instance.
KNITRO.KN_free(kc)









#==========================================================
PART D
Incorporating demographic heterogeneity and the parameter es-
timates, what is the population distribution of inertial costs?
==========================================================#
#Evaluate η(α) at the Knitro solution,
# and estimate the distribution... for plot a kernal density?
A5 = transform_to_mat(a5...)
η5 = simulate_inertia(A5)
ηcomb = Dict(2 => zeros(nIs), 3 => zeros(nIs))
for plan ∈ plans, year ∈ [2,3]
    ηcomb[year] .+= η5["$plan$year"][:,1,1] ./ 1000
end
df2 = DataFrame(year = "year 2", η = ηcomb[2])
df3 = DataFrame(year = "year 3", η = ηcomb[3])
df = vcat(df2, df3)
p4d = @df df density(:η, group = (:year), legend = :topleft)
xlabel!("Interia Costs (\$\$1000)")
ylabel!("kernel density")
title!("Density Plots of Inertial Costs, by year")
savefig(p4d, "4-d.svg")




#==========================================================
PART E
Consider the mean of the CARA risk preference distribution of
random coefficients. Translate this coefficient into the value of X
that makes a family indifferent between no gamble and a gamble
where they win $100 with 50% probability and lose $X with 50%
probability.
==========================================================#
# mean of the calculated coefficients for each family
risks5 = simulate_risk_dist(α_mean = A5[1,1],
                            α_income = A5[1,2], 
                            α_age = A5[1,3], 
                            α_sd = A5[1,4], 
                            rand_coef = randcoef_risk
                            )[:,:,1] # last index is just a copy

