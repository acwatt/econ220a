#==============================================================================
file: ps1.jl
description: estimating model for Problem Set #1 of Econ 220A (UC Berkeley, 2022),
             a model from Adverse Selection and Inertia in Health Insurance 
             Markets: When Nudging Hurts by Ben Handel (2013). There are two
             companion PDFs describing the problem set and data:
             - ProblemSet_DemandEst_Handel2013-5.pdf
             - Data_Description_Handel_ASIN-1.pdf
author: Aaron C Watt (UCB Grad Student, Ag & Resource Econ)
notes:
    - Finished summary stats Q1 a & b
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
# Plotting
using StatsPlots
# Optimizing functions
using Optim, NLSolversBase
using LinearAlgebra: diag
using JuMP, KNITRO
# Parallelize!
# using Distributed
# addprocs(2)

# # Ensure BlackBoxOptim loaded on all workers
# @everywhere using BlackBoxOptim



#==============================================================================
                                 SETUP
==============================================================================#
dir = "/media/a/E/Programming/github/econ220a/problem set 1/"
file = "ASIN-ChoiceModelData-FINAL.mat"
# Open MATLAB datafile as a dictionary
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
#                        QUESTION 1: SUMMARY STATISTICS
==============================================================================#

#==========================================================
PART A & B: EXPENDITURES & TABULATIONS
==========================================================#
# Run the code in ps1_summarystats.jl to save summary stats to CSV
include(string(dir,"ps1_summarystats.jl"))
#! CANNOT RUN THIS AGAIN AFTER RUNNING ps1_transformdata.jl
# ps1_transformdata.jl transforms the variables
# Must reload using lines above if summary stats need to be rerun



#==========================================================
PART C & D: TRANSFORM DATA AND VARIABLE CREATION
==========================================================#
# Run the code in ps1_transformdata.jl to create objects in memory
# to use during the simulations / estimation
include(string(dir,"ps1_transformdata.jl"))



























#==========================================================
                        FUNCTIONS
==========================================================#

elmax(x, a) = [maximum([i, a]) for i in x]

function simulate_risk_dist(; α_mean, α_sd, α_income, α_age, rand_coef)
    #= Use one of the 1xSim matrices generated in 1c and multiply
    by alpha(1,4)/100. Fill in an (nIs,Sim) matrix with this same
    row for all families. This is the standard deviation of risk prefer-
    ence random coefficients for candidate parameter alpha(1,4)/100.=#
    sd = α_sd/100 * rand_coef  # dim = 1 x Sim
    sd = repeat(sd, nIs)       # dim = nIs x Sim
    #= generate an (nIs,1) matrix of family specific means where
    each entry is the mean of the risk preference random coefficient
    distribution for each family. This should equal alpha(1,1)/100 +
    alpha(1,2)/100*Income + alpha(1,3)/100*MaxAge
    Transform into an (nIs,Sim) matrix which repeats this columns 
    vector Sim times.=#
    means = α_mean/100 .+ α_income/100*Inc_perm .+ α_age/100*Ages_max  # dim = nIs x 1
    means = repeat(means, outer=(1, Sim))                              # dim = nIs x Sim
    #=Now, add the ‘means’ matrix to the ’standard deviation’
    matrix: you’ve just simulated a distribution of risk preference for
    each family.
    Take this (nIs,Sim) matrix and transform it into an
    (nIs,K,Sim) matrix by repeating the (nIs,Sim) matrix K times=#
    dist = means .+ sd
    dist = repeat(dist, outer=(1, 1, K))  # dim = nIs x Sim x K
    # maximum of matrix and 0.000001, truncating coefficients above 0
    return max.(dist, 0.000001)
end


function simulate_ppo1200_dist(; αm_sgl, αsd_sgl, αm_fam, αsd_fam, rand_coef)
    # IND' = vector of 1=single, 0=family. Family Status
    # For each row (family), need to create mean + sd*rand_coef, depending on Family Status
    dist = IND' .* (αm_sgl .+ αsd_sgl*rand_coef) .+ (1 .- IND') .* (αm_fam .+ αsd_fam*rand_coef)  # dim = nIs x Sim
    dist = repeat(dist, outer=(1, 1, K))  # dim = nIs x Sim x K
    # maximum of matrix and 0.000001, truncating coefficients above 0
    return max.(dist, 0.000001)
end


"""Create ϵₖⱼₜ with sd = σ²_ϵⱼ(Yₖ) (vector of family-plan-time specific shocks)
    k = family, j = plan, t = year. 
    For each plan-year, create an nIs x Sim matrix of family-specific shocks,
    based on family status IND (nIs x 1) and ϵ matricies (1 x Sim)
    Then repeat this matrix K times
"""
function simulate_ϵ(α3)
    # m = 1  # indexing the ϵ matricies
    # for year ∈ 1:3
    #     for plan ∈ ["500", "1200"]
    #         ϵⱼ = ϵ[m, :]'  # dim = 1 x Sim
    #         # Get vector of α standard deviations for this plan / family status
    #         αidx = (plan == "1200") .+ 2(1 .- IND')  .+ 1
    #         αⱼ = getindex(α3, Int.(αidx))  # dim = nIs x 1
    #         # Create (nIs,K,Sim) matrix for this plan / family status
    #         # IND' = vector of 1=single, 0=family. Family Status
    #         # For each row (family), need to create mean + sd*rand_coef, depending on Family Status
    #         dist = αⱼ * ϵⱼ  # dim = nIs x Sim
    #         dist = repeat(dist, outer=(1, 1, K))  # dim = nIs x Sim x K
    #         # Assign the new matrix to the varname in memory
    #         varname = string("ϵ$(plan)_$(year)")
    #         @eval (($(Symbol(varname))) = ($dist))
    #         m += 1
    #     end
    # end
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
    # η₀, η₂: time-constant (for years 2,3) family vectors (nIs x 1)
    η₀ = α[4,2]*IND' .+ α[4,3]*(1 .- IND')
    # η₁ * X
    INCOME = α[5,1]*Inc2  # No Inc for year 3, assume Inc3 = Inc2
    SOPHIST = α[5,2]*QS
    MANAGER = α[5,3]*managerX
    CHRONIC = α[5,4]*CC2  # No CC for year 3, assume CC3 = CC2

    η = Dict()
    for year ∈ 2:3, plan ∈ plans
        FSA = α[4,4] * vars["FSAY$year"]
        ΔEXPEND = α[6,1] * vars["CSAL$year"]
        CHOICE = mat_dict["choice$plan$year"]

        # Linear combination (Eq on pg 2665) ()
        η_ = η₀ + INCOME + SOPHIST + MANAGER + CHRONIC + FSA + ΔEXPEND

        # Multiply each of the 6 matrices element-wise by 
        # their respective choice‘p’‘t’ matrix so the only non-zero
        # elements in the matrix for plan p in year t are the 
        # rows where the family chose plan p in year t-1
        η_ = repeat(η_, outer=(1,Sim)) .* CHOICE
        η["$plan$year"] = repeat(η_, outer=(1,1,K))
    end
    return η
end





#==============================================================================
QUESTION 2: THE LIKELIHOOD FUNCTION

Now we will move to the actual estimation. You will be estimating
the primary choice model in Handel (2013) referenced above. You
will need to reference section three of the paper in order to under-
stand the exercise. The first part of that section describes the choice
model while the Estimation subsection describes specific assumptions
made in estimation (specifically, distributional assumptions for unob-
served heterogeneity and functions linking observed demographics to
risk preferences, inertia, PPO1200 preferences, and epsilons). The
discussion of the model and estimation in the paper helps to provide
context for the detailed implementation that this problem set walks
you through. Note: the cost model described in that section is al-
ready estimated in your data and will be used as an input into choice
model estimation. You will also want to reference Online Appendix
B on the author’s website with details on the simulated maximum
likelihood estimation algorithm.

The remainder of the problem set has the following major components:
    - Set up the likelihood function that describes the likelihood that
      families make certain sequence of choices, subject to candidate
      model parameters
    - Set up and run a non-linear optimization routine that searches
      over the parameter space to find the model parameters that max-
      imize the value of this likelihood function. These are the final
      estimates.

First we’ll set up the likelihood function then the non-linear optimiza-
tion that takes that function as an input. For candidate parameters
considered in non-linear optimization, this likelihood function simu-
lates choices as if those candidate parameters are the true parameters,
then matches the predicted choices to actual choices made. This func-
tion is called in the optimizer you will build in the next question and
when it returns the best value, those are the parameters that simulate
choices that best match the observed choices. Before you complete
questions 2 and 3 setting up the likelihood function, you may want to
take a quick read through question 4 to see how the likelihood function
will be used. For example, section 4 defines an initial parameter value
‘guess’ alpha0 that will be used to start the non-linear optimization
routine. The likelihood function instructions below describe what each
parameter to be estimated is, and which entry in alpha that parameter
corresponds to.
==============================================================================#


#==========================================================
PART A: Define the function

Here, alpha will be your matrix of parameters that you will es-
timate, and after alpha you should include all variables you will
need to bring as fixed factors into the estimation, separated by
commas. This includes the matrices generated above in ’Esima-
tionCode’ and most of the demographic variables in the main
data.
==========================================================#
"""nLL(α)

    Return the negative log likelihood from the input parameter matrix α.

    α is a 6x4 matrix defined as:

    Row 1:
        Col 1: intercept of the mean of the normal distribution of risk preferences  α[1,1]=2.32e-4
        Col 2: amount this mean shifts by if income tier increases by 1              α[1,2]=2.9e-5
        Col 3: amount the mean shifts by if family max. age increases by 1                  2.27e-6
        Col 4: standard deviation of the risk preference random coefficient                 1.88e-4
            distribution around the mean formed with these prior three items

    Row 2: Used to generate random coefficient distributions for PPO1200 using the other vector of 1 X Sim draws in 1c
        Col 1: mean for singles                  -2912
        Col 2: standard deviation for singles    843
        Col 3: mean for families                 -2871
        Col 4: standard deviation for families   897

    Row 3: Used to create six (nIs,K,Sim) matrices for the mean 0 normal epsilon draws for PPO500 and PPO1200 over the three years.
        Col 1: standard deviations for individuals for PPO500   204
        Col 2: standard deviations for individuals for PPO1200  502
        Col 3: standard deviations for families for PPO500      329
        Col 4: standard deviations for families for PPO1200     811

    Row 4: Used to references for individuals with high total costs, and inertia for each family.
        Col 1: constant preference for all people with HTC=1                   856
        Col 2: inertia intercept for individuals (η₀+η₂)                α[4,2]=2480
        Col 3: inertia intercept for families (η₀)                      α[4,3]=1729
        Col 4: inertia coeficient for FSA Enrollment (η₁)               α[4,4]=-551

    Row 5: 
        Col 1: inertia coeficient for Income (η₁)                       α[5,1]=-32
        Col 2: inertia coeficient for Quantitative Sophistication Indicator (η₁)  α[5,2]=5
        Col 3: inertia coeficient for Manager Indicator (η₁)            α[5,3]=198
        Col 4: inertia coeficient for Chronic Conditions Indicator (η₁) α[5,4]=80

    Row 6: 
        Col 1: inertia coeficient for Large Expenditure Change in Prior Year (η₁)  α[6,1]=156
        Col 2: 
        Col 3: 
        Col 4: 

    The estimated α from the paper is
    α = [2.32e-4 2.90e-5 2.27e-6 1.88e-4
          -2912     843   -2871     897
            204     502     329     811
            856    2480    1729    -551
            -32       5     198      80
            156       0       0       0 ]

"""
function nLL(α)

    #==========================================================
    PART B: create matrices describing risk preferences, 
            PPO1200 coefficient, and epsilons

    The first step is to take the candidate parameters in alpha (which
    will be optimized over) and create matrices describing risk pref-
    erences, PPO1200 coefficient, and epsilons as functions of those
    parameters. These will be used later in the file to simulate choices
    conditional on those parameters, which will then be matched to
    the data. I’ll help you get started: define alpha(1,1)/100 as the
    intercept of the mean of the normal distribution of risk prefer-
    ences, alpha(1,2)/100 as the amount this mean shifts by if income
    tier increases by 1, alpha(1,3)/100 is the amount the mean shifts
    by if family max age increases by 1, while alpha(1,4)/100 is the
    standard deviation of the risk preference random coefficient dis-
    tribution around the mean formed with these prior three items.

    Use one of the 1xSim matrices generated in 1c and multiply
    by alpha(1,4)/100. Fill in an (nIs,Sim) matrix with this same
    row for each family. This is the standard deviation of risk prefer-
    ence random coefficients for candidate parameter alpha(1,4)/100.
    Next, generate an (nIs,1) matrix of family specific means where
    each entry is the mean of the risk preference random coefficient
    distribution for each family. This should equal alpha(1,1)/100 +
    alpha(1,2)/100*Income + alpha(1,3)/100*MaxAge. Transform
    into an (nIs,Sim) matrix which repeats this columns vector Sim
    times. Now, add the ‘means’ matrix to the ’standard deviation’
    matrix: you’ve just simulated a distribution of risk preference for
    each family. Take this (nIs,Sim) matrix and transform it into an
    (nIs,K,Sim) matrix by repeating the (nIs,Sim) matrix K times.
    The reason for this will make sense later when you simulate ex-
    pected utilities: it helps speed up those calculations. Finally,
    make sure to set this matrix equal to the maximum of itself and
    0.000001, truncating risk preference coefficients above 0 as the
    theory specifics.
    ==========================================================#
    # α = [2.23e-4 2.90e-5 2.27e-6 1.88e-4
    #        -2912     843   -2871     897
    #          204     502     329     811
    #          856    2480    1729    -551
    #          -32       5     198      80
    #          156       0       0       0]
    # Create (nIs,K,Sim) matrix of random coefficient distributions for risk preference
    γₖ = simulate_risk_dist(α_mean = α[1,1],
                            α_income = α[1,2], 
                            α_age = α[1,3], 
                            α_sd = α[1,4], 
                            rand_coef = randcoef_risk) # this is defined in ps1_transformdata.jl
    
    
    #==========================================================
    PART C: generate the random coefficient distributions for 
            PPO1200 using the other vector of 1 X Sim draws in 1c
    
    I’m going to be less helpful for the next few steps. First, gener-
    ate the random coefficient distributions for PPO1200 using the
    other vector of 1 X Sim draws in 1c and parameters alpha(2,1), al-
    pha(2,2) (mean and standard deviation for singles) and alpha(2,3),
    alpha(2,4) (mean and standard deviation for families). Create a
    (nIs,K,Sim) matrix similar to what you did for risk preferences for
    these draws. Note: YOU have to make sure to fill in this matrix
    appropriately as a function of family status, use the IND indica-
    tor variable to determine ‘single’ or ‘family’ (or Tier2 variable)
    NOT the size of Ages and Genders, which doesn’t map correctly.
    ==========================================================#
    # Create (nIs,K,Sim) matrix of random coefficient distributions for PPO1200
    δₖ₁₂₀₀ = simulate_ppo1200_dist(αm_sgl = α[2,1],
                                      αsd_sgl = α[2,2], 
                                      αm_fam = α[2,3], 
                                      αsd_fam = α[2,4], 
                                      rand_coef = randcoef_1200)

    
    #==========================================================
    PART D: Create six (nIs,K,Sim) matrices for the mean 0 
            normal epsilon draws for PPO500 and PPO1200 over 
            the three years.
    
    Create six (nIs,K,Sim) matrices for the mean 0 normal epsilon
    draws for PPO500 and PPO1200 over the three years. Let al-
    pha(3,1) and(3,2) be the standard deviations for individuals and
    alpha(3,3) and alpha(3,4) for families. Use the 6 X Sim matrix
    generated in 1c and procedure similar to that for the PPO1200
    random coefficients to generate the final six matrices. Each ma-
    trix should correspond to either PPO1200 or PPO500, and to
    either year 1,2,or 3 (since these draws vary for an individual over
    time).
    ==========================================================#
    # Create dictionary of 6 normal ϵ matrices with keys = plan and year 
    ϵₖⱼₜ = simulate_ϵ(α[3, :])

    
    #==========================================================
    PART E: Create a matrix that describes preferences for 
            individuals with high total costs

    Create a matrix that describes preferences for individuals with
    high total costs, as described by the matrix you transformed in
    1c. This preference is constant for all people with HTCi = 1:
    use parameter (4,1) and create an (nIs,K,Sim) matrix = 0 from
    people with HTCi=0 and = alpha(4,1) with HTCi=1. This will
    be added to PPO500 and PPO1250 money at stake up front later
    (so reflects the disutility of high-cost people for those plans).
    ==========================================================#
    H = α[4,1] .* HTCi


    #==========================================================
    PART F: fill in inertia for each family.

    Finally, before simulating choices, you need to fill in inertia for
    each family. Use alpha(4,2) and alpha(4,3) to be individual and
    family inertia intercepts. Then, using observables demographics,
    use alpha(4,4), alpha(5,1), alpha(5,2), alpha(5,3), alpha(5,4), and
    alpha(6,1) to index inertia as a function of observed demographics:
        FSA enrollment (1 or 0)                            α[4,4]
        Income (1-5)                                       α[5,1]
        Quantitative Sophistication Indicator (1 or 0)     α[5,2]
        Manager Indicator (1 or 0)                         α[5,3]
        Chronic Conditions Indicator (1 or 0)              α[5,4]
        Large Expenditure Change in Prior Year (1 or 0)    α[6,1]

    Generate inertia preference matrices for both year two and year
    three (some of these indicators changes over time!!) that are
    (nIs,Sim) using the above variables and parameters (use the lin-
    ear format described in the Estimation section in the paper).
    Then, you will create 6 (nIs,K,Sim) matrices describing inertia
    preferences for each plan in years 2 and 3 (the non-active choice
    years). Use the choice‘p’‘t’ matrices created in 1d (inputs into the
                    ^ nIs x Sim
    likelihood function!) to multiply these matrices such that inertia
    benefits the incumbent plan but has no benefit for plans which
    are alternative options in the observed data. Thus, only one plan
    for each family can get a benefit from incumbency in each year,
    and you should take that into account in these matrices, which
    will be added to utility.
    - repeat nIs x Sim matrices K times to get nIs x K x Sim
    ==========================================================#
    # Generate dictionary of inertia preference matrices
    η = simulate_inertia(α)

   





    




    #==============================================================================
    QUESTION 3: generate choice predictions conditional on candidate parameters
                within the likelihood function, and match them to observed choices.
    ==============================================================================#

    #==========================================================
    PART A
    Create state by state (of health draws) v-NM utility values. To
    do this compute u_k(x) as specified in the paper:
            uₖ(x) = -γₖ(Xᴬₖ)⁻¹ e^[-γₖ(Xᴬₖ) x]
    where
            x = Wₖ - Pₖⱼₜ - OOP + η(Xᴮₖₜ, Yₖ) 1ₖⱼ,ₜ₋₁ + δₖ(Yₖ)1₁₂₀₀ + αHₖ,ₜ₋₁ 1₂₅₀ + ϵₖⱼₜ(Yₖ)
    See page 15 of the paper for an explanation of notation and the
    estimation section for further implementation details.
    Note: In this problem set, α should enter as a coefficient on PPO500 and 1250 (same
    for both plans) rather than on PPO250, since this it how it enters in the sample code.

    You should use your EUvector‘p’‘t’ matrices from 1d, which are
    inputs into the likelihood function, here. Those (nIs,K,Sim) ma-
    trices describe W-OOP-P for each state / health realization and
    health plan, which is part of x above that is fixed regardless of
    choice model parameters. You should use these as inputs to com-
    puting the (nIs,K,Sim) matrices of v-NM utility values above.
    Try to do this without loops, which slow down programs sub-
    stantially!
    --> loops in Julia are fast :)
    
    x = Wₖ - Pₖⱼₜ - OOP + η(Xᴮₖₜ, Yₖ) 1ₖⱼ,ₜ₋₁ + δₖ(Yₖ)1₁₂₀₀ + αHₖ,ₜ₋₁ 1₂₅₀ + ϵₖⱼₜ(Yₖ)
    x = EUvector‘p’‘t’  + η["$plan$year"]     + δₖ₁₂₀₀      + H            + ϵₖⱼₜ[plan][year]

    uₖ(x) = -ℯ^[-γₖ(Xᴬₖ) x] / γₖ(Xᴬₖ)
    uₖ(x) = -ℯ.^(-γₖ .* x) ./ γₖ
    ==========================================================#
    # Compute the (nIs,K,Sim) matrices of v-NM utility values
    # Create one matrix per year-plan
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
        # Usim["$plan$year"] = -exp.(BigFloat.(-γₖ .* x)) ./ γₖ
        # println("$plan $year ", mean(Usim["$plan$year"]), " ", mean(x), " ", mean(γₖ))
    end

    
    #==========================================================
    PART B
    Take the matrices produced in the last subsection, and create
    three expected utility matrices (one for each year ) with dimen-
    sions (nIs,nPlans,Sim). These matrices should give the expected
    utility over the K simulated health draws for each person, plan,
    and simulation of random coefficients. Once computed, normalize
    these matrices by dividing by all entries by the expected utility
    for PPO1200 in each year and then taking the reciprocals of these
    normalized values. Taking the reciprocals is necessary because
    the CARA utility function is negative, so when you normalize,
    making all values positive, you have to take the reciprocal. The
    result is that expected utilities for PPO1200 will be 1 and those
    of the other plans will be defined in relation to that.
    ==========================================================#
    # For each year-plan utilty matrix, collapse the K matricies by
    # taking the mean -- resulting in one year-plan nIs x Sim matrix
    # Stack all year-plan matricies together for the same year,
    # resulting in three nIs x nPlans x Sim matricies
    Umean = Dict()
    for year ∈ 1:3
        # Preallocate an nIs x nPlans x Sim matrix with zeros (faster)
        Umean["$year"] = zeros(nIs, Sim, nPlans)
        # Umean["$year"] = zeros(BigFloat, nIs, Sim, nPlans)
        for (i, plan) ∈ enumerate(plans)
            # Fill in each of the nIs x Sim matrices with the mean of the K health simulations
            # Each element of the below nIs x Sim matrix is the mean of the K health simulation values 
            # Fill each of the nPlans (3) matrices in Umean["$year"] with the K-mean for that plan-year
            # End up with an nIs x nPlans x Sim tensor, a stack of nPlans matrices
            Umean["$year"][:, :, i] = mean(Usim["$plan$year"], dims=3)
        end
    end

    # Normalize by EU of PPO1200 for each year
    for year ∈ 1:3, i ∈ 1:3
        Umean["$year"][:, :, i] = Umean["$year"][:, :, 3] ./ Umean["$year"][:, :, i]
    end

    # for year ∈ 1:3; println("Umean $year: ", minimum(Umean["$year"]), " ", mean(Umean["$year"]), " ", maximum(Umean["$year"])); end

    # for year ∈ 1:3
    #     # Make sure that all elements of the bottom stack are 1 (for PPO1200)
    #     @assert  sum((Umean["$year"][:, :, 3] .≈ 1) .- 1) == 0
    # end
    

    
    #==========================================================
    PART C
    Create a ‘smoothed’ matrix of expected utilities that first raises
    the expected utilities you just calculated to some even power (say
    6) and second divides these transformed utilities for each family
    and each simulation in each year by the sum of their transformed
    utilities. Thus, for each year, simulated draw of preferences, and
    family, the numbers in your matrix should summed over the sec-
    ond dimension (plans) should add up to one. This ’smoothed’ ma-
    trix creates a continuous probabilistic function where the choice
    of the plan yielding highest expected utility approaches one and is
    increasing as the simulated utility gap between that and the other
    plans becomes larger. Your resulting matrices can thus be seen as
    a discrete probability distribution over all three choices for each
    simulated preference draw, person, and year. See Estimation Ap-
    pendix B for more details.
    ==========================================================#
    Usum = Dict()
    # Raise to power of 6
    for year ∈ 1:3
        s = sum(Umean["$year"].^6, dims=3)
        Usum["$year"] = Umean["$year"].^6 ./ s
    end
    # Make sure that the three plan values add up to 1 for all families
    # @assert sum((sum(Usum["3"], dims=3) .≈ 1) .- 1) == 0

    
    #==========================================================
    PART D
    Finally, it’s time to compute the log-likelihood function value,
    to do this you have to generate the probability for a given se-
    quence of three observed choices given the simulated choices for
    the candidate parameters. It is crucial to compute probabilities
    over the entire sequence of choices since inertia links choices over
    time for a given preference parameter simulation. This is the part
    where you match simulated choices to actual observed data. You
    should:
    - For each person, year, and simulated preference draw, multi-
      ply the probability of having made the observed year 1 choice
      by the observed year 2 choice by the observed year 3 choice.
      NOTE: In the limit, as the smoothing factor becomes large,
      this equals 1 for the sequence that contingently maximizes
      expected utility through the years, and 0 otherwise.
    - Find the average of the probabilities of observing the actual
      sequence of choices across all simulated draws (so along the
      Sim dimension). This should leave you with a (nIs,1) ma-
      trix that describes the probability of a given family making
      a given sequence of choices over time if the candidate param-
      eters, including those governing the distributions of random
      coefficients, are the true parameters.
    - Create the final log-likelihood function value equal to the
      sum of the log of these family choice sequence probabilities
      across all families. This is the log-likelihood function value.
      Make the output of the function the negative of this value,
      because non-linear optimization will search for the minimum
      of the likelihood function rather than the maximum (this
      convention is standard).
    ==========================================================#
    # Use the observed choice to select and multiply the probabilities
    probs = ones(nIs, Sim)
    # probs = ones(BigFloat, nIs, Sim)
    for year ∈ 1:3
        yearprobs = zeros(nIs, Sim)
        # yearprobs = zeros(BigFloat, nIs, Sim)
        # 
        for plan ∈ 1:3
            # Indicator for if the family chose this plan-year (row is 1 if so, 0 otherwise)
            yearplanindicator = repeat(Int.(choice[:, year] .== plan), outer=(1, Sim))
            # Each of the three stacked plan matricies will have zero rows in them,
            # and will not overlap non-zero rows with other plan matrices
            # so element-wise summing them will just squish all the non-zero rows
            # into the same nIs x Sim matrix, which is the correct rows of probability
            # that the plan choices were observed in this year
            yearprobs .+= Usum["$year"][:,:,plan] .* yearplanindicator
        end
        # Now we just need to multiply all the probability-of-observed-plan-choices together
        # This is an nIs x Sim matrix, where each row is the simulated joint probability that
        # the family would choose the plans they did in all three years, for each risk preference
        # simulation.
        probs .*= yearprobs
    end

    # Take the average over all simulations (along the Sim dimension)
    meanprobs = mean(probs, dims=2)
    
    # return the log-likelihood
    value = -sum(log.(meanprobs))
    println("Function Value: $(round(value))")
    return value
end





























#==============================================================================
QUESTION 4: OPTIMIZE OVER THE LIKELIHOOD

Now, you have your likelihood function all set up, which is the hard
part. To estimate the model parameters, all you have to do is run
a non-linear optimizer to find the parameter values that yield that
highest likelihood function value. This question takes you through
that process:
==============================================================================#
# Try adding manual large likelihood values for parameters outside the bounds
function nLL_bound(α)
    lower_violation = sum(α .< lb) > 0
    upper_violation = sum(α .> ub) > 0
    if lower_violation || upper_violation
        # Return large value so the algorithm keeps searching inside the bounds
        return Inf
    else
        return nLL(α)
    end
end

# Create an nLL wrapper to pass individual parameter elements to
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

function nLL_vec(a::Vector)
    return nLL_wrapper(a...)
end

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

function transform_to_21vec(x::AbstractMatrix)
    y = reshape(copy(x), 24, 1)
    # remove the (6,2), (6,3), and (6,4) elements
    y = deleteat!(vec(y), [12, 18, 24])
    return y
end



#==========================================================
PART A
You will do constrained minimization with either fmincon or ktr-
link from the KNITRO optimization package. To do this you
need to define the starting values and bounds for parameter ma-
trix alpha which goes into the likelihood function. Define your
bounds as follows (the bounds for the parameters alpha must be
the same dimension as alpha):
lb = [0, −.5, −.5, 0; −Inf, 0.1, −Inf, 0.1; 1, 1, 1, 1; −Inf, 0, 0, −5000; -5000, −5000, −5000, −5000; −5000, 0, 0, 0; ];
ub = [2, 1, 1, 1; Inf, Inf, Inf, Inf ; 20000, 20000, 20000, 20000; Inf, 20000, 20000, 5000; 5000, 5000, 5000, 5000; 5000, 0, 0, 0; ];

Also, you have to define starting values for the parameter matrix 
alpha. Use the following matrix to begin with for this:
alpha0 = [0.06, 0.003, 0, 0.04; −2500, 700, −2200, 800; 300, 800, 300, 800; -500, 1250, 1750, −500; 0, 0, 0, 0; 0, 0, 0, 0; ];
==========================================================#
lb = [
          0.0     -0.5     -0.5      0.0
        -Inf       0.1    -Inf       0.1
          1.0      1.0      1.0      1.0
        -Inf       0.0      0.0  -5000.0
      -5000.0  -5000.0  -5000.0  -5000.0
      -5000.0      0.0      0.0      0.0
]
lb_vec = transform_to_21vec(lb)

α₀ = [
        0.06     0.003      0.0     0.04
    -2500.0    700.0    -2200.0   800.0
      300.0    800.0      300.0   800.0
     -500.0   1250.0     1750.0  -500.0
        0.0      0.0        0.0     0.0
        0.0      0.0        0.0     0.0
]
α₀_vec = transform_to_21vec(α₀)


ub = [
        2.0      1.0      1.0      1.0
        Inf      Inf      Inf      Inf
    20000.0  20000.0  20000.0  20000.0
        Inf   20000.0  20000.0   5000.0
    5000.0   5000.0   5000.0   5000.0
    5000.0      0.0      0.0      0.0
]
ub_vec = transform_to_21vec(ub)

α₂ = [
    0.000232     2.9e-5      2.27e-6     0.000188
    -2912.0        843.0     -2871.0       897.0
      204.0        502.0       329.0       811.0
      856.0       2480.0      1729.0      -551.0
      -32.0          5.0       198.0        80.0
      156.0          0.0         0.0         0.0
]
α₂_vec = transform_to_21vec(α₂)

function stack_cols(x)
    rows, cols = size(x)
end

α₃ = reshape(α₂, 24, 1)[1:21]
reshape([α₃; [0,0,0]], 6, 4)




#==========================================================
PART B
Before beginning non-linear optimization compute your likelihood
function value at alpha0. What value do you get? If your
value is in the neighborhood of 2000-3000, you are in good shape
(the actual value is the negative of this, but the output will be
positive since you have defined your likelihood output as the neg-
ative of the true output). Make sure when evaluate your likelihood
function that you pass in not only alpha0, but also all fixed vari-
ables that go into that function, including demographics, the nor-
mal draws and objects from 1c, and the matrices from 1d which
help speed things up. Interpret this value of the likelihood
function, what does it mean?
==========================================================#
@time nLL(α₀)  # 6116.314354461145

@time nLL(α₂)  # 7428.797342352567

a2 = [
    0.304718      0.927983      0.514715      0.599309
    -4540.06       3875.96       8052.78       1024.7
     5074.24       3919.88        794.49      11247.8
     6998.8       14282.5       18565.5        1949.27
    -1942.63      -3566.87      -4394.86       2534.03
     1550.16          0.0           0.0           0.0
]
@time nLL(a2)


alpha0 = [
    0.06 0.003 0 0.04
     -2500 700 -2200 800
    300 800 300 800
    -500 1250 1750 -500
    0 0 0 0
    0 0 0 0
]

# Starting values also run for this simulated data code  to find best likelihood value. Actual alpha0 used here above was one with highest likelihood value. 

alpha1 = [
    0.1 0.04 0 0.03
     900 0.1 1500 500
    800 300 800 300
    -2500 700 1500 -200
    20 100 100 -100
    -100 0 0 0
]
alpha2 = [
    0.3 0.1 -0.0003 0.08
     600 0.1 2500 300
    500 200 500 200
    -1500 700 1200 -800
    200 0 0 -200
    -100 0 0 0
]
alpha3 = [
    0.5 0.02 -0.0003 0.03
     1500 0.1 1700 1
    400 150 400 150
    -1000 500 1500 -500
    0 50 200 -100
    100 0 0 0
]
alpha4 = [
    0.9 0.05 0.003 0.06
     1000 0.1 1500 100
    150 150 150 150
    -2000 800 1050 -700
    -50 -50 100 100
    -50 0 0 0
]































#==========================================================
PART C
Set up your non-linear optimization routine with fmincon or ktr-
link. This should take in starting value alpha0, but optimize over
alpha starting from that point to find the alpha that minimizes
the negative likelihood function (maximizes the true function)
subject to the bounds on alpha passed in. If you have trouble
with this part, after struggling for a while you can email
me for some hints. You will want to setup some options to be
passed into the optimization such as:

options = optimset( 0 M axF unEvals 0 , 40000, 0 M axIter 0 , 40000, 0 T olF un 0 , 10 − 4);

What do these options do for the optimization? How does vary-
ing them impact the results if at all? What other options could
be used and how would they impact estimation?

Run your non-linear optimization, and report out the values for
all 21 parameters in alpha that you’ve estimated (alpha(6,2), al-
pha(6,3) and alpha(6,4) are zeros filling out the matrix.). Tell me
what you get for:
    -Risk preference estimates
    -Inertia estimates
    -PPO1200 random coefficient estimates
    -Epsilon variance estimates
    -Estimated preference of high total cost people for PPO250
    -The final likelihood function value at the estimated parame-
        ters (i.e. the ’maximum’ likelihood function value)
==========================================================#

#= ####################################################### 
        BASIC JULIA OPTIM PACKAGE
- The function is not auto-differentiable, so needed to use
    a "derivative-free" algorithm with parameter bounds.
- Two algorithms in the optim.jl package works for this:
    1. NelderMead() is the default
    2. SAMIN() Simulated Annealing with bounds

After some testing, the Simulated Annealing algorithm seemed
to preform better -- faster iterations and lower objective
function values for the same number of iterations. This does
not mean that it's more likely than the NelderMead algorithm
to converge on the global minimum, but for the purposes of 
testing different starting parameter values with a relatively
small number of iterations (compared to what I will be
using with the Knitro package), Simulated Annealing can give
me a better quick idea of how starting parameter values
affect the estimated minimizing parameter values.

Below I show the results of minimizing the objective function
using the SAMIN algorithm with different starting parameters,
limiting to 500 iterations for each run. I also compare these
results to the Knitro package, run for X hours, from the first 
suggested starting values. Just due to lack of time, I end the
estimation there, with the understanding that for real research,
these may need to run for longer on a better computer with higher
precision data types (see below section).

There also exists the BlackBoxOptim package, which might be
better, but it does not work well when much of the domain
of the function returns NaN. The loglikelihood function 
returns NaN at some parameter
values that are far from the starting point alpha0. This
happens because some of the X values are large in the 
utility function, so when calculating e^(-γ*x), this 
evaluates to zero when using the Float64 data type.
The obvious solution would be to use a more precise data type,
and BigFloat is the next largest (without creating my own 
custom datatype) and has arbitrary precision. But when the 
numbers are very small (on the order of 10^-2000), this can
take up a large amount of memory. As a result, each evaluation
of the log-likelihood function takes a very long time to run
when using BigFloats.

This might be the best solution
to fully explore the parameter space, but sadly my computer
doesn't seem to have enough RAM to run this with BigFloats,
and it would seem to take an extremely long time. This could
be parallelized with the Distributed library, but I decided
to focus on other methods that would get me quicker results.
I am working on a maximum likelihood project in Julia over
the summer, however, my summer
project lends well to taking derivatives, so I will likely
stick to gradient decent and other more optimized packages.
I am keeping these methods in mind for more complex likelihood
functions that don't have easy gradients to derive (or they
don't exist ex-ante). If I needed to use BigFloats for a 
research project like this, I would probably spin-up an
AWS elastic compute instance to run something like this
in parallel, where the RAM and CPU can scale easily.
=# #######################################################

#! Optimize on the full 6x4 α matrix
resNM = Optim.optimize(nLL, lb, ub, α₀, NelderMead(), Optim.Options(iterations=50))
Optim.minimizer(resNM)
Optim.minimum(resNM)
Optim.iterations(resNM)
#= This results in some of the three zeros being non-zero,
which doesn't really matter, but wanted to use a 21-vector
version instead. =#


#! Compare between NelderMead and SAMIN

#= NELDERMEAD RESULT AFTER 50 ITERATIONS
minimum(resNM) = 5438.954383339348
    Seconds run:   967  (vs limit Inf)
    Iterations:    50
    f(x) calls:    567

SAMIN RESULT AFTER 50 ITERATIONS
minimum(resSM) = 4130.821817817749
    Seconds run:   85  (vs limit Inf)
    Iterations:    50
    f(x) calls:    50

I think the NelderMead run used the 6x4 matrix, so it may have taken longer
just because it had more parameters to change, even though those parameters
did not actually affect the function value. But the difference in the
likelihood value after 50 iterations made me think SAMIN would be good for
this short test. I did other longer comparisons, but didn't save the results
but in general, SAMIN minimized the function much faster than NelderMead.
I don't know if this is generalizable, but might be worth comparing
NelderMead, SAMIN, and Knitro over many iterations for any serious research project.
Perhaps running all in parallel. And need to search many points. 
=#

# Time limit for each of 5 iterations to finish in 8 hours (in seconds)
time_limit5 = round(8*3600/5)
max_iter = 40_000

#! Save results from 5 different starting parameter vectors using SAMIN
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

# Convert strings saved to vectors
# reordering index vector to match order in paper
 rr = [15, 10, 20, 5, 11, 16, 21, 6, 1, 7, 12, 17, 2, 8, 13, 18, 4, 3, 9, 14, 19]

# ax: final param values; fx: final function value; iterx: function iterations; hourx: hours run
alpha0
a0 = transform_to_21vec([0.09688801258129792 0.0031723346939378973 -0.00011720381932389255 0.0231123314570683; -2500.0 700.0 -2200.0 800.0; 1263.0280770307288 32.60199522848822 11822.503819486612 5843.938487718894; -500.0 9464.105148578532 19796.688827386457 2744.4285755036844; 4566.83601823626 2618.8911679089965 4945.547619297684 4599.025722723732; 4252.686338748319 0.0 0.0 0.0])
f0 = 2.636075e+03
iter0 = 3683
hour0 = 8/5
col0 = [a0[rr]; [f0]; [iter0]; [hour0]]

alpha1
a1 = transform_to_21vec([0.33002667495341054 -0.020415136894857633 -0.4506208593746673 0.07666154611931172; 600.0 0.1 2500.0 300.0; 10088.418769865137 5327.0770765116085 5641.200606354002 1470.1555917804735; -1500.0 16059.153059978 17954.55550591196 800.6473700094812; 4149.55130858559 -740.6309433652841 3230.3487060441985 730.9802312789325; 3658.9380548488753 0.0 0.0 0.0])
f1 = 7.006379e+03
iter1 = 3804
hour1 = 8/5
col1 = [a1[rr]; [f1]; [iter1]; [hour1]]

alpha2
a2 = transform_to_21vec([0.33002667495341054 -0.020415136894857633 -0.4506208593746673 0.07666154611931172; 600.0 0.1 2500.0 300.0; 10088.418769865137 5327.0770765116085 5641.200606354002 1470.1555917804735; -1500.0 16059.153059978 17954.55550591196 800.6473700094812; 4149.55130858559 -740.6309433652841 3230.3487060441985 730.9802312789325; 3658.9380548488753 0.0 0.0 0.0])
f2 = 7.006379e+03
iter2 = 3804
hour2 = 8/5
col2 = [a2[rr]; [f2]; [iter2]; [hour2]]

alpha3
a3 = transform_to_21vec([0.07077629277866049 0.007360715245031127 0.00024211588307838133 0.04606341465723755; 1500.0 0.1 1700.0 1.0; 803.4826373995114 143.5630946235941 9511.086878827799 4080.024004492574; -1000.0 15961.01258553203 17347.8532119649 2008.01646713804; 3260.3224609100957 -1108.4887563414381 3793.489435633259 1915.9255204139442; 2102.074815042707 0.0 0.0 0.0])
f3 = 2.749408e+03
iter3 = 3794
hour3 = 8/5
col3 = [a3[rr]; [f3]; [iter3]; [hour3]]

alpha4
a4 = transform_to_21vec([0.9 0.05 0.003 0.06; 1000.0 0.1 1500.0 100.0; 150.0 150.0 150.0 150.0; -2000.0 800.0 1050.0 -700.0; -50.0 -50.0 100.0 100.0; -50.0 0.0 0.0 0.0])
f4 = NaN
iter4 = 3709
hour4 = 8/5
col4 = [a4[rr]; [f4]; [iter4]; [hour4]]






#= ####################################################### 
        KNITRO PACKAGE
The Knitro package has changed a little bit since the 2013 paper
but I was able to get it working in Julia. Sadly, the Julia
wrapper has some limitations -- it cannot parallelize (there
is an open GitHub issue about this, and needs to be resolved
by the folks at Artelys / Knitro). This means that I could
not run the multi-start function efficiently, which would
have allowed me to explore effects of starting parameters
much better. So I chose to just use one run of knitro, and 
compare it to other runs of the Julia Optim package with
the other starting parameter values. Knitro seems to take
longer but more reliably minimizes the objective function.
(After just a few hours, it got to a lower value than
leaving the NelderMead alogrithm to run for 7.1 hours). 





Final Statistics
----------------
Final objective value               =   2.57664797201973e+03
Final feasibility error (abs / rel) =   0.00e+00 / 0.00e+00
Final optimality error  (abs / rel) =   6.66e-04 / 6.66e-04
# of iterations                     =        437 
# of CG iterations                  =         64 
# of function evaluations           =      11643
# of gradient evaluations           =          0
Total program time (secs)           =   25651.86133 ( 25372.193 CPU time)
Time spent in evaluations (secs)    =   25585.15820

Solution Vector
---------------
x[       0] =   5.81624945391e-03,   lambda[       0] =  -3.63790430279e-02
x[       1] =  -2.50000000000e+03,   lambda[       1] =   0.00000000000e+00
x[       2] =   1.87592278071e+02,   lambda[       2] =   2.49507016724e-09
x[       3] =  -5.00034030392e+02,   lambda[       3] =   0.00000000000e+00
x[       4] =   4.99999999058e+03,   lambda[       4] =   5.64592286220e-03
x[       5] =   1.29503272026e+03,   lambda[       5] =   1.71430243683e-07
x[       6] =  -4.76759555815e-04,   lambda[       6] =  -3.26926845896e-01
x[       7] =   7.01983506841e+02,   lambda[       7] =  -0.00000000000e+00
x[       8] =   1.00000257523e+00,   lambda[       8] =  -2.41087625472e-02
x[       9] =   1.26422116323e+04,   lambda[       9] =   1.47032612412e-07
x[      10] =  -2.42213712488e+03,   lambda[      10] =   3.32056214322e-08
x[      11] =  -1.25474724273e-05,   lambda[      11] =  -5.42619286981e+00
x[      12] =  -2.20000420498e+03,   lambda[      12] =   0.00000000000e+00
x[      13] =   7.65718286314e+03,   lambda[      13] =   6.34971662095e-09
x[      14] =   1.99999999789e+04,   lambda[      14] =   1.89841461382e-03
x[      15] =   4.99999947757e+03,   lambda[      15] =   4.59740794089e-04
x[      16] =   2.93511072139e-05,   lambda[      16] =   2.48833069306e-01
x[      17] =   8.01689955868e+02,   lambda[      17] =  -7.06225156035e-11
x[      18] =   1.00000330348e+00,   lambda[      18] =  -2.66111523059e-02
x[      19] =   4.63907933581e+03,   lambda[      19] =   1.12066338522e-06
x[      20] =   4.99999992143e+03,   lambda[      20] =   6.49951182089e-04
===============================================================================
=# #######################################################
include(string(dir,"knitro-applied.jl"))

# Convert above string results to vector of parameters
x=Dict()
x[ 0] =   5.81624945391e-03
x[ 1] =  -2.50000000000e+03
x[ 2] =   1.87592278071e+02
x[ 3] =  -5.00034030392e+02
x[ 4] =   4.99999999058e+03
x[ 5] =   1.29503272026e+03
x[ 6] =  -4.76759555815e-04
x[ 7] =   7.01983506841e+02
x[ 8] =   1.00000257523e+00
x[ 9] =   1.26422116323e+04
x[10] =  -2.42213712488e+03
x[11] =  -1.25474724273e-05
x[12] =  -2.20000420498e+03
x[13] =   7.65718286314e+03
x[14] =   1.99999999789e+04
x[15] =   4.99999947757e+03
x[16] =   2.93511072139e-05
x[17] =   8.01689955868e+02
x[18] =   1.00000330348e+00
x[19] =   4.63907933581e+03
x[20] =   4.99999992143e+03
# minimizing parameter vector
a5 = [x[i] for i in 0:20]
# minimized function value
f5 = 2.57664797201973e+03
# # of function evaluations
iter5 = 11643
# hours spent
hour5 = 25651.86133 / 3600
# Put into one column
col5 = [a5[rr]; [f5]; [iter5]; [hour5]]


# Table made by hand in ps1-tables.ods
# then copy-pasta'd to tablesgenerator.com for formatting














#=====================
Try using ConstrainedOptim.jl
=====================#

fun(x) =  (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
x0 = [0.0, 0.0]
df = TwiceDifferentiable(fun, x0)
lx = [-0.5, -0.5]; ux = [1.0, 1.0]
dfc = TwiceDifferentiableConstraints(lx, ux)
res = Optim.optimize(df, dfc, x0, IPNewton())
res = Optim.optimize(fun, lx, ux, x0, NelderMead())
resNM = Optim.optimize(nLL, lb, ub, α₀, NelderMead(), Optim.Options(iterations=50))
resSM = Optim.optimize(nLL_vec, transform_to_21vec(lb), transform_to_21vec(ub), transform_to_21vec(α₀), SAMIN(), Optim.Options(iterations=50))


#= BlackBoxOptim does not seem to work. 
All function evaluations return NaN, 
seems that all -Inf values convert to Nan, but Inf values are fine.
Could not find any documentation on this.
Changed the Inf values to a large number.
A different issue: nLL() seems to return NaN for some α in the domain.
Should track down what is blowing up... probably the ℯ^γ utility function?
=#
using BlackBoxOptim
infreplace = 10000
searchrange = replace.(collect(zip(lb_vec, ub_vec)), -Inf => -infreplace, Inf => infreplace)

resbb = bboptimize(nLL_vec; SearchRange = searchrange, Method = :random_search, MaxTime = 10.0)

function rosenbrock2d(x)
    return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end
res = bboptimize(rosenbrock2d; SearchRange = (-5.0, 5.0), NumDimensions = 2)

#=====================
Parallelize?
=====================#
# First run without any parallel procs used in eval
opt1 = bbsetup(slow_rosenbrock; Method=:xnes, SearchRange = (-5.0, 5.0),
               NumDimensions = 50, MaxFuncEvals = 5000)
el1 = @elapsed res1 = bboptimize(opt1)
t1 = round(el1, digits=3)

# When Workers= option is given, BlackBoxOptim enables parallel
# evaluation of fitness using the specified worker processes
opt2 = bbsetup(slow_rosenbrock; Method=:xnes, SearchRange = (-5.0, 5.0),
               NumDimensions = 50, MaxFuncEvals = 5000, Workers = workers())
el2 = @elapsed res2 = bboptimize(opt2)
t2 = round(el2, digits=3)





# Example 2
f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
lower = [1.25, -2.1]
upper = [Inf, Inf]
initial_x = [2.0, 2.0]
inner_optimizer = GradientDescent()
results = optimize(f, lower, upper, initial_x, Fminbox(inner_optimizer))




############################################# Runs but returns the same initial matrix
# Try with our example
inner_optimizer = GradientDescent()
results = optimize(nLL, lb, ub, α₀, Fminbox(inner_optimizer))
results2 = optimize(nLL, lb, ub, α₂, Fminbox(inner_optimizer))
# Doesn't work because the gradient is evaluated as NaN (i think)




# Custom printing function
start_time = time()
time_to_setup = zeros(1)
function advanced_time_control(x)
    println(" * Iteration:       ", x.iteration)
    so_far =  time()-start_time
    println(" * Time so far:     ", so_far)
    if x.iteration == 0
        time_to_setup[:] = time()-start_time
    else
        expected_next_time = so_far + (time()-start_time-time_to_setup[1])/(x.iteration)
        println(" * Next iteration ≈ ", expected_next_time)
        println()
        # return expected_next_time < 13 ? false : true
    end
    println()
    false
end
results3 = optimize(nLL_vec, transform_to_21vec(lb), transform_to_21vec(ub), transform_to_21vec(α₂), Fminbox(inner_optimizer), Optim.Options(callback = advanced_time_control))
############################################












#########################################











#=========
Try using KNITRO with JuMP
JuMP is a nice modeling syntax
But I was unable to get this to work and instead
just used the base-knitro julia handler
see knitro-applied.jl
=========#
# Example 1
model = Model(KNITRO.Optimizer)

initval = [-2., 1.]
@variable(model, x[i=1:2], start=initval[i])
@NLobjective(model, Min, 100*(x[1] - x[2]^2)^2 + (1 - x[1])^2)
@constraint(model, x[1] <= 0.5)
c1 = @constraint(model, x[1] * x[2] >= 1)
c1 = @constraint(model, x[1] + x[2]^2 >= 0)

a = JuMP.optimize!(model)

using MathOptInterface
const MOI = MathOptInterface
MOI.get(model, MOI.VariablePrimal(), x)
get(model, VariablePrimal(), x)


# Example 2
m = Model(optimizer_with_attributes(KNITRO.Optimizer,
                                    "honorbnds" => 1, "outlev" => 1, "algorithm" => 4)) # (1)
@variable(m, x, start = 1.2) # (2)
@variable(m, y)
@variable(m, z)
@variable(m, 4.0 <= u <= 4.0) # (3)

mysquare(x) = x^2
register(m, :mysquare, 1, mysquare, autodiff = true) # (4)

@NLobjective(m, Min, mysquare(1 - x) + 100 * (y - x^2)^2 + u)
@constraint(m, z == x + y)

optimize!(m)
(value(x), value(y), value(z), value(u), objective_value(m), termination_status(m)) # (5)





# Apply to our problem
model = Model(KNITRO.Optimizer)

register(model, :nLL_wrapper, 21, nLL_wrapper, autodiff=false)
@variable(model, a[i=1:21], start=transform_to_21vec(α₀)[i])
@NLobjective(model, Min, nLL_wrapper(a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],a[20],a[21]))
@constraint(model, transform_to_21vec(lb) .<= a .<= transform_to_21vec(ub))

JuMP.optimize!(model)
solution_summary(model)




#=========
Try using KNITRO without JuMP
=========#
#! try to find the replication package from the published paper
# see knitro-applied.jl








#=========
Try using KNITRO with NLopt for black-box and derivative free
=========#
nLL_nlopt(x::Vector, grad::Vector) = nLL_vec(x)
using NLopt, Optim
opt = Opt(:LD_MMA, 21)
opt.min_objective = nLL_nlopt
opt.lower_bounds = transform_to_21vec(lb)
opt.upper_bounds = transform_to_21vec(ub)
opt.maxeval = 1000 # Increase this eventually
(optf,optx,ret) = NLopt.optimize(opt, transform_to_21vec(α₀))
# The above code just uses NLopt, not knitro yet. Was just trying to get NLopt to run
















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





#==========================================================
PART F
Finally, re-run the optimization with a few different starting val-
ues that you choose. Do the parameter estimates change? If so,
which ones are ’better’ i.e. which ones would you select to report
as estimates? Why?
==========================================================#
























#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#!                      NOTES SECTION
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#==========================================================
OTHER OPTIMIZATION FUNCTIONS using derivatives

See below link for Optim solvers and use cases.
https://julianlsolvers.github.io/Optim.jl/stable/#user/config/#solver-options

For a low-dimensional problem with analytic gradients and 
Hessians, use the Newton method with trust region. For larger 
problems or when there is no analytic Hessian, use LBFGS, 
and tweak the parameter m if needed. If the function is 
non-differentiable, use Nelder-Mead. Use the HagerZhang 
linesearch for robustness and BackTracking for speed.

Also, these methods could be good for large search domains
when there is no way to autdiff the gradient,
but all the usual features are not fully implented:
- Simulated Annealing
- Simulated Annealing with bounds SAMIN
- Particle Swarm
==========================================================#

# Create a gradient (and hessian?)
func = TwiceDifferentiable(nLL, α₀; autodiff=:forward);
opt = optimize(func, α₀)


# Nelder-Mead, no Gradient
opt = optimize(nLL, α₀)
# Did not finish after >24 hours. Uses NelderMead, 
# but because I never gave it stopping time/iterations, it just keeps
# going until it coverges
results = optimize(f, lower, upper, initial_x, NelderMead())

# No gradient? Maybe it caluculates a gradient
Optim.minimizer(optimize(nLL, α₀, BFGS()))

# No gradient 3 (tries to find analytical-ish gradient using autodiff, but doesn't work for this function)
Optim.minimizer(optimize(nLL, α₀, BFGS(); autodiff = :forward))

# If I had a gradient g!
# Add upper and lower bounds
inner_optimizer = GradientDescent()
results = optimize(f, g!, lower, upper, initial_x, Fminbox(inner_optimizer))

# Limited memory -- can often converge faster, but may be less accurate.
# see 
optimize(f, x0, LBFGS(); autodiff = :forward)


# Could try Simulated Annealing with bounds (SAMIN) https://julianlsolvers.github.io/Optim.jl/stable/#algo/samin/

# Could try Particle Swarm https://julianlsolvers.github.io/Optim.jl/stable/#algo/particle_swarm/

#! Which methods use a numerical gradient?





#==========================================================
NOTE ON KNITRO SOLVER INSTALLATION
1. Go to www.artelys.com/solvers/knitro/ to download a trial
    license. You need to create an account. They will send you
    an activation email, mine was kinda buggy. When I clicked
    the link in the email, it said error, but when I went back
    to the website, I was able to login.
2. After your account is activated, go back to www.artelys.com/solvers/knitro/
    and click Download a Trial License. Login, and click
    the link to download Knitro (of whatever the latest version is)
3. This will take you to a page with download links on the left
    and a form on the right. If it's in a different language, click
    the "en" at the very top right of the screen to change to
    English. Then click Installer or Archive for the Artelys Knitro 13.0
    (or whatever version is latest). This will download an .exe or .tar.gz
    file. 
4. On the same page, click the fill out the form on the right that says
    My License. They will email you the license .txt file.
5. Copy the license file into one of the locations specified in
    the INSTALL file. I put it in my $HOME directory.
6. If you downloaded the .tar.gz file, decompress it and put the entire
    "knitro-XX.X.X..." folder into whatever directory you like. (like a
    programs directory or something. I put it in $HOME).
7. Now we need to tell Julia where the knitro folder is, then
    install knitro:
- open a juila prompt or wherever you are executing julia
    scripts. Add the knitro folder location to the environment
    variables dictionary:
        `ENV["KNITRODIR"] = "path/to/knitro-XX.X.X-..."`
    where "path/to/knitro-XX.X.X-..." is your path. For example:
        `ENV["KNITRODIR"] = "/home/a/knitro-13.0.1-Linux-64"`
    was my path.
--> instead of doing this in each session, you can permently add
    KNITRODIR as an environment variable at startup by editing the
    bashrc file. I added this line to the end of ~/.bashrc:
        `export KNITRODIR="/home/a/knitro-13.0.1-Linux-64"`
- Add the knitro package using
    ```julia
    julia> ]
    (v1.7) pkg> add KNITRO
    ```
    If it says "Unable to locate KNITRO installation", you might need
    to restart your julia session if you were running things before.
    Also, this seems to have issues when running from a jupyter notebook.
    I suggest you do the installation process from a julia REPL.
- If that works, next build KNITRO using
    ```julia
    julia> ]
    (v1.7) pkg> build KNITRO
    ```
- Test this quickly by `using KNITRO`.
- If that runs without erorr, test more using KNITRO's test:
    ```julia
    julia> ]
    (v1.1) pkg> test KNITRO
    ```
    It will spit out a lot of info, but look for four Test Summary's. 
    You should see "Test C API", "Test examples",
    "Test MathOptInterface", and "Test C API License".
    For me, "Test MathOptInterface" took longer than the other tests,
    you might need to be patient. At the very end, you should see
    "Testing KNITRO tests passed", success!
- I suggest also installing JuMP, it makes the syntax a lot nicer:
    ```julia
    julia> ]
    (v1.7) pkg> add JuMP
    ```
    
- Now it's time to test one of the examples. 

```julia 
using JuMP, KNITRO
model = Model(with_optimizer(KNITRO.Optimizer))

initval = [-2., 1.]
@variable(model, x[i=1:2], start=initval[i])
@NLobjective(model, Min, 100*(x[1] - x[2]^2)^2 + (1 - x[1])^2)
@constraint(model, x[1] <= 0.5)
c1 = @constraint(model, x[1] * x[2] >= 1)
c1 = @constraint(model, x[1] + x[2]^2 >= 0)

JuMP.optimize!(model)
```
==========================================================#







