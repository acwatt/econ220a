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
# Optimizing functions
using Optim, NLSolversBase
using LinearAlgebra: diag



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


#==========================================================
PART C & D: TRANSFORM DATA AND VARIABLE CREATION
==========================================================#
# Run the code in ps1_transformdata.jl to create objects in memory
# to use during the simulations / estimation
include(string(dir,"ps1_transformdata.jl"))





























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
    return elmax(dist, 0.000001)
end


function simulate_ppo1200_dist(; αm_sgl, αsd_sgl, αm_fam, αsd_fam, rand_coef)
    # IND' = vector of 1=single, 0=family. Family Status
    # For each row (family), need to create mean + sd*rand_coef, depending on Family Status
    dist = IND' .* (αm_sgl .+ αsd_sgl*rand_coef) .+ (1 .- IND') .* (αm_fam .+ αsd_fam*rand_coef)  # dim = nIs x Sim
    dist = repeat(dist, outer=(1, 1, K))  # dim = nIs x Sim x K
    # maximum of matrix and 0.000001, truncating coefficients above 0
    return elmax(dist, 0.000001)
end


"""Create ​ ϵₖⱼₜ with sd = σ​²_ϵⱼ(​Yₖ) (vector of family-plan-time specific shocks)
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
        Col 1: intercept of the mean of the normal distribution of risk preferences  α[1,1]=2.23e-4
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
    α = [2.23e-4 2.90e-5 2.27e-6 1.88e-4
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
                              rand_coef = randcoef_risk)
    
    
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
    for year ∈ 1:3
        yearprobs = zeros(nIs, Sim)
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
    return -sum(log.(meanprobs))
end





























#==============================================================================
QUESTION 4: OPTIMIZE OVER THE LIKELIHOOD

Now, you have your likelihood function all set up, which is the hard
part. To estimate the model parameters, all you have to do is run
a non-linear optimizer to find the parameter values that yield that
highest likelihood function value. This question takes you through
that process:
==============================================================================#



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

ub = [
        2.0      1.0      1.0      1.0
        Inf      Inf      Inf      Inf
    20000.0  20000.0  20000.0  20000.0
        Inf   20000.0  20000.0   5000.0
    5000.0   5000.0   5000.0   5000.0
    5000.0      0.0      0.0      0.0
     ]

α₀ = [
        0.06     0.003      0.0     0.04
    -2500.0    700.0    -2200.0   800.0
      300.0    800.0      300.0   800.0
     -500.0   1250.0     1750.0  -500.0
        0.0      0.0        0.0     0.0
        0.0      0.0        0.0     0.0
    ]






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
@time nLL(α₀)





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
# Create a gradient (and hessian?)
func = TwiceDifferentiable(nLL, α₀; autodiff=:forward);
opt = optimize(func, α₀)

# No gradient
opt = optimize(nLL, α₀)

# No gradient? 2
Optim.minimizer(optimize(nLL, α₀, BFGS()))

# No gradient 3
Optim.minimizer(optimize(nLL, α₀, BFGS(); autodiff = :forward))

# Add upper and lower bounds
inner_optimizer = GradientDescent()
results = optimize(f, g!, lower, upper, initial_x, Fminbox(inner_optimizer))

optimize(f, x0, LBFGS(); autodiff = :forward)




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
# Nelder-Mead, no Gradient
opt = optimize(nLL_bound, α₀)

# Could try Simulated Annealing with bounds (SAMIN) https://julianlsolvers.github.io/Optim.jl/stable/#algo/samin/

# Could try Particle Swarm https://julianlsolvers.github.io/Optim.jl/stable/#algo/particle_swarm/



#==========================================================
PART D
Incorporating demographic heterogeneity and the parameter es-
timates, what is the population distribution of inertial costs?
==========================================================#






#==========================================================
PART E
Consider the mean of the CARA risk preference distribution of
random coefficients. Translate this coefficient into the value of X
that makes a family indifferent between no gamble and a gamble
where they win $100 with 50% probability and lose $X with 50%
probability.
==========================================================#






#==========================================================
PART F
Finally, re-run the optimization with a few different starting val-
ues that you choose. Do the parameter estimates change? If so,
which ones are ’better’ i.e. which ones would you select to report
as estimates? Why?
==========================================================#


















