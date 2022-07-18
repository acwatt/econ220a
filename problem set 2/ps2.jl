#==============================================================================
file: ps2.jl
description: estimating model for Problem Set #2 of Econ 220A (UC Berkeley, 2022),
             Demand Estimation and Merger Simulation. The problem set file is
             . The two original files are:
             - PSET1-2019-2.pdf
             - ps1data.csv
author: Aaron C Watt (UCB Grad Student, Ag & Resource Econ)
notes:
    - setting up problems
    - Nevo citation: Aviv Nevo, “A Practitioner’s Guide to Estimation of Random-Coefficients Logit Models of Demand,” Journal of Economics & Management Strategy 9, no. 4 (Winter 2000): 513–48, https://doi.org/10.1162/105864000567954.
==============================================================================#

#?==============================================================================
#?                                PACKAGES
#?==============================================================================
# For Stats stuff
using Statistics, StatsBase, StatsPlots, StatsFuns, GLM
using CovarianceMatrices, StatsModels, FixedEffectModels
# For dataframes
using DataFrames, CSV, DataFramesMeta, Chain
# Random Variable distributions
using Random, Distributions
Random.seed!(0);
# Output to files
using Latexify, TikzPictures
# Optimizing functions
# using Optim, NLSolversBase
# using LinearAlgebra: diag
# using JuMP, KNITRO




#?=============================================================================
#?                                 SETUP
#?=============================================================================
# Load data
df = DataFrame(CSV.File(string(dirname(@__FILE__), "/ps1data.csv")))








#?=============================================================================
#?=============================================================================
#?=============================================================================
#?                        QUESTION 1: LOGIT
#?=============================================================================
#?=============================================================================
#?=============================================================================
############# PART A
# Estimate an aggregate Logit model using OLS (Product 3 is reference group)

# Add δ (mean utility levels)
df = add_logit_depvar!(df)
# Regress δ on X to estimate α,β that minimizes error term
# δ = mean "observed" utility, Xβ-αp = mean predicted utility (α,β are mean utility parameters)
reg1a = reg(df, @formula(δjn ~ pjn + d1 + d2  + d4 + d5 + x1 + x2), Vcov.cluster(:market))
# Predict the shares using Xβ-αp (Nevo Eq. 6)
df = add_predicted_shares!(reg1a, df)
# Calculate estimated own-price elasticities for each product 
df = add_ownprice_elasticities!(reg1a, df)
# Calculate estimated cross-price elasticity with respect to product 3.
df = add_crossprice_elasticities!(reg1a, df)
df[!, [:market, :sjn, :sjn_hat, :δjn, :ηjjn, :ηj3n, :prodid]]
# Calculate average (across markets) elasticities



#logit = glm(fm, train, Binomial(), ProbitLink())


#*#############################################################################
#*####################### QUESTION 1 FUNCTIONS ################################

"""Create the logit dependent variable: δⱼ,ₙ = log(sjn) - log(sj0);  sj0 = outside share (1-Σsjn)"""
function add_logit_depvar!(df)
    df.δjn = @chain df begin
        groupby(:market)
        combine(:sjn => (x -> 1-sum(x)) => :sj0)
        leftjoin(df, _, on = :market)
        log.(_[!,:sjn]) - log.(_[!,:sj0])
    end
    return df
end


eq6(x) = exp.(x) ./ (1 + sum(exp.(x)))

"""Estimates shares using linear paremeters (Nevo Eq 6) using regress results."""
function add_predicted_shares!(reg, df_)
    df = copy(df_)
    # Predict δ for all observations (xβ-αp)
    df.δ = predict(reg1a, df)
    # Estimate share using Eq (6)
    df_.sjn_hat = @chain df begin
        groupby(:market)
        combine(:δ => eq6 => :sjn_hat)
        _[!,:sjn_hat]
    end
    return df_
end


"""Estimate own-price elasticities of the predicted market shares"""
function add_ownprice_elasticities!(reg, df)
    # Get price coefficient
    αidx = findfirst(==("pjn"), reg.coefnames)
    α = -reg.coef[αidx]
    df.ηjjn = -α * df.pjn .* (1 .- df.sjn_hat)
    return df
end


"""Estimate cross-price elasticities compared to product 3"""
function add_crossprice_elasticities!(reg, df)
    # Get price coefficient
    αidx = findfirst(==("pjn"), reg.coefnames)
    α = -reg.coef[αidx]
    df.ηjjn = -α * df.pjn .* (1 .- df.sjn_hat)
    df = @chain df begin
        filter(:prodid => ==(3), _)
        combine(:market, [:pjn,:sjn_hat] => ((p,s) -> α*p.*s) => :ηj3n)
        leftjoin(df[df.prodid .!=3,[:market,:prodid]], _, on = :market)
        leftjoin(df, _, on = [:market, :prodid])
        sort([:market, :prodid])
    end
    return df
end


#! 2SLS Weight matrix for logit and nested logit (BLP): 
# Nevo Footnote 22: That is, Φ=Z'Z, which is the “optimal” weight matrix under the assumption of homoskedastic errors.







"""Calculate average of other products' characteristics; add columns."""
function add_avg_characteristics!(df)
	df.x1_avg = @chain df begin
		groupby(:market)
		# Avg of other products for j = ((sum of all) - (product j)) / 4
		# collect(:x1 => (x -> (sum(x) .- x) ./ 4) => :x1_avg)
		collect(:x1 => sum => :x1_avg)
        _.x1_avg
	end
end
add_avg_characteristics!(copy(df))






#?==============================================================================
#?==============================================================================
#?==============================================================================
#?                        QUESTION 2: NESTED LOGIT
#?==============================================================================
#?==============================================================================
#?==============================================================================

############# PART A
# Had to run this first: `sudo apt install texlive-luatex`
TikzPictures.standaloneWorkaround(true)
tp = TikzPicture(L"""
	\draw (0,0) -- (10,10);
	\draw (10,0) -- (0,10);
	\node at (5,5) {tikz $\sqrt{\pi}$};"""
	, options="scale=0.25", preamble="")


















#?==============================================================================
#?==============================================================================
#?==============================================================================
#?                        QUESTION 3: MARGINAL COSTS
#?==============================================================================
#?==============================================================================
#?==============================================================================





















#?==============================================================================
#?==============================================================================
#?==============================================================================
#?                        QUESTION 4: PRICING EQUATIONS
#?==============================================================================
#?==============================================================================
#?==============================================================================























#?==============================================================================
#?==============================================================================
#?==============================================================================
#?                        QUESTION 5: RANDOM COEFFICIENTS MODEL
#?==============================================================================
#?==============================================================================
#?==============================================================================



























#*==============================================================================
#*==============================================================================
#*==============================================================================
#*                                   NOTES
#*==============================================================================
#*==============================================================================
#*==============================================================================










