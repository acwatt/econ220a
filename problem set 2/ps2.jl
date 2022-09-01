#==============================================================================
file: ps2.jl
description: estimating model for Problem Set #2 of Econ 220A (UC Berkeley, 2022),
             Demand Estimation and Merger Simulation. The problem set file is
             . The two original files are:
             - PSET1-2019-2.pdf
             - ps1data.csv
author: Aaron C Watt (UCB Grad Student, Ag & Resource Econ)
notes:
    Moved to pluto notebook ( ps2-notebook.jl):
    - problem 1
    Working on in this file:
    - problem 2
Citations:
    - [0] Julia: A Fresh Approach to Numerical Computing. Jeff Bezanson, Alan Edelman, Stefan Karpinski, Viral B. Shah. (2017) SIAM Review, 59: 65–98. doi: 10.1137/141000671.
    - [1] Steven T. Berry, “Estimating Discrete-Choice Models of Product Differentiation,” The RAND Journal of Economics, 1994, 242–62.
    - [2] Aviv Nevo, “A Practitioner’s Guide to Estimation of Random-Coefficients Logit Models of Demand,” Journal of Economics & Management Strategy 9, no. 4 (Winter 2000): 513–48, https://doi.org/10.1162/105864000567954.
    - [3] Daniel A Ackerberg and Gregory S Crawford, “Estimating Price Elasticities in Diﬀerentiated Product Demand Models with Endogenous Characteristics,” Working Paper, 2009.
    - [4] Ackerberg and Crawford, “Estimating Price Elasticities in Diﬀerentiated Product Demand Models with Endogenous Characteristics.” (2009) Working Paper
    - [5] Kenneth E Train, Discrete Choice Methods with Simulation, 2nd ed., 2009.
==============================================================================#

#?==============================================================================
#?                                PACKAGES
#?==============================================================================
# For Stats stuff
using Statistics, StatsBase, StatsPlots, StatsFuns, GLM
using CovarianceMatrices, StatsModels, FixedEffectModels
using FiniteDifferences, NLsolve
# For dataframes
using DataFrames, CSV, DataFramesMeta
using ShiftedArrays
using OrderedCollections
# Random Variable distributions
using Random, Distributions
Random.seed!(0);
# Output to files
using Latexify
using RegressionTables
using Luxor: Point,arrow,@draw,ellipse,circle,fontsize,sethue,label
# Optimizing functions
# using Optim, NLSolversBase
# using LinearAlgebra: diag
# using JuMP, KNITRO




#?=============================================================================
#?                                 SETUP
#?=============================================================================
# Load data
cd(dirname(@__FILE__))
df = DataFrame(CSV.File(string(dirname(@__FILE__), "/ps1data.csv")))








#?=============================================================================
#?=============================================================================
#?=============================================================================
#?                        QUESTION 1: LOGIT
#?=============================================================================
#?=============================================================================
#?=============================================================================

#*#############################################################################
#*####################### QUESTION 1 FUNCTIONS ################################
assign(df, args...) = DataFramesMeta.transform(df, args...)


"""Create the logit dependent variable: δⱼ,ₙ = log(sjn) - log(sj0);  sj0 = outside share (1-Σsjn)"""
function add_logit_depvar!(df)
    @chain df begin
        groupby(:market)
        @transform!(:sj0 = 1 .- sum(:sjn))
        @transform!(:δjn = log.(:sjn) .- log.(:sj0))
    end
end


eq6(x) = exp.(x) ./ (1 + sum(exp.(x)))

"""Estimates shares using linear paremeters (Nevo Eq 6) using regress results."""
function add_predicted_shares!(reg, df_)
    df = copy(df_)
    # Predict δ for all observations (xβ-αp)
    df.δ = predict(reg, df)
    # Estimate share using Eq (6)
    df_.sjn_hat = @chain df begin
        groupby(:market)
        combine(:δ => eq6 => :sjn_hat)
        _[!,:sjn_hat]
    end
    return df_
end


"""Estimate own-price and cross-price elasticities of the predicted market shares for Conditional Logit; add columns.

    Only adds cross-price elasticities relative to product 3.
    ηjjn and ηj3n
    Nevo (2000) pg 522
    Need to use the predicted market share (̂sⱼₙ = Xβ)
"""
function add_price_elasticities!(reg, df)
	sort!(df, [:market, :prodid])
    # Get price coefficient
    α = -get_coef(reg, "pjn")
    df.ηjjn = -α * df.pjn .* (1 .- df.sjn_hat)
    df.ηj3n = @chain df begin
        filter(:prodid => ==(3), _)
        combine(:market, [:pjn,:sjn_hat] => ((p,s) -> α*p.*s) => :ηj3n)
        leftjoin(df[df.prodid .!=3,[:market,:prodid]], _, on = :market)
        leftjoin(df, _, on = [:market, :prodid])
        sort([:market, :prodid])
		_[!,:ηj3n]
	end
    add_marginal_costs_NL!(df)
    return df
end


"""Calculate average of other products' characteristics; add columns."""
function add_avg_characteristics!(df)
	@chain df begin
		groupby(:market)
		# Avg of other products for j = ((sum of all) - (product j)) / 4
		combine(:x1 => (x -> (sum(x) .- x) ./ 4) => :x1a,
				:x2 => (x -> (sum(x) .- x) ./ 4) => :x2a)
		insertcols!(df, :x1_avg => _.x1a, :x2_avg => _.x2a) 
	end
end


"""Return coefficient from regression for variable var"""
function get_coef(reg, var)
    idx = findfirst(==(var), reg.coefnames)
    return reg.coef[idx]
end


"""Create the within-group share variable ̄sⱼₕ"""
function add_ingroup_share!(df)
    # Add group identifier (problem 2)
    group_dict = Dict(1=>1, 2=>2, 3=>2, 4=>3, 5=>3)
    # Calculate in-group share within each market-group
    @chain df begin
        @transform!(:group = get.(Ref(group_dict), :prodid, missing))
        groupby([:market, :group])
        @transform!(:sjhn = :sjn ./ sum(:sjn))
    end
end


"""Add estimated marginal costs to dataframe.
    Requires the own-price elasticity :ηjjn column.
"""
add_marginal_costs_NL!(df) = @transform!(df, :mc = :pjn .* (1 .+ 1 ./ :ηjjn))


# Add δ (mean utility levels)
df = add_logit_depvar!(df)
# Add within group-market shares ̄sⱼₕ
df = add_ingroup_share!(df)



#*#############################################################################
#*####################### QUESTION 1 WORK #####################################

#!############ PART A
#! Estimate an aggregate Logit model using OLS (Product 3 is reference group)
df1a = copy(df)
# Regress δ on X to estimate α,β that minimizes error term
# δ = mean "observed" utility, Xβ-αp = mean predicted utility (α,β are mean utility parameters)
reg1a = reg(df1a, @formula(δjn ~ pjn + d1 + d2  + d4 + d5 + x1 + x2), Vcov.cluster(:market))
# Predict the shares using Xβ-αp (Nevo Eq. 6)
add_predicted_shares!(reg1a, df1a)
# Calculate estimated own-price elasticities for each product 
# Calculate estimated cross-price elasticity with respect to product 3.
add_price_elasticities!(reg1a, df1a)
df1a[!, [:market, :sjn, :δjn, :ηjjn, :ηj3n, :prodid]]



#!############ PART B
#! Estimate the same Logit model using IV with cost-shifters as instruments
df1b = copy(df)
# IV Regress δ on X to estimate α,β that minimizes error term, Cost Shifter istruments (w's)
# δ = mean "observed" utility, Xβ-αp = mean predicted utility (α,β are mean utility parameters)
reg1b_ss = reg(df1b, @formula(δjn ~ (pjn ~ w1 + w2) + d1 + d2  + d4 + d5 + x1 + x2), Vcov.cluster(:market))
# Predict the shares using Xβ-αp (Nevo Eq. 6)
add_predicted_shares!(reg1b_ss, df1b)
# Calculate estimated own-price elasticities for each product 
# Calculate estimated cross-price elasticity with respect to product 3.
add_price_elasticities!(reg1b_ss, df1b)


#!############ PART C
#! Estimate the same Logit model using IV with the average of the characteristics of other products as instruments
df1c = copy(df)
# Create instrument: avg characteristics of other products
add_avg_characteristics!(df1c)
# Regress δ on X to estimate α,β that minimizes error term, avg other product x's (x_avg's)
# δ = mean "observed" utility, Xβ-αp = mean predicted utility (α,β are mean utility parameters)
reg1c = reg(df1c, @formula(δjn ~ (pjn ~ x1_avg + x2_avg) + d1 + d2  + d4 + d5 + x1 + x2), Vcov.cluster(:market))
# Predict the shares using Xβ-αp (Nevo Eq. 6)
add_predicted_shares!(reg1c, df1c)
# Calculate estimated own-price elasticities for each product 
# Calculate estimated cross-price elasticity with respect to product 3.
add_price_elasticities!(reg1c, df1c)



# Average (across markets) elasticities of OLS model
@chain df1a begin
	groupby(:prodid)
	combine(:ηjjn => mean => :ηjj, :ηj3n => mean => :ηj3)
	latexify(_, fmt=FancyNumberFormatter(3))
end

# Average (across markets) elasticities of IV model
@chain df1b begin
	groupby(:prodid)
	combine(:ηjjn => mean => :ηjj, :ηj3n => mean => :ηj3)
	# rename(:ηjj => :ηⱼⱼ, :ηj3 => :ηⱼ₃)
	latexify(_, fmt=FancyNumberFormatter(3))
end



####################################################################
# Below are latex table outputs. These are not currenlty used,
# because the default output in the Pluto notebook is good enough.
# May come back to use these if not submitting the pluto notebook.
####################################################################
# First stage cost-shift instruments of IV model
reg1b_fs = reg(df, @formula(pjn ~ w1 + w2 + d1 + d2  + d4 + d5 + x1 + x2), Vcov.cluster(:market))

# Fist stage average characteristic instruments of IV model
reg1c_fs = reg(df1c, @formula(pjn ~ x1_avg + x2_avg + d1 + d2  + d4 + d5 + x1 + x2), Vcov.cluster(:market))

#=
repl_dict = Dict("δjn" => "\$\\delta_{jn}\$", "_avg" => "\$_{\\text{avg}}\$")
regtable(reg1a,reg1b_ss,reg1c; renderSettings = asciiOutput(), regression_statistics = [:nobs, :r2, :f_kp], transform_labels = repl_dict)
regtable(reg1a,reg1b_ss,reg1c; renderSettings = latexOutput("1-secondstages.tex"), regression_statistics = [:nobs, :r2, :f_kp], transform_labels = repl_dict)
regtable(reg1b_fs,reg1c_fs; renderSettings = asciiOutput(), regression_statistics = [:nobs, :r2, :f, :f_kp], transform_labels = repl_dict)
regtable(reg1b_fs,reg1c_fs; renderSettings = latexOutput("1-firststages.tex"), regression_statistics = [:nobs, :r2, :f, :f_kp], transform_labels = repl_dict)

# Testing putting first and second stages in same table
regtable(reg1a,reg1b_ss,reg1c,reg1b_fs,reg1c_fs;
         renderSettings = asciiOutput(),
         regression_statistics = [:nobs, :r2, :f], 
         transform_labels = repl_dict,
         groups=["2nd Stage" "2nd Stage" "2nd Stage" "1st Stage" "1st Stage"])
regtable(reg1a,reg1b_ss,reg1c,reg1b_fs,reg1c_fs;
         renderSettings = latexOutput("1-bothstages.tex"),
         regression_statistics = [:nobs, :r2, :f], 
         transform_labels = repl_dict,
         groups=["2nd Stage" "2nd Stage" "2nd Stage" "1st Stage" "1st Stage"])
=#












#?==============================================================================
#?==============================================================================
#?==============================================================================
#?                        QUESTION 2: NESTED LOGIT
#?==============================================================================
#?==============================================================================
#?==============================================================================


#!############ PART A
function draw_arrow_diagram()
    s = 50  # overall scale factor (including font size)
    y = 2.5  # vertical scale factor
    x = 1.5  # horizontal scale factor
    p0 = Point(0s*x, -1s*y)
    pairport = Point(-2s*x, 0s*y)
    pnone = Point(0s*x, 0s*y)
    pcity = Point(2s*x, 0s*y)
    p2 = Point(-3s*x, s*y)
    p3 = Point(-s*x, s*y)
    p4 = Point(s*x, s*y)
    p5 = Point(3s*x, s*y)
    lw = 3; ahl = 20;
    myarrow(x,y;ahl=ahl) = arrow(x, y, linewidth=lw, arrowheadlength = ahl)
    @draw begin
        fontsize(0.5s)
        # Arrows
        myarrow(p0, pnone)  # to firm one
        myarrow(p0, pairport; ahl=0)
        myarrow(p0, pcity)
        myarrow(pairport, p2)
        myarrow(pairport, p3)
        myarrow(pcity, p4)
        myarrow(pcity, p5)
        
        # Market circle
        sethue("white")
        ellipse(p0, 2.5s, 1.2s, :fill)
        sethue("black")
        ellipse(p0, 2.5s, 1.2s, :stroke)
        label("Market n", :S, p0, offset=-10)
        
        # Airport circle
        sethue("white")
        circle(pairport, 0.9s, :fill)
        sethue("black")
        circle(pairport, 0.9s, :stroke)
        label("Airport", :S, pairport, offset=-10)
        
        # City circle
        sethue("white")
        circle(pcity, 0.9s, :fill)
        sethue("black")
        circle(pcity, 0.9s, :stroke)
        label("City", :S, pcity, offset=-10)

        # Firm Text
        label("Firm 1", :S, pnone)
        label("Firm 2", :S, p2)
        label("Firm 3", :S, p3)
        label("Firm 4", :S, p4)
        label("Firm 5", :S, p5)
    end 10s*x 3s*y
end



#*#############################################################################
#*####################### QUESTION 2 FUNCTIONS ################################

"""Create the new observed δ analog from BLP Eq. 27: δⱼ(s,σ) = ln(sⱼ) - σ ln(̄sⱼₕ) - ln(s₀)"""
function add_mean_utility_NL!(df::DataFrame, σ::Real)
    # Calculate analytic mean utility from observed shares and estimated σ
    @chain df begin
        groupby(:market)
        @transform!(:δ_obs = log.(:sjn) .- σ*log.(:sjhn) - log.(:sj0))
    end
end


"""Create the new predicted δ from below BLP Eq. 27: δ = xⱼβ -αpⱼ"""
function add_mean_utility_NL!(df::DataFrame, reg::FixedEffectModel)
    # Calculate predicted mean utility from observed characteristics and estimated β,α
    df.δ_hat = @chain df begin
        @transform(:sjhn = 1) # so log(sjhn) = 0, removed σ⋅log(sⱼₕₙ) from prediction
        predict(reg, _)
    end
end


"""Estimate Nested Logit in-group share denominator Dₕ (BLP Eq 23).
    Dₕ = ∑(j∈Jₕ) exp(δⱼ / (1-σ))
    requires δ_hat = X'β  in df
"""
function add_group_denominator_NL!(df, σ)
    @chain df begin
        groupby([:market, :group])
        @transform!(:Dh = sum(exp.(:δ_hat / (1-σ))))
    end
end


"""Estimate Nested Logit group share denominator ∑ₕ Dₕ (BLP Eq 24).
    D_sum = ∑ₕ Dₕ in each market
    requires δ_hat = X'β  in df
    Note: h=0 is outside group, and D₀ = 1 (since δ₀ is normalized to 0, and e⁰=1)
"""
function add_market_demoninator_NL!(df, σ)
    @chain df begin
        groupby(:market)
        @combine(:D_sum = sum(:Dh.^(1-σ)) + 1)
        @select(:market, :D_sum)
        leftjoin!(df, _, on=:market)
    end
end


"""Estimate Nested Logit predicted in-group share (BLP Eq 23).
    sⱼ|ₕ = exp(δⱼ / (1-σ)) / Dₕ
    requires δ_hat = X'β  in df
"""
function add_ingroup_share_NL!(df, σ)
    @chain df begin
        groupby(:market)
        @transform!(:sjhn_hat = exp.(:δ_hat / (1-σ)) ./ :Dh)
    end
end


"""Estimate Nested Logit predicted market share (BLP Eq 25).
    sⱼ = exp(δⱼ / (1-σ)) / (Dₕ^σ * D_sum)
    requires δ_hat = X'β  in df
"""
function add_market_share_NL!(df, σ)
    @chain df begin
        groupby(:market)
        @transform!(:sjn_hat = exp.(:δ_hat / (1-σ)) ./ (:Dh.^σ .* :D_sum))
    end
end


"""Expand prices and market shares for product 3 to all rows to help with elasticity computation"""
function add_product3_price_and_shares!(df)
    @chain df begin
        @subset(:prodid .== 3)
        rename(:sjn_hat => :s3n_hat, :pjn => :p3n)  # Creates a new dataframe
        @select(:market, :s3n_hat, :p3n)
        leftjoin!(df, _, on=:market)        # merge just :s3n and :p3n onto orginal dataframe
    end
end


"""Return Nested Logit cross-price elasticity w.r.t. product 3; checks if in same nest as product 3.
    From [4], pg 15: (supressing the market subscript) 
    ∂sⱼ/∂pₖ = -α sⱼ (1/(1-σ) - σ/(1-σ)*sⱼ|ₕ - sⱼ)  if j = k (same product)
        =   α sₖ (σ/(1-σ)*sⱼ|ₕ + sⱼ)            if j & k in same nest
        =   α sₖ sⱼ                             otherwise
    ηⱼₖ = ∂sⱼ/∂pₖ * pₖ / sⱼ
        = -α pₖ (1/(1-σ) - σ/(1-σ)*sⱼ|ₕ - sⱼ)      if j = k (same product)
        =   α sₖ (σ/(1-σ) * sⱼ|ₕ * pₖ / sⱼ + pₖ)    if j & k in same nest
        =   α sₖ pₖ                                otherwise

    Need to use the predicted market share ̂sⱼₙ(δ_hat) BLP Eq 25
"""
function cross_price_elasticity_NL(α, σ, p3, s3, sjh, sj, productID, nestID)
    if productID == 3  # do not want to return own-price elasticity
        return missing
    elseif nestID == 2  # product 3 is in next 2: return in-nest cross price elasticity
        return α*s3*p3 * (σ/(1-σ) * sjh/sj + 1)
    else
        return α*s3*p3
    end
end
"""Nested Logit cross-price-derivative of market share (within nest)"""
∂sⱼ∂pₖ(α, σ, sⱼ, sⱼ₍ₕ₎, sₖ) = α*sₖ*(σ/(1-σ) * sⱼ₍ₕ₎ + sⱼ)


"""Return Nested Logit own-price elasticity"""
own_price_elasticity_NL(α, σ, pj, sjh, sj) = -α*pj * (1/(1-σ)  - σ/(1-σ)*sjh - sj)
"""Nested Logit own-price-derivative of market share"""
∂sⱼ∂pⱼ(α, σ, sⱼ, sⱼ₍ₕ₎) = -α*sⱼ*(1/(1-σ) - σ/(1-σ)*sⱼ₍ₕ₎ - sⱼ)


"""Estimate own-price and cross-price elasticities of the predicted market shares for Nested Logit; add columns.
    Only adds cross-price elasticities relative to product 3.
    ηjjn and ηj3n
    Requires that p3 (price of product 3) and s3 (market share of product 3) be added as new columns.
"""
function add_price_elasticities_NL!(reg, df)
	sort!(df, [:market, :prodid])
    # Get coefficients
    α = -get_coef(reg, "pjn")
    σ = get_coef(reg, "log(sjhn)")
    @chain df begin
        @transform!(:ηj3n = cross_price_elasticity_NL.(α, σ, :p3n, :s3n_hat, :sjhn_hat, :sjn_hat, :prodid, :group))
        @transform!(:ηjjn = own_price_elasticity_NL.(α, σ, :pjn, :sjhn_hat, :sjn_hat))
    end
end


function report_avg_price_elasticities(df)
    @chain df begin
        groupby(:prodid)
        @combine(:ηjj = mean(:ηjjn), :ηj3 = mean(:ηj3n))
        latexify(fmt=FancyNumberFormatter(3))
    end
end


"""Return Nested Logit price elasticities using dataframe and BLP Eq 28 regression results."""
function estimate_price_elasticities_NL!(df, reg)
    # Create the new observed δ analog from BLP Eq. 27: δⱼ(s,σ) = ln(sⱼ) - σ ln(̄sⱼₕ) - ln(s₀)
    σ = get_coef(reg, "log(sjhn)")
    add_mean_utility_NL!(df, σ)  # adds δ_obs
    # Create predicted δ: δⱼ = xⱼβ - αpⱼ (removing the σ ln(̄sⱼₕ) term)
    add_mean_utility_NL!(df, reg) # adds δ_hat

    # Predict the market and within group shares using predicted δ, σ (using α,β,σ in Nevo Eq's 23, 24, 25)
    # We need predicted market shares to calculate elasticities
    add_group_denominator_NL!(df, σ)
    add_market_demoninator_NL!(df, σ)
    add_ingroup_share_NL!(df, σ)
    add_market_share_NL!(df, σ)

    # Calculate estimated own-price elasticities for each product 
    # Calculate estimated cross-price elasticity with respect to product 3.
    add_product3_price_and_shares!(df)
    add_price_elasticities_NL!(reg, df)
    # Calculate marginal costs
    add_marginal_costs_NL!(df)

    # Report average price elasticities across markets
    report_avg_price_elasticities(df)
end


""""Add 'average' in-nest in-market other-product characteristics.
    Because there is at most 2 products in each market-nest, this 'average'
    is just the product characteristic value of the other product in the nest.
"""
function add_avg_characteristics_NL!(df)
    @chain df begin
        # groupby market and nest: results in list of dataframes with one or two rows
        groupby([:market, :group])
        # In dataframe with two rows, the average characteristic of the other row (the other product)
        # is retrieved by circularly shifting the characteristics down one row
        # (so the top row moves down one and the bottom row moves to the top)
        @transform!(:x1_avg = ShiftedArrays.circshift(:x1, 1))
        @transform!(:x2_avg = ShiftedArrays.circshift(:x2, 1))
        # If product is in nest 1, replace with 0
        @transform!(:x1_avg = ifelse.(:group .== 1, 0, :x1_avg))
        @transform!(:x2_avg = ifelse.(:group .== 1, 0, :x2_avg))
    end
end







#!############ PART B
#= Notes on Nested Logit
    Berry 1994: (referring to groups as h, because there is no unicode subscript g)
    1. group the products into H + 1 exhaustive and mutually exclusive sets, h = 0, 1, ... , H (H=2 here)
    2. Denote set of products in group h as Iₕ: I₀={0 outside product}, I₁={1}, I₂={2,3}, I₃={4,5}
    3. For product j∈Iₕ, utilty of consumer i = uᵢⱼ = δⱼ + ζᵢₕ + (1-σ)ϵᵢⱼ; δⱼ = xⱼβ - αpⱼ + ξⱼ; ϵᵢⱼ iid T1EV
        ζᵢₕ is common to all products in group h for consumer i, ~ dist. depends on σ ∈ [0,1)
        can interpret as Random Coefficients model with ζᵢₕ random coefficients on only group-specific dummy variables
        dⱼₕ = 1 if j ∈ Iₕ
        Then uᵢⱼ = δⱼ + Σₕ [dⱼg ζᵢₕ] + (1-σ)ϵᵢⱼ
    ... skipped some steps here
    4. ln(sⱼ) - ln(s₀) = xⱼβ - αpⱼ + σln(̄sⱼₕ) + ξⱼ; ̄sⱼₕ = observed share of product j in group h
        ̄sⱼₕ = sⱼ / Σⱼ∈ₕ sⱼ
=#

# The group, and in-group share variables are created in question 1 functions
df2 = copy(df)
# We keep the same dependent variable ln(sⱼ)-ln(s₀), but we need to call it something else
# In Nested Logit, ln(sⱼ)-ln(s₀) = δ + σ ln(̄sⱼₕ)
# Change the name of the dependent variable from δ to ln sⱼs₀ = ln(sⱼ)-ln(s₀)
rename!(df2, :δjn => :lnsjs0)
df2b = copy(df2)
# Regress ln(sⱼ)-ln(s₀) on X and ln(̄sⱼₕ) to estimate α,β,σ that minimizes error term ( BLP Eq 28)
reg2b = reg(df2b, @formula(lnsjs0 ~ pjn + log(sjhn) + d1 + d2  + d4 + d5 + x1 + x2), Vcov.cluster(:market))
# Calculate own and cross-price elasticities
estimate_price_elasticities_NL!(df2b, reg2b)










#!############ PART C
#! Instrumental Variables (IV) with cost-shifters
#= Notes on Nested Logit
    2SLS Weight matrix for logit and nested logit (BLP): 
        Nevo Footnote 22: That is, Φ=Z'Z, which is the “optimal” weight matrix under the assumption of homoskedastic errors.
        I think this is taken care of when using 2SLS
        Beta note: β = beta in this font (depends on the IDE used to read this file)
=#

df2c = copy(df2)
# Regress ln(sⱼ)-ln(s₀) on X and ln(̄sⱼₕ) to estimate α,β,σ that minimizes error term (BLP Eq 28)
reg2c_fs1 = reg(df2c, @formula(pjn ~ w1 + w2 + d1 + d2  + d4 + d5 + x1 + x2), Vcov.cluster(:market))
reg2c_fs2 = reg(df2c, @formula(log(sjhn) ~ w1 + w2 + d1 + d2  + d4 + d5 + x1 + x2), Vcov.cluster(:market))
reg2c_ss = reg(df2c, @formula(lnsjs0 ~ (pjn + log(sjhn) ~ w1 + w2) + d1 + d2  + d4 + d5 + x1 + x2), Vcov.cluster(:market))
# Calculate own and cross-price elasticities
estimate_price_elasticities_NL!(df2c, reg2c_ss)










#!############ PART D
#! Instrumental Variables (IV) with average of characteristics of other products within the group
#= From the problem set instructions:
You will have to create this set of instruments from the raw
data. Note the average x1 for product 2 is just the x1
for product 3 (and vice versa). There is no average other
characteristic for product 1 as it is in its own group, so just
give it a value of 0.
=#
df2d = copy(df2)

# Add the instruments: average other-nest product characteristics
add_avg_characteristics_NL!(df2d)

# Regress ln(sⱼ)-ln(s₀) on X and ln(̄sⱼₕ) to estimate α,β,σ that minimizes error term (BLP Eq 28)
reg2d_fs1 = reg(df2d, @formula(pjn ~ x1_avg + x2_avg + d1 + d2  + d4 + d5 + x1 + x2), Vcov.cluster(:market))
reg2d_fs2 = reg(df2d, @formula(log(sjhn) ~ x1_avg + x2_avg + d1 + d2  + d4 + d5 + x1 + x2), Vcov.cluster(:market))
reg2d_ss = reg(df2d, @formula(lnsjs0 ~ (pjn + log(sjhn) ~ x1_avg + x2_avg) + d1 + d2  + d4 + d5 + x1 + x2), Vcov.cluster(:market))
# Calculate own and cross-price elasticities
estimate_price_elasticities_NL!(df2d, reg2d_ss)




#!############ PART E
#! Instrumental Variables (IV) with with both cost-shifters and 
#! within-group average of rivals' characteristics as in-struments.
df2e = copy(df2)
add_avg_characteristics_NL!(df2e)

# Regress ln(sⱼ)-ln(s₀) on X and ln(̄sⱼₕ) to estimate α,β,σ that minimizes error term (BLP Eq 28)
reg2e_fs1 = reg(df2e, @formula(pjn ~ w1 + w2 + x1_avg + x2_avg + d1 + d2  + d4 + d5 + x1 + x2), Vcov.cluster(:market))
reg2e_fs2 = reg(df2e, @formula(log(sjhn) ~ w1 + w2 + x1_avg + x2_avg + d1 + d2  + d4 + d5 + x1 + x2), Vcov.cluster(:market))
reg2e_ss = reg(df2e, @formula(lnsjs0 ~ (pjn + log(sjhn) ~ w1 + w2 + x1_avg + x2_avg) + d1 + d2  + d4 + d5 + x1 + x2), Vcov.cluster(:market))
# Calculate own and cross-price elasticities
estimate_price_elasticities_NL!(df2e, reg2e_ss)
























#?==============================================================================
#?==============================================================================
#?==============================================================================
#?                        QUESTION 3: MARGINAL COSTS
#?==============================================================================
#?==============================================================================
#?==============================================================================
#=
Berry 1994: MC ≈ pⱼ * (1 + 1 / ηjj)

- need to estimate marginal costs for each product in each market for all 4 nested logit Models
--> Added to problem 2 functions so MC are estimated right after share elasticities
- need to calculate % of MC that are <0 for each product in each model
- need to combine these in to single dataframe and latexify
=#
"""Return vector of shares of marginal cost estimates that are negative, by product."""
function get_share_neg(df)
    @chain df begin
        groupby(:prodid)
        @combine(:share_neg = sum(:mc .< 0) / length(:mc))
        _.share_neg
    end
end


"""Return vector of shares of marginal cost estimates that are above price, by product."""
function get_share_above_price(df)
    @chain df begin
        groupby(:prodid)
        @combine(:share_above_p = sum(:mc .> :pjn) / length(:mc))
        _.share_above_p
    end
end


mc_df = DataFrame(prodid = 1:5)
df_dict = OrderedDict("b" => df2b, "c" => df2c, "d" => df2d, "e" => df2e)

for (part, d) ∈ df_dict
    @transform!(mc_df, $("MC<0_"*part) = get_share_neg(d), $("MC>p_"*part) = get_share_above_price(d))
end
latexify(mc_df, fmt=FancyNumberFormatter(3))

@chain df2b begin
    @transform(:posMC = :mc .> 0)
    groupby(:posMC)
    @combine(:avg_p = mean(:mc))
    latexify(fmt=FancyNumberFormatter(3))
end

#!#####################################################################################################################
#!              NEXT
#!#####################################################################################################################
#=
- Check if I included group 0 (outside) in the estimation of s_g (or maybe just the denomintator)
- editing add_market_demoninator_NL! to add outside group in denominator
=#





















#?==============================================================================
#?==============================================================================
#?==============================================================================
#?                        QUESTION 4: PRICING EQUATIONS
#?==============================================================================
#?==============================================================================
#?==============================================================================
#=
- need to solve for p in system of equations (can group then subset to just goods 2 and 3, then add new prices as vector solution to matrix equation)
- Given α,σ,β,X,p₁,p₄,p₅, can write market and within nest shares as functions of p₂ and p₃ (this includes replacing p2,p3 in DF and using -predict-)
- Need full vector of X's used in IV reg
=#



#! Goal: predict new prices after firms 2 and 3 merge, for each market
#=BLP 1995 Eq 3.3 has the FOC for multi-product firms
  Using the estimated Nested Logit parameters, predict new p2,p3 price equilibrium
  assuming p1,p4,p5, marginal costs, product characteristics, and estimated parameters are fixed
  Need to solve a system of 2 equations (2 FOC's for the merged firm, w.r.t. p2 and p3)
  These are non-linear equations -- taking p2 and p3 derivatives of s2 and s3 (market shares)

  Have data on p1, p3, p4, X, W, D for all products
  Have reg estimates of α, β, σ
  Predict δⱼ:     Function of p2, p3, reg: replace p2, p3 in DF and predict δ
  Predict Dh:     Function of δⱼ, σ
  Predict Dh sum: Function of Dh, σ
  Predict sj|h:   Function of δⱼ, Dh, σ
  Predict Sj:     Function of δ, Dh, Dh sum, σ
  Estimate ∂sⱼ/∂pₖ for use in BLP equation 3.3
  Predict new prices and % change from old prices

  I have the regression coefficients from IV part E
  Write a function of p2,p3 that equals 0, solved by NLsolve, using observed p2,p3 as starting points
  Try using nlsolve(f!, initial_x, autodiff = :forward)
  where f!(F, x*) = 0 and f!(F,x) = F[1],F[2] = somefunctionof(x[1], x[2])
   --> autodiff did not work because using `predict()`. Would probably work
       if I extracted the matrix and used matrix multiplication
       Sticking with finite difference numerical differentiation instead
       Should be able to apply to any sⱼ function of p
=#

#! LOGIT -- use reg1b_ss for δ predictions

"""Return new p2 or p3 if product id is 2 or 3, respectively. Otherwise, return old price p."""
p2p3(prodid, p2, p3, p) = ifelse(prodid == 2, p2, ifelse(prodid == 3, p3, p))
"""Replace observed prices for products 2 and 3 in single-market dataframe (5 rows) with input p2 and p3"""
replace_p2p3(df, p2, p3) = @transform(df, :pjn = p2p3.(:prodid, p2, p3, :pjn))

"""Estimate Logit mean utility as a function of new p2, p3, given data and IV estimates"""
δ(df, p2, p3) = @chain df replace_p2p3(p2, p3) predict(reg1b_ss, _)

"""Estimate market shares s2,s3 based on p2, p3 -> δ(p2,p3)
    type = "logit", "nested logit"
    Need to set df.type
"""
function s(df, p2, p3)
    type = first(df.type)
    if type == "logit"
        eq6(δ(df, p2, p3))[2:3]
    elseif type == "nested logit"
        merger_NL(df, p2, p3) 
    end
end
s(df, p) = s(df, p[1], p[2])
s2(p, df) = s(df, p)[1];  s3(p, df) = s(df, p)[2]

"""Estimate price derivatives of market shares at p2,p3"""
mygrad(f,x,df) = FiniteDifferences.grad(central_fdm(5, 1), x -> f(x, df), x)
∂s2∂p2(p, df) = mygrad(s2, p, df)[1][1];  ∂s2∂p3(p, df) = mygrad(s2, p, df)[1][2]
∂s3∂p2(p, df) = mygrad(s3, p, df)[1][1];  ∂s3∂p3(p, df) = mygrad(s3, p, df)[1][2]

"""Residuals function to find root -- find p2,p3 that sets F[1] = 0 = F[2]"""
function price_eqations!(F, p, df)
    mc = df.mc[2:3]
    F[1] = s2(p, df) + (p[1] - mc[1])*∂s2∂p2(p, df) + (p[2] - mc[2])*∂s3∂p2(p, df)
    F[2] = s3(p, df) + (p[1] - mc[1])*∂s2∂p3(p, df) + (p[2] - mc[2])*∂s3∂p3(p, df)
end

"""Return vector of 5 product prices after merger of firms 2 and 3"""
function p_soln(df)
    sort!(df, :prodid)
    p0 = df.pjn[2:3]
    # Solve system of equations for (p2,p3), using finite differencing
    f!(F, p) = price_eqations!(F, p, df)
    p2p3 = nlsolve(f!, p0).zero
    return [df.pjn[1]; p2p3; df.pjn[4:5]]
end

"""Return percentage point increase from old to new"""
percent_inc(p_old, p_new) = (p_new / p_old - 1) * 100

df1z = copy(df1b)
df1z[!, :type] .= "logit"
# Add new prices for products 2 and 3 after merger of firms 2 and 3
df1z.pjn_merge = combine(p_soln, groupby(df1z, :market)).x1
# Add percentage point increases in prices
@chain df1z @transform!(:pjn_incr = percent_inc.(:pjn, :pjn_merge))
# Calculate average percentage point increase for products 2 and 3
@chain df1z begin
    @subset(:prodid .∈ [2:3])
    groupby(:prodid)
    @combine($("Mean Price Increase (ppt)") = mean(:pjn_incr))
    latexify(fmt=FancyNumberFormatter(3))
end



#! NESTED LOGIT -- use reg2e_ss for δ predictions

"""Estimate Nested Logit δ as a function of new p2, p3, given data and IV estimates.
    Removing the σ⋅log(sⱼₕₙ) term from the regression prediction to get δ instead of ln(sj/s0)
"""
δNL(df, p2, p3) = @chain df replace_p2p3(p2, p3) @transform(:sjhn = 1) predict(reg2e_ss, _)

"""Estimate Nested Logit market shares of goods 2 and 3 for a single-market dataframe"""
function merger_NL(d, p2, p3)
    σ = get_coef(reg2e_ss, "log(sjhn)")
    df = copy(d)
    df.δ_hat = δNL(d, p2, p3)
    @chain df begin
        groupby(:group)
        @transform!(:Dh = sum(exp.(:δ_hat / (1-σ))))
        @transform!(:D_sum = sum(:Dh.^(1-σ)) + 1)
        @transform!(:sjn_hat = exp.(:δ_hat / (1-σ)) ./ (:Dh.^σ .* :D_sum))
    end
    df.sjn_hat[2:3]
end


df2z = copy(df2e)
df2z[!, :type] .= "nested logit"
# Add new prices for products 2 and 3 after merger of firms 2 and 3
@time df2z.pjn_merge = combine(p_soln, groupby(df2z, :market)).x1
# Add percentage point increases in prices
@chain df2z @transform!(:pjn_incr = percent_inc.(:pjn, :pjn_merge))
# Calculate average percentage point increase for products 2 and 3
@chain df2z begin
    @subset(:prodid .∈ [2:3])
    groupby(:prodid)
    @combine($("Mean Price Increase (ppt)") = mean(:pjn_incr))
    latexify(fmt=FancyNumberFormatter(3))
end




















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


#logit = glm(fm, train, Binomial(), ProbitLink())





# USING R mlogit PACKAGE
# https://cran.r-project.org/web/packages/mlogit/vignettes/e2nlogit.html
R"""
library("mlogit")
data("HC", package = "mlogit")
HC <- dfidx(HC, varying = c(2:8, 10:16), choice = "depvar")
cooling.modes <- idx(HC, 2) %in% c('gcc', 'ecc', 'erc', 'hpc')
room.modes <- idx(HC, 2) %in% c('erc', 'er')
# installation / operating costs for cooling are constants, 
# only relevant for mixed systems
HC$icca[! cooling.modes] <- 0
HC$occa[! cooling.modes] <- 0
# create income variables for two sets cooling and rooms
HC$inc.cooling <- HC$inc.room <- 0
HC$inc.cooling[cooling.modes] <- HC$income[cooling.modes]
HC$inc.room[room.modes] <- HC$income[room.modes]
# create an intercet for cooling modes
HC$int.cooling <- as.numeric(cooling.modes)
# estimate the model with only one nest elasticity
nl <- mlogit(depvar ~ ich + och +icca + occa + inc.room + inc.cooling + int.cooling | 0, HC,
             nests = list(cooling = c('gcc','ecc','erc','hpc'), 
             other = c('gc', 'ec', 'er')), un.nest.el = TRUE)  # shared nest elast.
summary(nl)
"""

R"""
# estimate the model with two nest elasticities
nl2 <- mlogit(depvar ~ ich + och +icca + occa + inc.room + inc.cooling + int.cooling | 0, HC,
             nests = list(cooling = c('gcc','ecc','erc','hpc'), 
             other = c('gc', 'ec', 'er')), un.nest.el = FALSE)  # diff. nest elast.
summary(nl2)
"""

R"""
# estimate the model with different nest elasticities
nl3 <- mlogit(depvar ~ ich + och , HC,
             nests = list(cooling = c('gcc','ecc','erc','hpc'), 
             other = c('gc', 'ec', 'er')), un.nest.el = TRUE)  # different nest elast.
summary(nl3)
"""











#! DATAFRAMES MANIPULATION

# Example 1: create 
d = DataFrame(A=1:8, j=repeat(1:4,2), t=repeat(1:2, inner=4), g=repeat(1:2, inner=2, outer=2))
push!(d,     [9,     5,               1,                      3])
@chain d begin
    groupby([:t, :g])
    @combine :D_mean = mean(:A)
    @aside(println("after combine\n", _))
    groupby(:t)
    @combine :D_sum = sum(:D_mean)
    @select(:t, :D_sum)
    leftjoin!(d, _, on=:t)
end
@chain d groupby([:t, :g]) @transform!(:D_mean = mean(:A))


d = DataFrame(A=1:8, j=repeat(1:4,2), t=repeat(1:2, inner=4), g=repeat(1:2, inner=2, outer=2))
push!(d,     [9,     5,               1,                      3])
using ShiftedArrays
@chain d begin
    # groupby market and nest: results in list of dataframes with one or two rows
    groupby([:t, :g])
    # In dataframe with two rows, the average characteristic of the other row (the other product)
    # is retrieved by circularly shifting the characteristics down one row
    # (so the top row moves down one and the bottom row moves to the top)
    @transform!(:x = ShiftedArrays.circshift(:A, 1))
    @aside(println("After Shift, before replacing group 3\n", _))
    # If product is in nest 1, replace with 0
    @transform!(:x = ifelse.(:g .== 1, 0, :x))
end





#! USING RCALL TO COMPARE REGRESSIONS

Random.seed!(107)
d_ = @chain DataFrame(y=1:10) begin
    @transform(:x = 3*:y .+ rand(10)*5)
    @transform(:z = 2*:x .+ rand(10)*10)
    @transform(:w = 5*:y .+ rand(10)*10)
end
println("\n\nStarting")
println("x ~ z + w")
println(reg(d_, @formula(x ~ z + w)))
println("x ~  z")
println(reg(d_, @formula(x ~  z)))
println("y ~ (x ~ z) + w")
println(reg(d_, @formula(y ~ (x ~ z) + w)))
println("y ~ z + w")
println(reg(d_, @formula(y ~ z + w)))

using RCall
# Had to install lfe library in the native R console. RCall install.packages seems to be broken.
R"""
library(fixest)
m1 = feols(x ~ z + w, data=$d_)
m2 = feols(y ~ w | x ~ z, data=$d_)
print(summary(m1))
print(fitstat(m1, "f"))
print("")
print(summary(m2))
"""

"""To use R package `mlogit`: Needed to run `install.packages("mlogit")` in R terminal. Would not work from here or in Julia terminal."""




#! SOLVING NONLINEAR EQUATIONS
using NLsolve
function f!(F, x)
    F[1] = (x[1]+3)*(x[2]^3-7)+18
    F[2] = sin(x[2]*exp(x[1])-1)
end
nlsolve(f!, [0.1; 1.2], autodiff = :forward)

# Solving nonlinear equation with extra function arguments
function g!(F, x, p)
    F[1] = x[1]^2 * (x[2]-3)*p[1]
    F[2] = x[1]^3 * 8*x[2]*p[2]
end
# anonomous function
nlsolve((F,x) -> g!(F,x,p1), [0.1; 1.2], autodiff = :forward)

# wrapper function
newg!(F, x) = g!(F, x, p)
nlsolve(newg!, [0.1; 1.2], autodiff = :forward)




