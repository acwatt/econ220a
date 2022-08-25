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
# For dataframes
using DataFrames, CSV, DataFramesMeta
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






#*#############################################################################
#*####################### QUESTION 1 WORK #####################################

# Add δ (mean utility levels)
df = add_logit_depvar!(df)
# Add within group-market shares ̄sⱼₕ
df = add_ingroup_share!(df)


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


# Create table of price elasticities

















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
"""
function add_market_demoninator_NL!(df, σ)
    @chain df begin
        groupby(:market)
        @combine(:D_sum = sum(:Dh.^(1-σ)))
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


"""Return Nested Logit own-price elasticity"""
own_price_elasticity_NL(α, σ, pj, sjh, sj) = -α*pj * (1/(1-σ)  - σ/(1-σ)*sjh - sj)


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
function estimate_price_elasticities_NL(df_, reg)
    df = copy(df_)
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

    # Report average price elasticities across markets
    report_avg_price_elasticities(df)
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
df2b = copy(df)
# We keep the same dependent variable ln(sⱼ)-ln(s₀), but we need to call it something else
# In Nested Logit, ln(sⱼ)-ln(s₀) = δ + σ ln(̄sⱼₕ)
# Change the name of the dependent variable from δ to ln sⱼs₀ = ln(sⱼ)-ln(s₀)
rename!(df2b, :δjn => :lnsjs0)
# Regress ln(sⱼ)-ln(s₀) on X and ln(̄sⱼₕ) to estimate α,β,σ that minimizes error term ( BLP Eq 28)
reg2b = reg(df2b, @formula(lnsjs0 ~ pjn + log(sjhn) + d1 + d2  + d4 + d5 + x1 + x2), Vcov.cluster(:market))
# Calculate own and cross-price elasticities
estimate_price_elasticities_NL(df2b, reg2b)










#!############ PART C
#! Instrumental Variables (IV) with cost-shifters
#= Notes on Nested Logit
    2SLS Weight matrix for logit and nested logit (BLP): 
        Nevo Footnote 22: That is, Φ=Z'Z, which is the “optimal” weight matrix under the assumption of homoskedastic errors.
        I think this is taken care of when using 2SLS
        Beta note: β = beta in this font (depends on the IDE used to read this file)
=#

df2c = copy(df2b)
# Regress ln(sⱼ)-ln(s₀) on X and ln(̄sⱼₕ) to estimate α,β,σ that minimizes error term (BLP Eq 28)
reg2c_fs = reg(df2c, @formula(pjn ~ w1 + w2 + log(sjhn) + d1 + d2  + d4 + d5 + x1 + x2), Vcov.cluster(:market))
reg2c_ss = reg(df2c, @formula(lnsjs0 ~ (pjn ~ w1 + w2) + log(sjhn) + d1 + d2  + d4 + d5 + x1 + x2), Vcov.cluster(:market))
# Calculate own and cross-price elasticities
estimate_price_elasticities_NL(df2c, reg2c_ss)
# Report both first stage regressions, for pⱼ and ln(̄sⱼₕ)








#!#####################################################################################################################
#!              NEXT
#!#####################################################################################################################
#=
- correct all statements about first stage F stat to the f stat reported in the IV output
- finish part c
=#


#!############ PART D
#! Instrumental Variables (IV) with average of characteristics of other products within the group
#=
You will have to create this set of instruments from the raw
data. Note the average x1 for product 2 is just the x1
for product 3 (and vice versa). There is no average other
characteristic for product 1 as it is in its own group, so just
give it a value of 0.

gropuby(:market)

if :prodid .== 1
    return 0
elseif :prodid .== 2
    return $var

Not sure how to get e.g. prod 3 x1 for prod 2... maybe need to google "vertical lookup"?
=#






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
@chain d 






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
















