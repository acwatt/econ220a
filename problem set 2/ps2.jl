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
    df.δ = predict(reg, df)
    # Estimate share using Eq (6)
    df_.sjn_hat = @chain df begin
        groupby(:market)
        combine(:δ => eq6 => :sjn_hat)
        _[!,:sjn_hat]
    end
    return df_
end


"""Estimate own-price and cross-price elasticities of the predicted market shares; add columns.

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


assign(df, args...) = DataFramesMeta.transform(df, args...)



"""Create the within-group share variable ̄sⱼₕ"""
function add_ingroup_share!(df)
    # Calculate in-group share within each market-group
    df.sjhn = @chain df begin
        groupby([:market, :group])
        combine(:sjn => (x -> (x ./ sum(x))) => :sjhn)
        _.sjhn
    end
    return df
end






#*#############################################################################
#*####################### QUESTION 1 WORK #####################################

# Add δ (mean utility levels)
df = add_logit_depvar!(df)
# Add group identifier (problem 2)
group_dict = Dict(1=>1, 2=>2, 3=>2, 4=>3, 5=>3)
df = @chain df assign(:prodid => (x -> get.(Ref(group_dict), x, missing)) => :group)
# Add within group-market shares ̄sⱼₕ
df = add_ingroup_share!(df)


############# PART A
# Estimate an aggregate Logit model using OLS (Product 3 is reference group)
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



############# PART B
# Estimate the same Logit model using IV with cost-shifters as instruments
df1b = copy(df)
# IV Regress δ on X to estimate α,β that minimizes error term, Cost Shifter istruments (w's)
# δ = mean "observed" utility, Xβ-αp = mean predicted utility (α,β are mean utility parameters)
reg1b_ss = reg(df1b, @formula(δjn ~ (pjn ~ w1 + w2) + d1 + d2  + d4 + d5 + x1 + x2), Vcov.cluster(:market))
# Predict the shares using Xβ-αp (Nevo Eq. 6)
add_predicted_shares!(reg1b_ss, df1b)
# Calculate estimated own-price elasticities for each product 
# Calculate estimated cross-price elasticity with respect to product 3.
add_price_elasticities!(reg1b_ss, df1b)


############# PART C
# Estimate the same Logit model using IV with the average of the characteristics of other products as instruments
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


############# PART A
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

############# PART B
#! 2SLS Weight matrix for logit and nested logit (BLP): 
# Nevo Footnote 22: That is, Φ=Z'Z, which is the “optimal” weight matrix under the assumption of homoskedastic errors.

#= Notes on Nested Logit
2SLS Weight matrix for logit and nested logit (BLP): 
Nevo Footnote 22: That is, Φ=Z'Z, which is the “optimal” weight matrix under the assumption of homoskedastic errors.

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


df2b = copy(df)
# Regress δ on X and ln(̄sⱼₕ) to estimate α,β,σ that minimizes error term
reg2b = reg(df2b, @formula(δjn ~ pjn + log(sjhn) + d1 + d2  + d4 + d5 + x1 + x2), Vcov.cluster(:market))
# Predict the market and within group shares using predicted δ, σ (using α,β,σ in Nevo Eq's 23, 24, 25)
add_predicted_NL_shares!(reg2b, df2b)
# Calculate estimated own-price elasticities for each product 
# Calculate estimated cross-price elasticity with respect to product 3.
add_price_elasticities!(reg1a, df1a)
df1a[!, [:market, :sjn, :δjn, :ηjjn, :ηj3n, :prodid]]








############# PART C
# 


# Report both first stage regressions, for pⱼ and ln(̄sⱼₕ)











using StatsModels: predict
using DataFrames: groupby

"""Estimate Nested Logit market shares and within-group shares using Nevo Eq's 23,24,25"""
function add_predicted_NL_shares!(reg, df)
    df1 = copy(df)
    σ = get_coef(reg2b, "log(sjhn)")
    # Predict mean utility δ (xβ-αp)
    df1.δ = @chain df begin
        assign(:sjhn => (x->1) => :sjhn) # so log(sjhn) = 0, removed σ⋅log(sⱼₕₙ) from prediction
        predict(reg2b,_)
    end
    # Add denominator Dₕ from Eq 23
    # df.sjn_hat, df.sjhn_hat = 
    @chain df1 begin
        # Add denominator
        groupby([:market, :group]) 
        @transform(:Dₕ = sum(exp.(:δjn / (1-σ)))) # Denominator
        @transform(:s_within = exp.(:δjn / (1-σ)) ./ :Dₕ)     # within group share
        groupby([:market, :group])
        @transform(:s_group = :Dₕ.^(1-σ) ./ sum(unique(:Dₕ.^(1-σ))) )  # group share
        @transform(:s_market = :s_within .* :s_group)
        _[!, [:market, :group, :δjn, :δ, :Dₕ, :s_within, :s_group, :s_market, :sjn]]
    end

    @chain df1 begin
        # Add within group share Eq 23
        groupby([:market, :group])
        combine(:, [:δ,:Dₕ] => (δ,D) -> δ ./ D => :sjhn_hat)
        # Add group share Eq 24
        groupby([:market, :group])
        combine(:, :Dₕ => (D -> D.^(1-σ) ./ sum(D.^(1-σ)) => :sjhn_hat))
        # Add market share Eq 25
        groupby([:market, :group])
        combine(:, :Dₕ => (D -> D.^(1-σ) ./ sum(D.^(1-σ)) => :sjhn_hat))
    end




    @chain df1 begin
        # Add denominator
        groupby([:market, :group])
        combine(:, :δ => (δ -> sum(exp.(δ / (1-σ)))) => :Dₕ)
        # Add within group share Eq 23
        groupby([:market, :group])
        combine(:, [:δ,:Dₕ] => (δ,D) -> δ ./ D => :sjhn_hat)
        # Add group share Eq 24
        groupby([:market, :group])
        combine(:, :Dₕ => (D -> D.^(1-σ) ./ sum(D.^(1-σ)) => :sjhn_hat))
        # Add market share Eq 25
        groupby([:market, :group])
        combine(:, :Dₕ => (D -> D.^(1-σ) ./ sum(D.^(1-σ)) => :sjhn_hat))
    end
    # Estimate within-group share using Eq (23)
    df.sjhn_hat = @chain df1 begin
        groupby([:market, :group])
        combine(:δ => (δ -> eq23(δ, σ)) => :sjn_hat)
        _[!,:sjn_hat]
    end
    return df_
end











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



