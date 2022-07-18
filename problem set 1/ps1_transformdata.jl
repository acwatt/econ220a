#==============================================================================
file: ps1_transformdata.jl
description: steps 1.c and 1.d for problem set. Must be run inside ps1.jl.
author: Aaron C Watt (UCB Grad Student, Ag & Resource Econ)
==============================================================================#



#==========================================================
PART C: TRANSFORM THE DATA

Now, let’s transform / add to the data to get it ready for choice
model estimation. To do this start a code file called ‘Estimation-
Code.m’ in Matlab and start by loading the data.
    - Generate a 6 X Sim matrix of standard normal random draws.
      These will be used to index the family-plan-time specific ep-
      silons for two plans over three years in estimation (the draws
      multiplied by the candidate standard deviations will give a
      normal distribution with that standard deviation.)
    - Generate two additional 1 X Sim matrices with standard nor-
      mal draws. These will index the risk preference and PPO1200
      normally distributed random coefficients.
    - The out-of-pocket variables in the simulated data are in nIs X
      K matrices. Create (nIs, K, Sim) three dimensional matrices
      that extend the data by adding a third dimension with Sim
      columns (this repeats the same two-dimensional matrix for
      each family Sim times). This will be useful to speed up your
      calculations.
    - Take the nIs X 1 HTCi variable and put it into an (nIs, K,
      Sim) matrix (you repeat the HTCi matrix many times)
    - Set a ’permanent’ income variable equal to the income in year
      2 of the data (nIs,1).
    - Set a family ’age’ variable equal to the maximum age in each
      family (nIs, 1).
==========================================================#
# 6xSim matrix of Std. Normal draws, for family-plan-time 
# specific epsilons for two plans over three years in estimation
ϵ = randn((6, Sim))

# 1xSim risk preference and PPO1200 normally distributed random coefficients
randcoef_risk = randn((1,Sim))
randcoef_1200 = randn((1,Sim))

# Dictionary to store new matrices
mat_dict = Dict()

# nIs x K x Sim matrix of OOP costs for each family (repeated Sim times)
for year ∈ 1:3
    for plan ∈ ["250", "500", "1200"]
        # Get original matrix from the imported Matlab dictionary
        varname = "PPO$(plan)OOP$year"
        old = vars[varname]  # nIs x K
        # Create new matrix, copying columns Sim times into a 3rd dimension (nIs x K x Sim)
        new = reshape(repeat(old, inner=(1,Sim,1)), nIs, K, Sim)
        # Assign the new matrix to the varname in memory
        @eval (($(Symbol(varname))) = ($new))
        mat_dict[varname] = new
    end
end

# Copy nIs x 1 HTCi into nIs x K x Sim
HTCi = repeat(HTCi, outer=(1, K, Sim))

# Permanent Income set to year 2 income (nIs,1)
Inc_perm = Inc2

# Max age in each family (nIs,1)
Ages_max = maximum(Ages, dims=2)




#==========================================================
PART D: COMPUTE QUANTITIES NECESSARY FOR CHOICE MODEL

There are two final steps to prepare the data for estimation.
These steps compute quantities necessary for choice model es-
timation, but that are ’fixed’ throughout estimation. This means
it is better to compute them up front and use them as inputs into
estimation.
    - Define a wealth variable equal to 75,000 (just a scalar).
    - Generate nine matrices (nIs,K,Sim) that describe the basic
      ’money at stake’ for each person and each potential health
      state realization and each year. Here, a person-state is on the
      nIs,K dimensions while we repeat this Sim times for ease in
      estimation. The nine matrices are for each of the three plans
      over each of the three years. Each entry (i,k,s) in year t for
      plan p should be (Wealth - OOP Expense draw (i,k) (for plan
      p in year t) - Price (i) (for plan p year t)). NOTE: prices
      depend on family status and are entered accordingly into
      price matrices already. Label these matrices EUvector‘p’‘t’
      where p is for plan and ‘t’ is for time.
    - Define six (nIs, Sim) matrices that define the actual default
      plans for families in non-active choice years. Label each ma-
      trix choice‘p’‘t’. Here, year 1 is an active choice year, so t
      should only be 2 or 3. Use the observed choices to fill in these
      matrices with binary values: if a family chose plan ‘p’ in year
      ’t-1’ then choice‘p’‘t’ = 1. This leads to (nIs,1) matrices that
      you should repeat Sim times on second dimension to use in
      estimation.
==========================================================#

wealth = 75000

# Money at stake, each person, year, health state realization
price_map = Dict("250"=>P1, "500"=>P2, "1200"=>P3);
for plan ∈ plans
    price_mat = price_map[plan]
    for year in 1:3
        # Get out of pocket expenses (already copied Sim times above to get nIs x K x Sim)
        OOP = mat_dict["PPO$(plan)OOP$year"]
        # Get plan-year prices for each family (nIs x Sim)
        prices = price_mat[year, :, :]
        # Repeat prices K times (prices don't change based on health risk or risk preference)
        prices = repeat(prices, outer=(1, 1, K))
        # Wealth - OOP - plan price (.- is element-wise subtraction)
        mat_dict["EUvector$plan$year"] = wealth .- OOP .- prices
        # This changes by-family (nIs row) and by-health-simulation (K matrix stack), but not by-risk-simulation (Sim column)
    end
end


# Default plans for families in non-active choice years (nIs x Sim)
plan_dict = Dict("250" => 1, "500" => 2, "1200" => 3)
for plan ∈ plans
    for year in 2:3  # only non-active choice years
        vec = Int.(choice[:, year] .== plan_dict[plan])
        mat_dict["choice$plan$year"] = repeat(vec, outer=(1, Sim))
    end
end


# Add new matricies to variables in memory
dict_to_mem(mat_dict)