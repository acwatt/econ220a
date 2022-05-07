




#==============================================================================
                               HELPER FUNCTIONS
==============================================================================#
summary_file = "summary_stats.csv"
summary_file = string(dir, summary_file)
function append_summary(df, title; replace=false)
    # Append title row to CSV
    if !isfile(summary_file) | replace  # Write new file if no file exists, or `replace` is true
        CSV.write(summary_file, DataFrame(title=title), writeheader=false)
    else
        CSV.write(summary_file, DataFrame(title=["", title]), writeheader=false, append=true)
    end
    # Append dataframe to CSV
    CSV.write(summary_file, df, append=true, writeheader=true)
end






#==============================================================================
#                        QUESTION 1: SUMMARY STATISTICS
==============================================================================#

#==========================================================
PART A: EXPENDITURES

For each family calculate the mean and standard deviation of
their out-of-pocket expenditures for each of the three health plans
in year 1.
==========================================================#
# What is the population mean of out-of-pocket expenditures in each plan?
title = "Year 1 population mean of all draws of out-of-pocket expenditures, by plan"
plans = ["250", "500", "1200"]; means = []; counts = []
for plan ∈ plans
    var = eval(Symbol(string("PPO", plan, "OOP1")))
    println("$title $plan: $(mean(var))")
    append!(means, mean(var))
    append!(counts, mean(var))
end
df = DataFrame(Plan=plans, Mean=means)
append_summary(df, title; replace=true)

# What is the population mean of the standard deviation of out-of-pocket 
# expenditures in each plan?
title = "Year 1 population mean of standard deviation of all draws of out-of-pocket expenditures, by plan"
sds = []
for plan ∈ plans
    var = eval(Symbol(string("PPO", plan, "OOP1")))
    println("$title $plan: $(mean(std(var, dims=2)))")
    append!(sds, mean(var))
end
df = DataFrame(Plan=plans, Mean=means)
append_summary(df, title)

# Provide the above two statistics,  broken down by family 
# status in year 1, given in the variable ’Tier1’. 
#
# Since all families have exactly K draws for each year, each plan
# I can first save the mean of each family over all K draws,
# then groupby family status and take the mean of each group
# i.e., the subset mean is the mean of the means within the subset
title = "Year 1 population mean of out-of-pocket expenditures, by plan, by family status"
tier_lookup = DataFrame("tier"=>[1,2,6,8,11,12], "Family Status"=>["Single", "With Spouse", "With Children", "With Spouse & Children", "Unknown (11)", "Unknown (12)"])
df_list = []
# convert Tier 
for plan ∈ ["250", "500", "1200"]
    var = eval(Symbol(string("PPO", plan, "OOP1")))
    families = DataFrame(mean = vec(mean(var, dims=2)),
                         tier = vec(Tier1),
                         std = vec(std(var, dims=2)))
    families = leftjoin(families, tier_lookup, on=:tier)
    gd = groupby(families, :"Family Status")
    println("\n$title $plan and family status:")
    for key ∈ keys(gd)
        mean_ = round(mean(gd[key].mean), digits=3)
        std_ = round(mean(gd[key].std), digits=3)
        println("$(key[1]): \t$mean_ \t$std_")
        # Create dataframe of summary information
        local df = DataFrame("Plan"=>plan, "Family Status"=>key[1], "Mean"=>mean_, "Std Dev"=>std_, "Count"=>nrow(gd[key]))
        push!(df_list, df)
    end
end
df = reduce(vcat, df_list)
append_summary(df, title)









#==========================================================
PART B: TABULATIONS

For year 1 provide tabulations for the following variables across
the population of families:
    - Indicator of chronic conditions in family
    - Enrollment in Flexible Spending Account
    - Employee Manager Status
    - Family Status
    - Income
    - Ages (do this for every individual in the population, use the
      matrix ’Ages’)
    - Genders (do this for every individual in the population, use
      the matrix ’Genders’)
==========================================================#
df_list = []
# Build the dataframe
df = DataFrame("Chronic condition in family"=>vec(CC1),
               "Enrollment in Flexible Spending Account"=>vec(FSAY1),
               "Employee Manager Status"=>vec(vec(managerX)),
               "tier"=>vec(Tier1),
               "Income"=>vec(Inc1),
               "Family Size"=>vec(Famsize))
# Merge family status
df = leftjoin(df, tier_lookup, on="tier")
df = select(df, Not("tier"))
# Summarize Family Status
df_family = combine(groupby(df, "Family Status"), "Family Status" => length => :Count)
append_summary(df_family, "Family Status Tabulation")

# Summarize data
df = describe(df, 
               mean => :Mean,
               std => "Std Dev",
               median => :Median,
               minimum => :Min,
               maximum => :Max,
               length => :Count,
               cols=Not("Family Status"))
df = rename(df, :variable => :Variable)
df[!,:Variable] = String.(df[!,:Variable])

# Add ages
Ages_ = skipmissing(replace(Ages, 0 => missing))
df_ = DataFrame("Variable" => "Age (by individual)",
               "Mean" => mean(Ages_),
               "Std Dev" => std(Ages_),
               "Median" => median(Ages_),
               "Min" => minimum(Ages_),
               "Max" => maximum(Ages_),
               "Count" => sum(Ages .!= 0))
df = append!(df, df_)

# Add Genders
# First need to add missing values into Genders using Ages 0 values
Genders_ = skipmissing((replace(Ages .== 0, 1=>missing) + Genders))
df_ = DataFrame("Variable" => "Gender (by individual)",
                "Mean" => mean(Genders_),
                "Std Dev" => std(Genders_),
                "Median" => median(Genders_),
                "Min" => minimum(Genders_),
                "Max" => maximum(Genders_),
                "Count" => sum(Ages .!= 0))
df = append!(df, df_)
append_summary(df, "Family and Individual Tabulations")



