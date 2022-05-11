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
lb = [
          0.0     -0.5     -0.5      0.0
        -_inf       0.1    -_inf       0.1
          1.0      1.0      1.0      1.0
        -_inf       0.0      0.0  -5000.0
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
        _inf      _inf      _inf      _inf
    20000.0  20000.0  20000.0  20000.0
        _inf   20000.0  20000.0   5000.0
    5000.0   5000.0   5000.0   5000.0
    5000.0      0.0      0.0      0.0
]

KNITRO.KN_set_var_lobnds(kc,  transform_to_21vec(lb)) # not necessary since infinite
KNITRO.KN_set_var_upbnds(kc,  transform_to_21vec(ub))
# Define an initial point.  If not set, Knitro will generate one.
KNITRO.KN_set_var_primal_init_values(kc, transform_to_21vec(α₀))


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







#= MULTI-START info
The multi-start procedure generates new start points by 
randomly selecting components of x that satisfy lower 
and upper bounds on the variables. Knitro finds a local 
optimum from each start point using the same problem 
definition and user options. The final solution returned 
from KN_solve() is the local optimum with the best 
objective function value if any local optima have been 
found. If no local optimum has been found, Knitro 
will return the best feasible solution estimate it 
found. If no feasible solution estimate has been 
found, Knitro will return the least infeasible point.
=#


#= Parallel multi-start
The multi-start procedure can run in parallel on 
shared memory multi-processor machines by setting 
numthreads greater than 1. See Parallelism for more 
details on controlling parallel performance in Knitro.
https://www.artelys.com/docs/knitro/2_userGuide/parallelism.html#sec-parallelism

When the multi-start procedure is run in parallel, 
Knitro will produce the same sequence of initial 
points and solves that you see when running multi-start 
sequentially (though, perhaps, not in the same order).
=#






#=
Final Statistics
----------------
Final objective value               =   7.50199018429976e+03
Final feasibility error (abs / rel) =   0.00e+00 / 0.00e+00
Final optimality error  (abs / rel) =   1.34e-02 / 1.34e-02
# of iterations                     =          7 
# of CG iterations                  =          2 
# of function evaluations           =        222
# of gradient evaluations           =          0
Total program time (secs)           =     475.52383 (   471.476 CPU time)

Solution Vector
---------------
x[       0] =   1.35022844344e-01,   lambda[       0] =   1.84719518646e+01
x[       1] =  -2.50000000000e+03,   lambda[       1] =   0.00000000000e+00
x[       2] =   3.00000201855e+02,   lambda[       2] =  -1.89076078593e-04
x[       3] =  -5.00000000000e+02,   lambda[       3] =   0.00000000000e+00
x[       4] =   3.59395777118e-03,   lambda[       4] =   2.46801365248e-07
x[       5] =   4.01293749016e-04,   lambda[       5] =   2.72315495294e-08
x[       6] =  -9.60312710392e-02,   lambda[       6] =   5.07245277320e+01
x[       7] =   7.00000108357e+02,   lambda[       7] =  -0.00000000000e+00
x[       8] =   7.99999930924e+02,   lambda[       8] =  -7.25860578923e-05
x[       9] =   1.25000063970e+03,   lambda[       9] =   0.00000000000e+00
x[      10] =   1.69105891691e-04,   lambda[      10] =   1.09672887426e-08
x[      11] =  -2.43233375477e-01,   lambda[      11] =   3.92739822611e+01
x[      12] =  -2.20000000000e+03,   lambda[      12] =   0.00000000000e+00
x[      13] =   3.00000331533e+02,   lambda[      13] =   0.00000000000e+00
x[      14] =   1.75000084943e+03,   lambda[      14] =   0.00000000000e+00
x[      15] =   4.31119323917e-04,   lambda[      15] =   2.85357084218e-08
x[      16] =   7.32597088355e-01,   lambda[      16] =   2.84982144363e+01
x[      17] =   8.00000094812e+02,   lambda[      17] =  -0.00000000000e+00
x[      18] =   8.00000131963e+02,   lambda[      18] =  -5.27681371714e-05
x[      19] =  -4.99999495048e+02,   lambda[      19] =   0.00000000000e+00
x[      20] =   6.37050985741e-04,   lambda[      20] =   4.15483482748e-08

===============================================================================

=#