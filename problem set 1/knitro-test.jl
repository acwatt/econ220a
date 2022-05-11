using KNITRO
function callbackEvalF(kc, cb, evalRequest, evalResult, userParams)
    x = evalRequest.x
    # Evaluate nonlinear objective
    dTmp = x[2] - x[1] * x[1]
    evalResult.obj[1] = 100.0 * (dTmp * dTmp) + ((1.0 - x[1]) * (1.0 - x[1]))

    return 0
end

function callbackEvalG!(kc, cb, evalRequest, evalResult, userParams)
    x = evalRequest.x

    # Evaluate gradient of nonlinear objective
    dTmp = x[2] - x[1] * x[1]
    evalResult.objGrad[1] = (-400.0 * dTmp * x[1]) - (2.0 * (1.0 - x[1]))
    evalResult.objGrad[2] = 200.0 * dTmp

    return 0
end

# Create a new Knitro solver instance.
kc = KNITRO.KN_new()

n = 2
KNITRO.KN_add_vars(kc, n)
KNITRO.KN_set_var_lobnds(kc,  [-KNITRO.KN_INFINITY, -KNITRO.KN_INFINITY]) # not necessary since infinite
KNITRO.KN_set_var_upbnds(kc,  [0.5, KNITRO.KN_INFINITY])
# Define an initial point.  If not set, Knitro will generate one.
KNITRO.KN_set_var_primal_init_values(kc, [-2.0, 1.0])

# Add the constraints and set their lower bounds
m = 2
KNITRO.KN_add_cons(kc, m)
KNITRO.KN_set_con_lobnds(kc, [1.0, 0.0])

# First load quadratic structure x0*x1 for the first constraint
KNITRO.KN_add_con_quadratic_struct(kc, 0, 0, 1, 1.0)

# Add linear term x0 in the second constraint
KNITRO.KN_add_con_linear_struct(kc, 1, 0, 1.0)

# Add quadratic term x1^2 in the second constraint
KNITRO.KN_add_con_quadratic_struct(kc, 1, 1, 1, 1.0)

# Set verbose printing level
KNITRO.KN_set_param(kc, "outlev", 6)

cb = KNITRO.KN_add_eval_callback_all(kc, callbackEvalF)
KNITRO.KN_set_cb_grad(kc, cb, callbackEvalG!)
nStatus = KNITRO.KN_solve(kc)

# Delete the Knitro solver instance.
KNITRO.KN_free(kc)