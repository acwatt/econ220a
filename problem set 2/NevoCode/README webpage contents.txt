Computer Code

Attached below is a Matlab script file and eight accompanying Matlab functions that compute the random coefficients discrete choice model described in "A Practitioner's Guide to Estimation of Random Coefficients Logit Models of Demand" (Journal of Economics & Management Strategy, 9(4), 513-548, 2000) and "Measuring Market Power in the Ready-to-Eat Cereal Industry" (Econometrica, 69(2), 307-342, 2001).

The code is provided for academic research. Users of this code (or a modified version of it) should reference the above papers. The code includes minimal documentation, and is provided with out any additional support. The code was successfully run using Matlab 5.1 on a Windows NT operating system. Questions regarding Matlab should be referred to MathWorks.

The program consists of the following files (all Matlab m-files):
rc_dc.m - A script file that reads in the data and calls the other functions;
gmmobj.m - This function computes the GMM objective function;
meanval.m - This function computes the mean utility level;
mufunc.m - This function computes the non-linear part of the utility (mu_ijt in the Guide);
mktsh.m - This function computes the market share for each product;
ind_sh.m - This function computes the "individual" probabilities of choosing each brand;
gradobj.m - This function computes the gradient of the objective function;
jacob.m - This function computes the Jacobian of the implicit function that defines the mean utility;
var_cov.m - This function computes the VCov matrix of the estimates.
cr_dum.m - This function creates a set of dummy variables