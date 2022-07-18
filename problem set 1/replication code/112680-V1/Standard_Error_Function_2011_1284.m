function output = Standard_Error_Function_2011_1284(Q,alpha,R,EPS,nIs,nPlans,nTs,HTC,Sim,K,EUvector11,EUvector21,EUvector31,EUvector12,EUvector22,EUvector32,EUvector13,EUvector23,EUvector33,choice,choice12,choice22,choice32,choice13,choice23,choice33,Income,Tier2,IND,FSAY2,FSAY3,QS,managerX,CC1,CC2,CSAL2,CSAL3,MAge)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   Sample Code to Compute Std.      %%%%%%%%%%%%
%%%%%%%%%%%%   Error for One Paramater in       %%%%%%%%%%%%
%%%%%%%%%%%%   Choice Model	 					%%%%%%%%%%%%
%%%%%%%%%%%%   								    %%%%%%%%%%%%
%%%%%%%%%%%%   This is function called from     %%%%%%%%%%%%
%%%%%%%%%%%%   primary Standard Error file		%%%%%%%%%%%%
%%%%%%%%%%%%   ‘Standard_Error_OneParam_2011_1284.m’ %%%%%%%
%%%%%%%%%%%%									%%%%%%%%%%%%
%%%%%%%%%%%%   Computes likelihood ratio test 	%%%%%%%%%%%%
%%%%%%%%%%%%   statistic which is put into 		%%%%%%%%%%%%
%%%%%%%%%%%%   optimization routine in main		%%%%%%%%%%%%
%%%%%%%%%%%%   called in main file. Finds 		%%%%%%%%%%%%
%%%%%%%%%%%%   parameter value for given param	%%%%%%%%%%%%
%%%%%%%%%%%%   that implies ones std. devitation  %%%%%%%%%%
%%%%%%%%%%%%   away in LR test for actual estimates. %%%%%%%
%%%%%%%%%%%%   Details described at end of main SE file%%%%%
%%%%%%%%%%%%   										 %%%%%%%
%%%%%%%%%%%%   Ben Handel                       %%%%%%%%%%%%
%%%%%%%%%%%%   handel@berkeley.edu              %%%%%%%%%%%%
%%%%%%%%%%%%									%%%%%%%%%%%%
%%%%%%%%%%%%   Adverse Selection and Inertia    %%%%%%%%%%%%
%%%%%%%%%%%%   in Health Insurane Markets:      %%%%%%%%%%%%
%%%%%%%%%%%%   When Nudging Hurts               %%%%%%%%%%%%
%%%%%%%%%%%%                                    %%%%%%%%%%%%
%%%%%%%%%%%%   March, 2013               	    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Standard error function is called from file ‘Standard_Error_OneParam_2011_1284.m’ provided. 
% This function inverts a standard likelihood ratio test to find a simple way to find the 
% standard errors in this simulated maximum likelihood model. First the code is given / described then 
% the logic / process is outlined in notes below the likelihood function. 

% alpha represents vector of paramater estimates that standard errors are being computed for
% This is an input into this optimization from the file ‘Standard_Error_OneParam_2011_1284.m’
% that calls this function 

% alphaX is a new set of parameters where all but 1 parameters are the same as the estimated values 
% and 1 value (the one we're finding the standard error 4, equals Q which is different than the estimated value.
% The non-linear optimization searches for the Q (value of this parameter) that sets that implies the likelihood 
% ratio test statistic value for that parameter is one standard deviation away from the estimted value. This 
% is then the standard error. The logic is described in more detail below at the end of ‘Standard_Error_OneParam_2011_1284.m’  

alphaX = alpha;
alphaX(1,1) = Q;

% Equation is (Log-Likelihood at actual values Alpha - Log-Likelihood with AlphaX) + 0.5
% We want to set this equation = 0, so we square the equation and minimize it to find the value of Q that sets it closest to 0. 

output = (Likelihood_2011_1284(alpha,R,EPS,nIs,nPlans,nTs,HTC,Sim,K,EUvector11,EUvector21,EUvector31,EUvector12,EUvector22,EUvector32,EUvector13,EUvector23,EUvector33,choice,choice12,choice22,choice32,choice13,choice23,choice33,Income,Tier2,IND,FSAY2,FSAY3,QS,managerX,CC1,CC2,CSAL2,CSAL3,MAge) - Likelihood_2011_1284(alphaX,R,EPS,nIs,nPlans,nTs,HTC,Sim,K,EUvector11,EUvector21,EUvector31,EUvector12,EUvector22,EUvector32,EUvector13,EUvector23,EUvector33,choice,choice12,choice22,choice32,choice13,choice23,choice33,Income,Tier2,IND,FSAY2,FSAY3,QS,managerX,CC1,CC2,CSAL2,CSAL3,MAge) + 0.5)^2



