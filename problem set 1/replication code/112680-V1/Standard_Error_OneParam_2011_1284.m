
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   Sample Code to Compute Std.      %%%%%%%%%%%%
%%%%%%%%%%%%   Error for One Paramater in       %%%%%%%%%%%%
%%%%%%%%%%%%   Choice Model	 					%%%%%%%%%%%%
%%%%%%%%%%%%   								    %%%%%%%%%%%%
%%%%%%%%%%%%  As discussed in README file, this %%%%%%%%%%%%
%%%%%%%%%%%%  code provides example of how to   %%%%%%%%%%%%
%%%%%%%%%%%%  compute standard errors for one 	%%%%%%%%%%%%
%%%%%%%%%%%%  parameter. The code can be repeated %%%%%%%%%%
%%%%%%%%%%%%  for each paramater separately  	%%%%%%%%%%%%
%%%%%%%%%%%%  to compute standard errors for all %%%%%%%%%%%
%%%%%%%%%%%%  parameters as is done in paper. 	%%%%%%%%%%%%
%%%%%%%%%%%%  As discussed in README / paper 	%%%%%%%%%%%%
%%%%%%%%%%%%  this is just one possible method 	%%%%%%%%%%%%
%%%%%%%%%%%%  to compute SE and it only estimates %%%%%%%%%%
%%%%%%%%%%%%  diagonal elements in the variance	%%%%%%%%%%%%
%%%%%%%%%%%%  covariance matrix, which appear in  %%%%%%%%%%
%%%%%%%%%%%%  table 5 in the paper. For tests 	%%%%%%%%%%%%
%%%%%%%%%%%%  like overall goodness-of-fit or  	%%%%%%%%%%%%
%%%%%%%%%%%%  multi-parameter significance  	%%%%%%%%%%%%
%%%%%%%%%%%%  the whole matrix would have to 	%%%%%%%%%%%%
%%%%%%%%%%%%  be computed using alternative   	%%%%%%%%%%%%
%%%%%%%%%%%%  methods. 						  	%%%%%%%%%%%%
%%%%%%%%%%%%									%%%%%%%%%%%%
%%%%%%%%%%%%   Ben Handel                       %%%%%%%%%%%%
%%%%%%%%%%%%   handel@berkeley.edu              %%%%%%%%%%%%
%%%%%%%%%%%%									%%%%%%%%%%%%
%%%%%%%%%%%%   Adverse Selection and Inertia    %%%%%%%%%%%%
%%%%%%%%%%%%   in Health Insurane Markets:      %%%%%%%%%%%%
%%%%%%%%%%%%   When Nudging Hurts               %%%%%%%%%%%%
%%%%%%%%%%%%                                    %%%%%%%%%%%%
%%%%%%%%%%%%   March, 2013               	    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% FIRST part of code repeats first part of code from 'EstimationCode_2011_1284.m'
%%%%%%%%%%% This simply processes simulated data into version used for estimation and, here, 
%%%%%%%%%%% used to compute SEs for those estimates. New Code for SEs is at end of file, starting on line 209. 

load 'ASIN-ChoiceModelData-FINAL_2011_1284.mat'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% See corresponding data desription document for detailed description  %%%%%%%%%%
%%%%%%%%%% of variables in simulated data. 										%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Generate standard normal matrices that epsilon standard    %%
%%%%%%%% deviation parameters will scale up or down in estimation   %%
%%%%%%%% Simulations generated in advance of likelihood function    %%
%%%%%%%% so draws are fixed throughout estimation. 					%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

muE = [0,0,0,0,0,0];
sigmaE = eye(6) %Identity matrix
EPS = mvnrnd(muE, sigmaE, Sim); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Generate standard normal vectors that risk preference and PPO1200	%%%%%%%
%%%%%%%% random coefficient standard deviations will scale up or down		%%%%%%%
%%%%%%%% NOTE: these will just impct the variances of risk preferences, 	%%%%%%%
%%%%%%%% means will be added later in the likelihood file					%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu = [0,0];
sigma = eye(2);
R = mvnrnd(mu, sigma, Sim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Put out of pocket simulated draw matrices in form to quickly add to expected utility calculations		%%%%%
%%%%%%%%%%%%% Add third dimension on random coefficient simulations in order to expedite calculations				%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PPO250OP1 = zeros(nIs,K,Sim);
PPO250OP2 = zeros(nIs,K,Sim);
PPO250OP3 = zeros(nIs,K,Sim);

PPO500OP1 = zeros(nIs,K,Sim);
PPO500OP2 = zeros(nIs,K,Sim);
PPO500OP3 = zeros(nIs,K,Sim);

PPO1200OP1 = zeros(nIs,K,Sim);
PPO1200OP2 = zeros(nIs,K,Sim);
PPO1200OP3 = zeros(nIs,K,Sim);

%%%%%%%%%%%%%%%% Fill in OOP Claims by Person and Year and Simulation %%%%%

for i = 1:Sim
    PPO250OP1(:,:,i) = PPO250OOP1;
    PPO250OP2(:,:,i) = PPO250OOP2;
    PPO250OP3(:,:,i) = PPO250OOP3;

    PPO500OP1(:,:,i) = PPO500OOP1;
    PPO500OP2(:,:,i) = PPO500OOP2;
    PPO500OP3(:,:,i) = PPO500OOP3;
    
    PPO1200OP1(:,:,i) = PPO1200OOP1;
    PPO1200OP2(:,:,i) = PPO1200OOP2;
    PPO1200OP3(:,:,i) = PPO1200OOP3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% High Total Claims Indicator: Set up earlier in simulation code. When 		%%%%%%%
%%%%%%%%%%%%%%% People have high mean spending in a year, estimate additional preference 	%%%%%%%
%%%%%%%%%%%%%%% factor for PPO250 vs. other two plans here. This part just translates 		%%%%%%%
%%%%%%%%%%%%%%% Indicator for each family from simulated data into matrix easy to use 		%%%%%%%
%%%%%%%%%%%%%%% in likelihood function to expedite calculation. 							%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HTC = zeros(nIs,K,Sim);

for k = 1:K
    for s = 1:Sim
    HTC(:,k,s) = HTCi;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Compute money at stake in each state (realization of health status ex		%%%%%%%
%%%%%% post) ahead of estimation so this doesn't have to be computed in 			%%%%%%%
%%%%%% likelihood function loop. Saves a lot of time. These will be inputs			%%%%%%%
%%%%%% into likelihood function in estimation. 										%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EUvector11 = zeros(nIs,K,Sim);
EUvector21 = zeros(nIs,K,Sim);
EUvector31 = zeros(nIs,K,Sim);
EUvector12 = zeros(nIs,K,Sim);
EUvector22 = zeros(nIs,K,Sim);
EUvector32 = zeros(nIs,K,Sim);
EUvector13 = zeros(nIs,K,Sim);
EUvector23 = zeros(nIs,K,Sim);
EUvector33 = zeros(nIs,K,Sim);

%%%% Wealth: doesn't matter since preferences are assumed to be CARA 

W = 75000;

%%% max is here to protect against extremely high draw, since CARA can't have negative $, could also just set wealth higher

for i = 1:nIs
    for k = 1:K
        EUvector11(i,k,:) = max(0.1,(W - reshape(PPO250OP1(i,k,:),1,Sim) - reshape(P1(1,i,:),1,Sim)));
        EUvector21(i,k,:) = max(0.1,(W - reshape(PPO500OP1(i,k,:),1,Sim) - reshape(P2(1,i,:),1,Sim)));
        EUvector31(i,k,:) = max(0.1,(W - reshape(PPO1200OP1(i,k,:),1,Sim) - reshape(P3(1,i,:),1,Sim)));
        EUvector12(i,k,:) = max(0.1,(W - reshape(PPO250OP2(i,k,:),1,Sim) - reshape(P1(2,i,:),1,Sim)));
        EUvector22(i,k,:) = max(0.1,(W - reshape(PPO500OP2(i,k,:),1,Sim) - reshape(P2(2,i,:),1,Sim)));
        EUvector32(i,k,:) = max(0.1,(W - reshape(PPO1200OP2(i,k,:),1,Sim) - reshape(P3(2,i,:),1,Sim)));
        EUvector13(i,k,:) = max(0.1,(W - reshape(PPO250OP3(i,k,:),1,Sim) - reshape(P1(3,i,:),1,Sim)));
        EUvector23(i,k,:) = max(0.1,(W - reshape(PPO500OP3(i,k,:),1,Sim) - reshape(P2(3,i,:),1,Sim)));
        EUvector33(i,k,:) = max(0.1,(W - reshape(PPO1200OP3(i,k,:),1,Sim) - reshape(P3(3,i,:),1,Sim)));
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Develop indicators for incumbent plans for years 2 and 3 so that we know ahead of estimation 			%%%
%%%%%%%%%% which plans are incumbent plans for each person in each year. This is just determined by 'data'			%%%
%%%%%%%%%% here, of course, this is the simulated data																%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

choice12 = zeros(nIs,Sim);
choice22 = zeros(nIs,Sim);
choice32 = zeros(nIs,Sim);
choice13 = zeros(nIs,Sim);
choice23 = zeros(nIs,Sim);
choice33 = zeros(nIs,Sim);

for i = 1:nIs
    for s = 1:Sim
       
	   if choice(i,1)==1
        choice12(i,s) = 1;
       end
       if choice(i,1)==2
        choice22(i,s) = 1;
       end
       if choice(i,1)==3
        choice32(i,s) = 1;
       end

       if choice(i,2)== 1
        choice13(i,s) = 1;
       end
       if choice(i,2)== 2
        choice23(i,s) = 1;
       end
       if choice(i,2)==3
        choice33(i,s) = 1;
       end
       
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Set income equal to year 2 income. Done since random coefficient for risk preferences, 					%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% which depends on income, is constant over time in panel within person by assumption.  					%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Income = Inc2;

%%%%% Set age equal to maximum age in family

MAge = max(Ages,[],2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Computing Standard Errors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% This file / computation shows how to compute standard error for one parameter from estimation. 
%%%% In actual analysis, SEs in table 5 were calculated one parameter at a time using this routine for each paramater. 
%%%% See notes at beginning of this file for extended discussion of this, and how alternative methods would be needed to 
%%%% test multiple paramaters at once. This is in essence a 'shortcut' to easily get the diagonals of the parameter
%%%% variance / covariance matrix. 
%%%%
%%%% Here, we illustrate methodology for alpha(1,1), which is mean of PPO1200 / CDHP random coefficient distribution for families. 
%%%%
%%%% NOTE: Detailed description of logic behind the SE computation that involves this file and the function called below is given at the 
%%%% end of this file. See those notes below for logical behind these one-paramater standard error calculations. Of course, feel 
%%%% free to contact author Ben Handel at handel@berkeley.edu with further questions / clarifications. 

options = optimset('MaxFunEvals',40000,'MaxIter',40000);

%%%%%% Alpha here represents actual parameter estiamtes of choice model, output of file 'EstimationCode_2011_1284.m'

alpha = [-5344,-8.97,0.07,1.37; 2179,0,3006,0;676,90,0,-1386;-2905,989,0,0;2430,161,50,-723;-8,-537,875,-221;61,0.28,-327,0;]

%%%%% Q0 is initial guess for value of alpha(1,1) which is one standard deviation away from true estimated value. Non-linear optimizer will 
%%%%% solve for Q which is exactly one standard deviation away. Then, abs(Q - alpha(1,1)) is standard error for that parameter. This computation can 
%%%%% then be repeated for each parameter. 

Q0 = *********;

[qsol,fval] = ktrlink( @(Q) Standard_Error_Function_2011_1284(Q,alpha,R,EPS,nIs,nPlans,nTs,HTC,Sim,K,EUvector11,EUvector21,EUvector31,EUvector12,EUvector22,EUvector32,EUvector13,EUvector23,EUvector33,choice,choice12,choice22,choice32,choice13,choice23,choice33,Income,Tier2,IND,FSAY2,FSAY3,QS,managerX,CC1,CC2,CSAL2,CSAL3,MAge), Q0,[],[],[],[],[],[],[],options);

qsol

SE_CDHPMEAN_FAMILY = abs(qsol-alpha(1,1)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Description of underlying logic (examine together with function called above and this main SE file):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Likelihood ratio statistic for accepting / rejecting null is D = -2(LL(THETA*)-LL(THETA**))
%
% LL(THETA*) is likelihood function at estimated parameter values, null
% LL(THETA**) is alternative hypothesis alphaX (changing one paramater) likelihood 
% 
% 2. D is asymptotically distributed chi-squared with 1 degree of freedom
%
% 3. Distribution of paramater we're testing is asymptotically normal, so one standard deviation in this 
% paramater corresponds to approx. .33 P-Value point in normal (similar to how two Std. Dev used for two-sided .05 test)
%
% 4. Translated to CHI-SQUARED TEST, .33 P-Value occurs at approx. D = 1 for test statistic. 
%
% 5. Given this, one standard deviation movement to Q for paramater in question (alpha(1,1) is example here) from original value implies equation:
%
% 1= -2(LL(THETA*)-LL(THETA**))
%
% 6. Solving this equation for Q in THETA** and then taking absolute value of (Q-alpha(1,1)), or whatever parameter you're doing SE for, computes the standard error. 
% 
% 7. This solving is implemented in the main SE file  ‘Standard_Error_OneParam_2011_1284.m’ which calls this function to find the right value of Q with non-linear optimization. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% END OF FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


















            
            
            