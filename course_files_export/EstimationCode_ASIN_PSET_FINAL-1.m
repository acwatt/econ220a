%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   Estimation Code for Primary      %%%%%%%%%%%%
%%%%%%%%%%%%   Choice Model Estimation Using    %%%%%%%%%%%%
%%%%%%%%%%%%   Simulated Data 					%%%%%%%%%%%%
%%%%%%%%%%%%   								    %%%%%%%%%%%%
%%%%%%%%%%%%  (PRIMARY CHOICE MODEL-ESTIMATION) %%%%%%%%%%%%
%%%%%%%%%%%%   						            %%%%%%%%%%%%
%%%%%%%%%%%%  Data are simulated with code that %%%%%%%%%%%%
%%%%%%%%%%%%  is available from Benjamin Handel %%%%%%%%%%%% 
%%%%%%%%%%%%  upon request by emailing at		%%%%%%%%%%%%
%%%%%%%%%%%%  handel@berkeley.edu. 				%%%%%%%%%%%%
%%%%%%%%%%%%  									%%%%%%%%%%%%
%%%%%%%%%%%%  Final simulated data file produced %%%%%%%%%%%
%%%%%%%%%%%%  is 'ASIN-ChoiceModelData-FINAL.mat'%%%%%%%%%%%
%%%%%%%%%%%%  and is available on American      %%%%%%%%%%%% 
%%%%%%%%%%%%  Economic Review Data Website.     %%%%%%%%%%%%
%%%%%%%%%%%%  Simulated Data are intended to    %%%%%%%%%%%%
%%%%%%%%%%%%  be similar in spirit to actual    %%%%%%%%%%%%
%%%%%%%%%%%%  data used, which are both         %%%%%%%%%%%%
%%%%%%%%%%%%  proprietary and protected from    %%%%%%%%%%%%
%%%%%%%%%%%% public availability by HIPAA       %%%%%%%%%%%%
%%%%%%%%%%%% regulations and IRB protocols		%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% NOTE: The exercise here with simulated data is %%%%
%%%%%%%% meant to illustrate the methodology and is not %%%%
%%%%%%%% directly linked to the underlying data.		%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load 'ASIN-ChoiceModelData-FINAL.mat'

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

%%%%%%%%%%%%%%%%%%%%%%%%%% ESTIMATION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Now data are set up for simulated maximum likelihood estimation. To do estimation 			 		%%%%%%
%%%%%%%%%% we do constrained minimzation where the constraints are bounds										%%%%%%
%%%%%%%%%% on the parameters, wihch is especially important for those that cannot be							%%%%%%
%%%%%%%%%% less than 0 in theory. Otherwise, only accept estimates if parameters are 'far' from bounds 			%%%%%%
%%%%%%%%%% to ensure that estimation / nonlinear optimization is running correctly. 							%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% For constrained minimization we use KNITRO component 'ktrlink' which is similar in syntax 			%%%%%%
%%%%%%%%%% and purpose to the native MATLAB function fmincon. However, both use different underlying 			%%%%%%	
%%%%%%%%%% algorithms and ktrlink / KNITRO are generally more effective. If a user wants to translate the 		%%%%%%
%%%%%%%%%% code that follows to fmincon (since ktrlink is associated with the KNITRO package, which only		%%%%%%
%%%%%%%%%% some universities have) it is straightforward. Please contact me with any questions on this. 		%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% NOTE: sample code supplied here estimates choice model for one set of starting values. To get best   %%%%%%
%%%%%%%%%% estimates one shuold start from multiple starting values and choose the solution with the highest    %%%%%%
%%%%%%%%%% likelihood function value. This is the best way to ensure that you're getting a global optimum as    %%%%%%
%%%%%%%%%% well as a local optimum. See ZIENA / KNITRO 8 documention online for further discussion of this point.%%%%%
%%%%%%%%%%																										%%%%%%
%%%%%%%%%% There are multiple potential ways to employ multiple starting values:								%%%%%%
%%%%%%%%%%																										%%%%%%
%%%%%%%%%% (1) Just run the program here for many different initial parameters alpha0 (or set up longer 		%%%%%%
%%%%%%%%%%     program to do these runs). Here, this is what we did for maximum transparency. We started 		%%%%%%
%%%%%%%%%%	   at 5 different starting values, where four are commented out in the code below the alpha0		%%%%%%
%%%%%%%%%%	   actually used in this code. Note: the alpha0 supplied in this code had the best likelihood		%%%%%%
%%%%%%%%%%	   function value of the six candidate starting values. Note: for actual estimation like this 		%%%%%%
%%%%%%%%%%	   I would suggest more like 20-50 starting values to be confident that you're getting the best 	%%%%%%
%%%%%%%%%%	   parameters. 																						%%%%%%
%%%%%%%%%%																										%%%%%%
%%%%%%%%%% (2) KTRLINK / KNITRO has an automatic algorithm to start a program from multiple values, which can 	%%%%%%
%%%%%%%%%%     also be parallelized if desired. See pages 39-43 of ZIENA / KNITRO 8.0 manual for details. The 	%%%%%%
%%%%%%%%%%     program selects starting values randomly in between the lower and upper bounds, many options 	%%%%%%
%%%%%%%%%%	   such as number of starting values of sub regions they are drawn from. Here, we stick with option %%%%%%
%%%%%%%%%%	   (1) for transparency, since it works with fmincon, but this is a good option for ktrlink. 		%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: specific factors in estimation users can look to vary include:
% 1. Starting values (as mentioned above)
% 2. Smoothing factor used in likelihood function to keep it continuous (tests show that factor of power 6 used is sufficient, using higher numbers makes minimal difference)
% 3. Different non-linear optimizer (e.g. fmincon or some other optimizer instead of ktrlink) 
% 4. Different algorithm within given optimizer (e.g. active set vs. some other optimizer)
% 5. TolFun option: determines fineness of when optimizer stops, lower numbers imply longer time, but potentially more precise value
% 6. Bounds if any parameter estimates hitting bounds or you want to test likelihood function when certain paramaters are more constrained than is natural

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Set uppper and lower bounds for parameters to guide constrained optimization %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lb = [0,-1,-1,0;-Inf,0.1,-Inf,0.1; 1,1,1,1; -Inf,0,0,-5000;-5000,-5000,-5000,-5000;-5000,0,0,0;];
ub = [2,1,1,1; Inf,Inf,Inf,Inf; 20000,20000,20000,20000; Inf,20000,20000,5000;5000,5000,5000,5000;5000,0,0,0;];

%%%%%% Set default options for estimation. Set MaxFunEvals and MaxIter high enough so they only stop if something is wrong. 
%%%%%% Set TolFun = 10^-4, since function is on the order of 2,000 in value. Could set this finer, tradeoff is estimation is longer, similar solution.  

options = optimset('MaxFunEvals',40000,'MaxIter',40000,'TolFun',10^-4);

%%%%%% alpha0 is starting value for non-linear optimization. Must start within lb and ub. 

alpha0 = [0.06,0.003,0,0.04; -2500,700,-2200,800;300,800,300,800;-500,1250,1750,-500;0,0,0,0;0,0,0,0;]

% Starting values also run for this simulated data code, to find best likelihood value. Actual alpha0 used here above was one with highest likelihood value. 

%alpha0 = [-1500,0.04,0,0.03; 900,-0.00001,1500,-300;800,300,800,300;-2500,700,1500,-200;20,100,100,-100;-100,0,0,0;]
%alpha0 = [-3000,0.1,-0.0003,0.08; 600,0.00001,2500,-500;500,200,500,200;-1500,700,1200,-800;200,0,0,-200;-100,0,0,0;]
%alpha0 = [-2000,0.02,-0.0003,0.03; 1500,0,1700,-800;400,150,400,150;-1000,500,1500,-500;0,50,200,-100;100,0,0,0;]
%alpha0 = [-3200,0.05,0.003,0.06; 1000,0,1500,-700;150,150,150,150;-2000,800,1050,-700;-50,-50,100,100;-50,0,0,0;]

%%%%%%% Evaluate likelihood function at alpha0, helps to test whether or not function is working appropriately. 

Likelihood_ASIN_PSET_1(alpha0,R,EPS,nIs,nPlans,nTs,HTC,Sim,K,EUvector11,EUvector21,EUvector31,EUvector12,EUvector22,EUvector32,EUvector13,EUvector23,EUvector33,choice,choice12,choice22,choice32,choice13,choice23,choice33,Income,Tier2,IND,FSAY2,FSAY3,QS,managerX,CC1,CC2,CSAL2,CSAL3,MAge)

%%%%%%% Run estimation with ktrlink non-linear optimization routine. Again, fmincon, included in MATLAB, could be easily substituted with a few simple syntax modifications
%%%%%%% though the two non-linear optimizers are different. 

[alphasol,fval] = ktrlink( @(alpha) Likelihood_ASIN_PSET_1(alpha,R,EPS,nIs,nPlans,nTs,HTC,Sim,K,EUvector11,EUvector21,EUvector31,EUvector12,EUvector22,EUvector32,EUvector13,EUvector23,EUvector33,choice,choice12,choice22,choice32,choice13,choice23,choice33,Income,Tier2,IND,FSAY2,FSAY3,QS,managerX,CC1,CC2,CSAL2,CSAL3,MAge), alpha0,[],[],[],[],lb,ub,[],options);

%%%%%%% Report parameter estimates once estimation is done running. NOTE: for this code with 50 simulated values for the risk preference distribution 
%%%%%%% and 50 simulated health risk draws and approximately 2500 families, users should expect code to take anywhere between 2 to 6 hours to run 
%%%%%%% depending on the starting values supplied and true parameters. For the starting values / simulated data here approximately 4 hours.  
%%%%%%% Risk preference distribution estimates 

MeanIntercept =  alphasol(1,1)/100
IncomeSlope = alphasol(1,2)/100
AgeSlope = alphasol(1,3)/100
RiskPrefSTDEV = alphasol(1,4)/100

%%%%%%% PPO1200 / CDHP random coefficients

% Single 

CDHPmeanE = alphasol(2,1)
CDHPVarE = alphasol(2,2)

% Family 

CDHPmean = alphasol(2,3)
CDHPVar = alphasol(2,4)

%%%%% Epsilon variance estimates 

% Single 

EPSSTDEVPPO500E= alphasol(3,1)
EPSSTDDEVPPO1200E = alphasol(3,2)

% Family 

EPSSTDEVPPO500 = alphasol(3,3)
EPSSTDEVPPO1200 = alphasol(3,4)

%%%%%%% Preference for PPO500 / PPO1200 for high-cost families vs. those not high-cost

HTC = alphasol(4,1)

%%%%%%% Inertia / Switching Cost Estimates 

% Individual Intercept

SCE = alphasol(4,2)

% Family Intercept

SCF = alphasol(4,3)

% Slope for link between Inertia / Switching Costs and Demographics

FSA = alphasol(4,4)
IncomeTier = alphasol(5,1)
QuantSoph = alphasol(5,2)
manager = alphasol(5,3)
Chronic = alphasol(5,4)
SalientChange = alphasol(6,1) 

% Report value of likelihood function to compare with that from other starting values

fval 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% END OF ESTIMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% RESULTS %%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

%%%%%%% The results one should get from estimation, when running the code above on the simulated data are given 
%%%%%%% below. These are provided so that the user can verify that he / she is running the code properly.  

%%%%%%% Moreover, this is an important metric to use if you're running fmincon in place of ktrlink (or some other
%%%%%%% non-linear optimzation algorithm in lieu of ktrlink). In that event you can see how your results compare 
%%%%%%% to the ktrlink results (since this program is not always available). 

%%%%%%% After the estimates are presented below we present the 'true' parameters that the simulated data were 
%%%%%%% created from. The estimated parameters are quite close to the true values, though a few do deviate
%%%%%%% (such as the epsilon/RP standard deviations) which also have high standard errors in the actual estimation. 
%%%%%%% To the extent that the parameters differ from the underlying factors, it reflects the fact that 
%%%%%%% choices are being simulated for ~2300 people (similar to that actually used in estimation) rather than
%%%%%%% a larger simulated sample, under which the estimates shuold converge to the true values. Overall, the 
%%%%%%% simulated data estimates are quite close to the true values, reflecting the fact that the sample size 
%%%%%%% is large enough to accurately capture the parameters modeled. 

% NOTE: Compared to 'true' values used to simulate data, switching cost estimates are exactly right, risk prefs and CDHP coefficient are
% 'nearby' and epsilon are a little off. The CDHP and risk prefs are slightly off in 'opposite' directions for their impact on choices
% suggesting that the goodness of fit for the model will be right one but that these parameters can substitute for one another with 
% finite data to some extent. 

%%%%%%%%%%%%%%% Estimates from this code on simulated data with ktrlink and primary starting value given above
%%%%%%%%%%%%%%% Run by Ben Handel in May 2013. Variable names correspond to  those given above for output of 
%%%%%%%%%%%%%%% code estimates, and should be self explanatory upon reading the paper / this and prior code.   

%%%%%%% PPO1200 / CDHP random coefficients

% Family 

%CDHPmean = -2,698
%CDHPVar = 675

% Single 

%CDHPmeanE = -2,625
%CDHPVarE = 613

%%%%%%% Risk preference distribution estimates 

%MeanIntercept = 0.000519  
%IncomeSlope = 0.000011
%AgeSlope = 0.0000013
%RiskPrefSTDEV = 0.000024

%%%%%%% Preference for PPO500 / PPO1200 for high-cost families vs. those not high-cost

%HTCS = -792

%%%%%%% Inertia / Switching Cost Estimates 

% Family Intercept 

%SCest = 1741

% Individual Intercept

%SCE = 1346

% Slope for link between Inertia / Switching Costs and Demographics

%FSA = -782
%IncomeTier = -32 
%QuantSoph = 77
%manager = 227
%Chronic = 78
%SalientChange = 219  

%%%%% Epsilon variance estimates 

% Family 

%EPSSTDEVPPO1200 = 510
%EPSSTDEVPPO500 = 76

% Single 

%EPSSTDDEVPPO1200E = 519
%EPSSTDEVPPO500E= 72

%* This estimate hit the natural lower bound, implying that other heterogeneity in the estiamtes is explaining 
% heterogeneity in choices better than this idiosyncratic shock. Note, these epsilons for PPO500 are the  parameters  
% which do not come with estimated value right in neighborhood of actual value, given below. When you replace the low 
% estimated value with the true value, the likelihood function value is worse, so the estimates are correct given 
% the simulated data, this element just has very low predictive power. 

%%%%% LIKELIHOOD FUNCTION VALUE

% Output of program is 1,854.4 which is negative of true likelihood function value, which is a negative number = -1,854.4
% This maximizes the log likelihood function value, which will always be negative since the output for each person is a probability < 1
% so the log terms will always be negative. This was the best likelihood value of the 5 starting values given above. As noted
% above, in actual estimation you either want to use the multi-start feature of ktrlink or manually use more starting values
% to make sure you're getting global optimum. In this case with the simulated data, the results for the different starting values
% are all in the same neighborhood of one another, despite the starting values being quite different.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% 'True' parameter values that were used to simulate the choices contained in the simulated dataset	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
%%%%%%% PPO1200 / CDHP random coefficients

% Family 

%CDHPmean = -2200
%CDHPVar = 800

% Single 

%CDHPmeanE = -2400
%CDHPVarE = 900

%%%%%%% Risk preference distribution estimates 

%MeanIntercept =  0.00075
%IncomeSlope = 0.00004
%AgeSlope = 0.0000025
%RiskPrefSTDEV = 0.0005

%%%%%%% Preference for PPO500 / PPO1200 for high-cost families vs. those not high-cost

%HTCS = -800

%%%%%%% Inertia / Switching Cost Estimates 

% Family Intercept

%SCest = 1750

% Individual Intercept

%SCE = 1250

% Slope for link between Inertia / Switching Costs and Demographics

%FSA = -600
%IncomeTier = -50
%QuantSoph = 10
%manager = 200
%Chronic = 80
%SalientChange = 150 

%%%%% Epsilon variance estimates 

% Family 

%EPSSTDEVPPO1200 = 800
%EPSSTDEVPPO500 = 400

% Single 

%EPSSTDDEVPPO1200E = 600
%EPSSTDEVPPO500E= 300
						
			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%		END OF FILE					%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%			
			
			
			
			
			