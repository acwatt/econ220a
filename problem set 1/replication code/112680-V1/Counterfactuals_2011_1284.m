%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   Sample Code for Counterfactual          %%%%%
%%%%%%%%%%%%   Analysis using Simulated Data           %%%%%
%%%%%%%%%%%%   Code presents analysis of specific      %%%%%
%%%%%%%%%%%%   counterfactual: 75% reduction in inertia %%%%
%%%%%%%%%%%%   and welfare analysis that presumes that 	%%%%
%%%%%%%%%%%%   inertia is not welfare relevant in and of %%%
%%%%%%%%%%%%   itself. In actual analysis / paper the  %%%%%
%%%%%%%%%%%%   author examines how the results vary along %%
%%%%%%%%%%%%   both of these dimenstions (e.g. inertia %%%%%
%%%%%%%%%%%%   reduction varies between 0 and 100% and %%%%%
%%%%%%%%%%%%   inertia ranges from 0% to 100% welfare  %%%%%
%%%%%%%%%%%%   relevant in and of itself. It is        %%%%%
%%%%%%%%%%%%   straightforward to extend the code here %%%%%
%%%%%%%%%%%%   to examine these other coutnerfactuals. %%%%%
%%%%%%%%%%%%   Note: analysis here studies counterfactual %%
%%%%%%%%%%%%   for 5 time periods from initial active  %%%%%
%%%%%%%%%%%%   choice year instead of 7 as done in paper. %%
%%%%%%%%%%%%   Analysis can easily be extended to 7 by    %%
%%%%%%%%%%%%   extending code here in obvious ways for   %%%
%%%%%%%%%%%%   additional years. 					   %%%%%
%%%%%%%%%%%%									%%%%%%%%%%%%
%%%%%%%%%%%%									%%%%%%%%%%%%
%%%%%%%%%%%%  Final simulated data file used here      %%%%%
%%%%%%%%%%%%  is 'ASIN-ChoiceModelData-FINAL.mat       %%%%%
%%%%%%%%%%%%  and is available on American             %%%%% 
%%%%%%%%%%%%  Economic Review Data Website.            %%%%%
%%%%%%%%%%%%  These are the same data used as an input %%%%%
%%%%%%%%%%%%  into choice model estimation, as done in code%
%%%%%%%%%%%%  'EsimationCode_2011_1284.m'			   %%%%%
%%%%%%%%%%%%  See that file and the README file for    %%%%%
%%%%%%%%%%%%  additional details on simulated data and %%%%%
%%%%%%%%%%%%  overall code available for this paper on %%%%%
%%%%%%%%%%%%  the American Economic Review website %%%%%%%%%
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

% NOTE: File gives example for main counterfactual with endogenous re-pricing. Coutnerfactual 
% with no re-pricing easy to derive from this code, seen in welfare calculations at end of file. 

load 'ASIN-ChoiceModelData-FINAL_2011_1284.mat'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% See corresponding data desription document for detailed description  %%%%%%%%%%
%%%%%%%%%% of variables in simulated data. 										%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. First part of code sets up all preferences apart from inertia in order to have baseline for all preferences when inertia is reduced
%% 2. Next, we do 'accept-reject' on random cofficients. This restricts unobserved heteroegeneity of family to be consistent with 
%%  observed choices made in simulated data, since this is correct set of unobserved heterogeneity draws to use for each family in counterfactual simulations
%% 3. Part 3 reduces inertia and simulates market outcomes and welfare over time. Done for 75% reduction in inertia and no direct welfare impact of inertia
%% as discussed above and in paper. All counterfactuals from the paper can be built with small variations on this code 
%% 4. NOTE: the counterfactual code that simulates new prices / choices with redueced inertia by 75% is example given here. To run welfare analysis user must repeat this
%%   code twice (once for baseline without reduced inertia, once for counterfactual code given with reduced inertia) and save output of each run. 
%%   Then, both these files are used as inputs to welfare code at end to compute change from baseline welfare due to reduced inertia. End of file gives instructions 
%%   before welfare example code given: to run baseline change SCRED to 1.0 from 0.25 (line 849 below)
%% 5. Code supplied so users can use code with new / different data and learn with simulated data. Results not meant to mimic paper, simulated data much simpler on 
%%    variety of dimensions and underlying foundations / parameters quite different than in actual data. However, code illustrates process used in actual analysis.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% Replace number of simulations for risk preference / epsilon parameters with larger number 
%%%%% 400 instead of 50 because this slows down estimation but counterfactual analysis 
%%%%% is fast regardless so use larger number of draws here. 

Sim = 400; 

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Set income equal to year 2 income. Done since random coefficient for risk preferences, 					%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% which depends on income, is constant over time in panel within person by assumption.  					%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Income = Inc2;

%%%%% Set age equal to maximum age in family

MAge = max(Ages,[],2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Extend 'observed' price matrices from simulated data to have larger Sim dimension %%%%%%%
%%%% used in this counterfactual analysis. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PPO250 Plan Prices over 3 'observed' years 

P1O = zeros(3,nIs,Sim);
P10(:,:,1:50)=P1;
P10(:,:,51:100)=P1;
P10(:,:,101:150)=P1;
P10(:,:,151:200)=P1;
P10(:,:,201:250)=P1;
P10(:,:,251:300)=P1;
P10(:,:,301:350)=P1;
P10(:,:,351:400)=P1;

%%%%%%%% Set up P1 to be counterfactual price matrix for 5 year time horizon

P1 = zeros(5,nIs,Sim);
P1(1:3,:,:) = P10;

%% PPO500 Plan Prices over 3 'observed' years 

P2O = zeros(3,nIs,Sim);
P20(:,:,1:50)=P2;
P20(:,:,51:100)=P2;
P20(:,:,101:150)=P2;
P20(:,:,151:200)=P2;
P20(:,:,201:250)=P2;
P20(:,:,251:300)=P2;
P20(:,:,301:350)=P2;
P20(:,:,351:400)=P2;

%%%%%%%% Set up P1 to be counterfactual price matrix for 5 year time horizon

P2 = zeros(5,nIs,Sim);
P2(1:3,:,:) = P20;

%% PPO1200 Plan Prices over 3 'observed' years 

P3O = zeros(3,nIs,Sim);
P30(:,:,1:50)=P3;
P30(:,:,51:100)=P3;
P30(:,:,101:150)=P3;
P30(:,:,151:200)=P3;
P30(:,:,201:250)=P3;
P30(:,:,251:300)=P3;
P30(:,:,301:350)=P3;
P30(:,:,351:400)=P3;

%%%%%%%% Set up P1 to be counterfactual price matrix for 5 year time horizon

P3 = zeros(5,nIs,Sim);
P3(1:3,:,:) = P30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulate Coefficients from Model Output

%%%%%%%%%%%%%%%%%%%%%%%%%%%% PPO1200 Random Coefficients

mu1200S = -2647;
sigma1200S = 601;
S1200 = normrnd(mu1200S, sigma1200S, Sim,1); 

mu1200F = -2731;
sigma1200F = 693;
F1200 = normrnd(mu1200F, sigma1200F, Sim,1); 

CDHPt = zeros(nIs,Sim);

for i = 1:nIs
    if IND(1,i)==0
		CDHPt(i,:) = (F1200)';
    end	
    if IND(1,i)==1
		CDHPt(i,:) = (S1200)';
    end    
end

CDHP = zeros(nIs,K,Sim);

for k =1:K
	CDHP(:,k,:) = CDHPt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Plan-Time Specific Errors %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create Shocks for Each Plan and Time Period relative to PPO250
% SD assumed similar to those estimated in paper

muE = [0,0,0,0,0,0];
sigmaE = eye(6);
EPS = mvnrnd(muE, sigmaE, Sim); 

ES500 = 8;
ES1200 = 512;
EF500 = 1;
EF1200 = 502;

Epsilon1200t1 = zeros(nIs,Sim);
Epsilon500t1 = zeros(nIs,Sim);
Epsilon1200t2 = zeros(nIs,Sim);
Epsilon500t2 = zeros(nIs,Sim);
Epsilon1200t3 = zeros(nIs,Sim);
Epsilon500t3 = zeros(nIs,Sim);
Epsilon1200t4 = zeros(nIs,Sim);
Epsilon500t4 = zeros(nIs,Sim);
Epsilon1200t5 = zeros(nIs,Sim);
Epsilon500t5 = zeros(nIs,Sim);

for i = 1:nIs
    if IND(1,i)==0
		Epsilon1200t1(i,:) = (EPS(:,1)*EF1200)';
		Epsilon500t1(i,:) = (EPS(:,2)*EF500)';
		Epsilon1200t2(i,:) = (EPS(:,3)*EF1200)';
		Epsilon500t2(i,:) = (EPS(:,4)*EF500)';
		Epsilon1200t3(i,:) = (EPS(:,5)*EF1200)';
		Epsilon500t3(i,:) = (EPS(:,6)*EF500)';
    end
    if IND(1,i)==1
		Epsilon1200t1(i,:) = (EPS(:,1)*ES1200)';
		Epsilon500t1(i,:) = (EPS(:,2)*ES500)';
		Epsilon1200t2(i,:) = (EPS(:,3)*ES1200)';
		Epsilon500t2(i,:) = (EPS(:,4)*ES500)';
		Epsilon1200t3(i,:) = (EPS(:,5)*ES1200)';
		Epsilon500t3(i,:) = (EPS(:,6)*ES500)';
    end    
end

muE = [0,0,0,0];
sigmaE = eye(4); 
EPS = mvnrnd(muE, sigmaE, Sim); 

for i = 1:nIs
    if IND(1,i)==0
		Epsilon1200t4(i,:) = (EPS(:,1)*EF1200)';
		Epsilon500t4(i,:) = (EPS(:,2)*EF500)';
		Epsilon1200t5(i,:) = (EPS(:,3)*EF1200)';
		Epsilon500t5(i,:) = (EPS(:,4)*EF500)';
    end
    if IND(1,i)==1
		Epsilon1200t4(i,:) = (EPS(:,1)*ES1200)';
		Epsilon500t4(i,:) = (EPS(:,2)*ES500)';
		Epsilon1200t5(i,:) = (EPS(:,3)*ES1200)';
		Epsilon500t5(i,:) = (EPS(:,4)*ES500)';
	end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Generate epsilon matrix orthogonal with orthogonal health shocks included   %%%%%%%%%%%
%%%%%%%%% to make computations faster later in code and include fewer loops 		  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Epsilon12001 = zeros(nIs,K,Sim);
Epsilon5001 = zeros(nIs,K,Sim);
Epsilon12002 = zeros(nIs,K,Sim);
Epsilon5002 = zeros(nIs,K,Sim);
Epsilon12003 = zeros(nIs,K,Sim);
Epsilon5003 = zeros(nIs,K,Sim);
Epsilon12004 = zeros(nIs,K,Sim);
Epsilon5004 = zeros(nIs,K,Sim);
Epsilon12005 = zeros(nIs,K,Sim);
Epsilon5005 = zeros(nIs,K,Sim);

for k =1:K
	Epsilon12001(:,k,:) = Epsilon1200t1;
	Epsilon5001(:,k,:) = Epsilon500t1;
	Epsilon12002(:,k,:) = Epsilon1200t2;
	Epsilon5002(:,k,:) = Epsilon500t2;
	Epsilon12003(:,k,:) = Epsilon1200t3;
	Epsilon5003(:,k,:) = Epsilon500t3;
	Epsilon12004(:,k,:) = Epsilon1200t4;
	Epsilon5004(:,k,:) = Epsilon500t4;
	Epsilon12005(:,k,:) = Epsilon1200t5;
	Epsilon5005(:,k,:) = Epsilon500t5;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Switching Cost Simulations with Actual Switching Cost Estimates %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Not Reduced Here, reduced for actual counterfactual simulation later in code %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SC2 = ones(nIs,Sim)*1786;
SC3 = ones(nIs,Sim)*1786;
SCE = 1376*ones(1,Sim);

for i = 1:nIs
    if IND(1,i)==1
        SC2(i,:) = SCE;
        SC3(i,:) = SCE;
    end   
end

SC2 = SC2 + (FSAY2*ones(1,Sim)*(-749)) + (Income*(-40)*ones(1,Sim))+ (QS*ones(1,Sim)*35) + (managerX*ones(1,Sim)*192) + (CC1*ones(1,Sim)*41) + (CSAL2*ones(1,Sim)*179);
SC3 = SC3 + (FSAY2*ones(1,Sim)*(-749)) + (Income*(-40)*ones(1,Sim))+ (QS*ones(1,Sim)*35) + (managerX*ones(1,Sim)*192) + (CC2*ones(1,Sim)*41) + (CSAL3*ones(1,Sim)*179);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Risk Preferences %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Risk preference mean intercept %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RAint = .000513;

%%%%%%%%%%% Risk preference mean slope for income %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RAinc = .000012;

%%%%%%%%%%% Risk preference mean slope for age %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RAage = .0000013; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RAmean = RAint*ones(1,nIs) + RAinc*(Income') + RAage*(MAge');
RAmean = (RAmean'*ones(1,Sim))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create Standard Dev. of Risk Preference Random Coefficient %%%%%%%%%%%%%%

muRA = 0;
sigmaRA = .00025;
R = normrnd(muRA, sigmaRA, Sim,1);
RAvar = R*ones(1,nIs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create matrix with Sim CARA coefficient simulations for nIs families 

RAi = (RAmean + RAvar)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create matrix where we have K health shock draws as orthogonal dimension %%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% to individual risk preference estimations. This makes estiamtion / computation %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% faster later in the code.  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RA = zeros(nIs,K,Sim);

for k = 1:K
RA(:,k,:) = RAi;
end

RA = max(RA,0.00000001);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Include Preference by HTC Families %%%%%%%%
%%%%%%%%%%%%%%% added as Negative to Other Plans %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CHTC = -773*HTC; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Generate expected utility baseline monetary values 
%%%%%%%%% for each state plan and year. Note that years after 'observed' 
%%%%%%%%% years, years 4-7, take on same health risk distributions as 
%%%%%%%%% year 3 to hold this factor constant in counterfactual simulations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

EUvector11 = zeros(nIs,K,Sim);
EUvector21 = zeros(nIs,K,Sim);
EUvector31 = zeros(nIs,K,Sim);
EUvector12 = zeros(nIs,K,Sim);
EUvector22 = zeros(nIs,K,Sim);
EUvector32 = zeros(nIs,K,Sim);
EUvector13 = zeros(nIs,K,Sim);
EUvector23 = zeros(nIs,K,Sim);
EUvector33 = zeros(nIs,K,Sim);
EUvector14 = zeros(nIs,K,Sim);
EUvector24 = zeros(nIs,K,Sim);
EUvector34 = zeros(nIs,K,Sim);
EUvector15 = zeros(nIs,K,Sim);
EUvector25 = zeros(nIs,K,Sim);
EUvector35 = zeros(nIs,K,Sim);

EUvector11x = zeros(nIs,K,Sim);
EUvector21x = zeros(nIs,K,Sim);
EUvector31x = zeros(nIs,K,Sim);
EUvector12x = zeros(nIs,K,Sim);
EUvector22x = zeros(nIs,K,Sim);
EUvector32x = zeros(nIs,K,Sim);
EUvector13x = zeros(nIs,K,Sim);
EUvector23x = zeros(nIs,K,Sim);
EUvector33x = zeros(nIs,K,Sim);

W = 75000;

SwC12 = zeros(nIs,K,Sim);
SwC22 = zeros(nIs,K,Sim);
SwC32 = zeros(nIs,K,Sim);
SwC13 = zeros(nIs,K,Sim);
SwC23 = zeros(nIs,K,Sim);
SwC33 = zeros(nIs,K,Sim);
SwC14 = zeros(nIs,K,Sim);
SwC24 = zeros(nIs,K,Sim);
SwC34 = zeros(nIs,K,Sim);
SwC15 = zeros(nIs,K,Sim);
SwC25 = zeros(nIs,K,Sim);
SwC35 = zeros(nIs,K,Sim);

choice12 = zeros(nIs,Sim);
choice22 = zeros(nIs,Sim);
choice32 = zeros(nIs,Sim);
choice13 = zeros(nIs,Sim);
choice23 = zeros(nIs,Sim);
choice33 = zeros(nIs,Sim);
choice14 = zeros(nIs,Sim);
choice24 = zeros(nIs,Sim);
choice34 = zeros(nIs,Sim);
choice15 = zeros(nIs,Sim);
choice25 = zeros(nIs,Sim);
choice35 = zeros(nIs,Sim);

Tol = 0.999;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PART II: Simulate choices for all draws of unobserved heterogeneity and keep the values of 
%%%%%%%%%% unobserved heterogeneity for each person that are consistent with the actual / 'observed' 
%%%%%%%%%% choices that they make in the simulated data. Accept-reject draws for year 3. These are 
%%%%%%%%%% draws used in coutnerfactual simulations with reduced inertia and are extended to years 4-7
%%%%%%%%%% as consistent in years 1-3. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

EUvector11 = exp(-((EUvector11).*max(RA,0.0000001)));
EUvector21 = exp(-((EUvector21+CHTC+Epsilon5001).*max(RA,0.0000001)));
EUvector31 = exp(-((EUvector31 + CDHP +CHTC + Epsilon12001).*max(RA, 0.0000001)));

EUvectorZ = zeros(nIs,nPlans,Sim);

for n = 1:Sim
    EUvectorZ(:,1,n) = (1/K)*sum(reshape(EUvector11(:,:,n),nIs,K),2);
    EUvectorZ(:,2,n) = (1/K)*sum(reshape(EUvector21(:,:,n),nIs,K),2);
    EUvectorZ(:,3,n) = (1/K)*sum(reshape(EUvector31(:,:,n),nIs,K),2);
end

%%%%%%%%%%%% Now create mirror EUvector / Expected Utility Matrices to Use in Accept-Reject Calculation

for i = 1:nIs
    for k = 1:K
        EUvector11x(i,k,:) = max(0.1,(W - reshape(PPO250OP1(i,k,:),1,Sim) - reshape(P1(1,i,:),1,Sim)));
        EUvector21x(i,k,:) = max(0.1,(W - reshape(PPO500OP1(i,k,:),1,Sim) - reshape(P2(1,i,:),1,Sim)));
        EUvector31x(i,k,:) = max(0.1,(W - reshape(PPO1200OP1(i,k,:),1,Sim) - reshape(P3(1,i,:),1,Sim)));
        EUvector12x(i,k,:) = max(0.1,(W - reshape(PPO250OP2(i,k,:),1,Sim) - reshape(P1(2,i,:),1,Sim)));
        EUvector22x(i,k,:) = max(0.1,(W - reshape(PPO500OP2(i,k,:),1,Sim) - reshape(P2(2,i,:),1,Sim)));
        EUvector32x(i,k,:) = max(0.1,(W - reshape(PPO1200OP2(i,k,:),1,Sim) - reshape(P3(2,i,:),1,Sim)));
        EUvector13x(i,k,:) = max(0.1,(W - reshape(PPO250OP3(i,k,:),1,Sim) - reshape(P1(3,i,:),1,Sim)));
        EUvector23x(i,k,:) = max(0.1,(W - reshape(PPO500OP3(i,k,:),1,Sim) - reshape(P2(3,i,:),1,Sim)));
        EUvector33x(i,k,:) = max(0.1,(W - reshape(PPO1200OP3(i,k,:),1,Sim) - reshape(P3(3,i,:),1,Sim)));
    end
end

EUvector11x = exp(-((EUvector11x).*max(RA,0.0000001)));
EUvector21x = exp(-((EUvector21x+CHTC+Epsilon5001).*max(RA,0.0000001)));
EUvector31x = exp(-((EUvector31x + CDHP +CHTC + Epsilon12001).*max(RA, 0.0000001)));

EUvectorZx = zeros(nIs,nPlans-1,Sim);

for n = 1:Sim
    EUvectorZx(:,1,n) = (1/K)*sum(reshape(EUvector11x(:,:,n),nIs,K),2);
    EUvectorZx(:,2,n) = (1/K)*sum(reshape(EUvector21x(:,:,n),nIs,K),2);
    EUvectorZx(:,3,n) = (1/K)*sum(reshape(EUvector31x(:,:,n),nIs,K),2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Accept Reject on the Random Coefficients %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% for year 1 'observed chocies' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AcceptReject = zeros(nIs,Sim);

for i = 1:nIs
   for s = 1:Sim
     
%%%%%%%%%% Accept random coefficients draw if expected utility it greater for actual plan chosen. Remeber, actual utility is negative of EUvector 
%%%%%%%%%% which is why we want this quantity to be lower in value for the plan actually chosen. Same in likelihood function / estimation. 
	 
	   if choice(i,1) == 1
           if (EUvectorZx(i,1,s)*Tol <= EUvectorZ(i,2,s)) & (EUvectorZx(i,1,s)*Tol  <= EUvectorZ(i,3,s))
           AcceptReject(i,s) =1;     
           end    
        end

        if choice(i,1) == 2
           if (EUvectorZx(i,2,s)*Tol  <= EUvectorZ(i,1,s)) & (EUvectorZx(i,2,s)*Tol  <= EUvectorZ(i,3,s))
           AcceptReject(i,s) =1 ;    
           end    
        end
        
        if choice(i,1) == 3
           if (EUvectorZx(i,3,s)*Tol  <= EUvectorZ(i,1,s)) & (EUvectorZx(i,3,s)*Tol  <= EUvectorZ(i,2,s))
           AcceptReject(i,s) =1  ;   
           end    
        end 
		
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Baseline AcceptReject for year 2 given actual prices %%%%%%%%%

%%%%%%%%%%%%% Bring in inertia estimates here as well since these apply to observed choices in years 2 and3

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
    end
end

for k = 1:K
    SwC12(:,k,:) = SC2.*choice12;
    SwC22(:,k,:) = SC2.*choice22;
    SwC32(:,k,:) = SC2.*choice32;
end

EUvector12 = exp(-((EUvector12 +SwC12).*max(RA,0.0000001)));
EUvector22 = exp(-((EUvector22 + SwC22 + CHTC + Epsilon5002).*max(RA,0.0000001)));
EUvector32 = exp(-((EUvector32 + CDHP + CHTC + SwC32 + Epsilon12002).*max(RA, 0.0000001)));

%%%%%%%%%%%%%%%%%%% Now sum CARA VNM to get expected utility %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Set up matrices for all forthcoming years, including %%%%%%
%%%%%%%%%%%%%%%%%%% those we will do later in counterfactuals %%%%%%%%%%%%%%%%%

EUvectorA = zeros(nIs,nPlans,Sim);
EUvectorAx = zeros(nIs,nPlans,Sim);
EUvectorB = zeros(nIs,nPlans,Sim);
EUvectorBx = zeros(nIs,nPlans,Sim);
EUvectorC = zeros(nIs,nPlans,Sim);
EUvectorD = zeros(nIs,nPlans,Sim);
EUvectorE = zeros(nIs,nPlans,Sim);
EUvectorF = zeros(nIs,nPlans,Sim);
EUvectorDx = zeros(nIs,nPlans,Sim);
EUvectorEx = zeros(nIs,nPlans,Sim);
EUvectorFx = zeros(nIs,nPlans,Sim);

for n = 1:Sim
    EUvectorA(:,1,n) = (1/K)*sum(reshape(EUvector12(:,:,n),nIs,K),2);
    EUvectorA(:,2,n) = (1/K)*sum(reshape(EUvector22(:,:,n),nIs,K),2);
    EUvectorA(:,3,n) = (1/K)*sum(reshape(EUvector32(:,:,n),nIs,K),2);
end

EUvector12x = exp(-((EUvector12x +SwC12).*max(RA,0.0000001)));
EUvector22x = exp(-((EUvector22x + SwC22 + CHTC + Epsilon5002).*max(RA,0.0000001)));
EUvector32x = exp(-((EUvector32x + CDHP + CHTC + SwC32 + Epsilon12002).*max(RA, 0.0000001)));

for n = 1:Sim
    EUvectorAx(:,1,n) = (1/K)*sum(reshape(EUvector12x(:,:,n),nIs,K),2);
    EUvectorAx(:,2,n) = (1/K)*sum(reshape(EUvector22x(:,:,n),nIs,K),2);
    EUvectorAx(:,3,n) = (1/K)*sum(reshape(EUvector32x(:,:,n),nIs,K),2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Accept reject on random coefficients for year 2 %%%%%%%%%%%%%%

AcceptReject2 = zeros(nIs,Sim);

for i = 1:nIs
   for s = 1:Sim
        
       if choice(i,2) == 1
           if (EUvectorAx(i,1,s)*Tol  <= EUvectorA(i,2,s)) & (EUvectorAx(i,1,s)*Tol  <= EUvectorA(i,3,s))
           AcceptReject2(i,s) =1 ;
           end    
        end

        if choice(i,2) == 2
           if (EUvectorAx(i,2,s)*Tol  <= EUvectorA(i,1,s)) & (EUvectorAx(i,2,s)*Tol  <= EUvectorA(i,3,s))
           AcceptReject2(i,s) =1 ; 
           end    
        end
        
        if choice(i,2) == 3
           if EUvectorAx(i,3,s)*Tol  <= EUvectorA(i,1,s)*Tol  & (EUvectorAx(i,3,s)*Tol  <= EUvectorA(i,2,s))
           AcceptReject2(i,s) =1 ;  
           end    
        end
    
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Baseline for year 3 AcceptReject given actual prices %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nIs
    for s = 1:Sim
       if choice(i,2)==1
        choice13(i,s) = 1;
       end
       if choice(i,2)==2
        choice23(i,s) = 1;
       end
       if choice(i,2)==3
        choice33(i,s) = 1;
       end
    end
end

for k = 1:K
    SwC13(:,k,:) = SC3.*choice13;
    SwC23(:,k,:) = SC3.*choice23;
    SwC33(:,k,:) = SC3.*choice33;
end

EUvector13 = exp(-((EUvector13 +SwC13).*max(RA,0.0000001)));
EUvector23 = exp(-((EUvector23 + SwC23 + CHTC + Epsilon5003).*max(RA,0.0000001)));
EUvector33 = exp(-((EUvector33 + CDHP + CHTC + SwC33 + Epsilon12003).*max(RA, 0.0000001)));

for n = 1:Sim
    EUvectorB(:,1,n) = (1/K)*sum(reshape(EUvector13(:,:,n),nIs,K),2);
    EUvectorB(:,2,n) = (1/K)*sum(reshape(EUvector23(:,:,n),nIs,K),2);
    EUvectorB(:,3,n) = (1/K)*sum(reshape(EUvector33(:,:,n),nIs,K),2);
end

EUvector13x = exp(-((EUvector13x +SwC13).*max(RA,0.0000001)));
EUvector23x = exp(-((EUvector23x + SwC23 + CHTC + Epsilon5003).*max(RA,0.0000001)));
EUvector33x = exp(-((EUvector33x + CDHP + CHTC + SwC33 + Epsilon12003).*max(RA, 0.0000001)));

for n = 1:Sim
    EUvectorBx(:,1,n) = (1/K)*sum(reshape(EUvector13x(:,:,n),nIs,K),2);
    EUvectorBx(:,2,n) = (1/K)*sum(reshape(EUvector23x(:,:,n),nIs,K),2);
    EUvectorBx(:,3,n) = (1/K)*sum(reshape(EUvector33x(:,:,n),nIs,K),2);
end

AcceptReject3 = zeros(nIs,Sim);

for i = 1:nIs
   for s = 1:Sim
        
       if choice(i,3) == 1
           if EUvectorBx(i,1,s)*Tol  <= EUvectorB(i,2,s) & EUvectorBx(i,1,s)*Tol  <= EUvectorB(i,3,s)
           AcceptReject3(i,s) =1  ;   
           end    
        end

        if choice(i,3) == 2
           if EUvectorBx(i,2,s)*Tol  <= EUvectorB(i,1,s) & EUvectorBx(i,2,s)*Tol  <= EUvectorB(i,3,s)
           AcceptReject3(i,s) =1 ;    
           end    
        end
        
        if choice(i,3) == 3
           if EUvectorBx(i,3,s)*Tol  <= EUvectorB(i,1,s) & EUvectorBx(i,3,s)*Tol  <= EUvectorB(i,2,s)
           AcceptReject3(i,s) =1;     
           end    
        end
    
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Accept-Reject FINAL 
%%%%%%%%%%%%%%%% Aggregate all three years 
%%%%%%%%%%%%%%%% Proceeds in several steps. First steps are strictest. If few/no draws 
%%%%%%%%%%%%%%%% satisfy those steps, use looser criteria so all are included in simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AcceptRejectF = zeros(nIs,Sim);

%%%%% Step 1: If all accept rejects in panel =1, then aggregate = 1

for i = 1:nIs
    for s = 1:Sim
    if AcceptReject(i,s) ==1 & AcceptReject2(i,s)==1 & AcceptReject3(i,s) == 1
       AcceptRejectF(i,s) = 1;
    end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Check some specific tests to see if people have any accept-reject draws = 1
%%%% Switchers are checked especially since switching may leads to more rejections
%%%% Indicator for whether or not you ever switch

Switch = zeros(nIs,1);

for i = 1:nIs
   if (choice(i,1)~=choice(i,2) | choice(i,2)~=choice(i,3))
       Switch(i,1)=1;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Check to see how many have positive final accept rejects vs. not. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OK = zeros(nIs,1);

for i = 1:nIs
   if mean(AcceptRejectF(i,:),2) == 0 
       OK(i,1)=1;
   end
end

%%% Num people with no positive accepts

mean(OK)

%%%%% Check to see how many people have 0 mean accept reject and have switched at some point in process. 
 
OK1 = zeros(nIs,1);

for i = 1:nIs
   if OK(i,1) ==1 & Switch(i,1)==1
      OK1(i,1)=1; 
   end
end

%%%%% For people with no accepts for criterion one above, set Accept-Reject for random coefficients = 1 if coefficients are accepted for at least one choice made. 

for i = 1:nIs
   for s = 1:Sim
        if (OK(i,1) == 1) & (AcceptReject(i,s)==1 | AcceptReject2(i,s)==1 | AcceptReject3(i,s)==1)
         AcceptRejectF(i,s)=1; 
         end
   end
end

%%%%% Recreate indicator. Now, if person still has 0 accepts, set all equal to accept, essentially saying we have no information about 
%%%%% this person above and beyond distributions estimated for population. 

OKN = zeros(nIs,1);

for i = 1:nIs
   if mean(AcceptRejectF(i,:),2) == 0 
       OKN(i,1)=1;
   end
end

for i = 1:nIs
    for s = 1:Sim
         if (OKN(i,1) == 1) 
         AcceptRejectF(i,s)=1; 
         end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  PART II calculating accepts-reject is complete. Now move ahead 
%%%%%%%%%% to part III which does actual counterfactual simulation for 
%%%%%%%%%% reduction in inertia. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% This is PART III, the actual counterfactual simulation to reduce inertia. As discussed above and in README file, here we do this %%%%%%%%%%%%%%%
%%%%% with one example where inertia is reduced by 75% and does inertia itself is not welfare relevant. Paper analyses many other options %%%%%%%%%%%%
%%%%% which can easily be looked at by modifying this code. 																				%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% Generate coutnerfactual inertia to use in counterfactual along with regular preferences simulated / given above directly from estimation. 

SCC2 = zeros(nIs,Sim);
SCC3 = zeros(nIs,Sim);

%%%%%% Inertia is reduced to this proportion of estimated inertia in counterfactual 
%%%%%% Thus in this example code inertia is reduced by 75% 

SCRED = 0.25;

%%%%%%
 SCC2 = SC2*SCRED;
 SCC3 = SC3*SCRED;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Simulate year 0 choices and expected plan costs
%%%%%%%%%%%%%%%%%%%%%%%% as a function of those choices, which will be passed through
%%%%%%%%%%%%%%%%%%%%%%%% to new prices for the year after. See discussion 
%%%%%%%%%%%%%%%%%%%%%%%% in section 6 of actual paper for more detail on assumed pricing model
%%%%%%%%%%%%%%%%%%%%%%%% exactly the same as that used within the firm itself
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CC0 = zeros(nIs,Sim);
MktShare0 = zeros(nPlans,1);

%%%%%%%%% Since plans are priced as a function of family status
%%%%%%%%% generate counterfactual enrollment and costs as a function 
%%%%%%%%% of family status

%%%%% Single cost and number of people enrolling

TOTALCOSTEZ = zeros(3,1);
TOTALNUMEZ = zeros(3,1);
TOTALNUMEZi = zeros(3,nIs);

%%%%% With spouse cost and number of people enrolling

TOTALCOSTSZ = zeros(3,1);
TOTALNUMSZ = zeros(3,1);
TOTALNUMSZi = zeros(3,nIs);

%%%%% With child(ren) cost and number of people enrolling

TOTALCOSTCZ = zeros(3,1);
TOTALNUMCZ = zeros(3,1);
TOTALNUMCZi = zeros(3,nIs);

%%%%% With spouse+children cost and number of people enrolling

TOTALCOSTFZ = zeros(3,1);
TOTALNUMFZ = zeros(3,1);
TOTALNUMFZi = zeros(3,nIs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% calculate choices made given preferences in year 1 
%%%%%%%%% Initial prices actually used by firm assumed as initial pricing state
%%%%%%%%% See section 6 in paper for additional discussion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nIs 
    PT = ones(nPlans,Sim);
    PT(1:nPlans-1,:) =  (reshape(EUvectorZ(i,3,:),1,Sim)'*ones(1,nPlans-1))'./reshape(EUvectorZ(i,[1:nPlans-1],:),nPlans-1,Sim);
    PTAgg = zeros(1,Sim);
    PTAgg = sum(PT,1);

    PT(1,:) = PT(1,:)./PTAgg;
    PT(2,:) = PT(2,:)./PTAgg;
    PT(3,:) = PT(3,:)./PTAgg;

    PT = ((PT'.^3)./(sum((PT'.^3),2)*ones(1,nPlans)))';
    
    UNI = unifrnd(0,1,1,Sim);
    
    for j = 1:Sim
        if PT(1,j) >= UNI(1,j)
           CC0(i,j)=1;
        elseif PT(1,j)+PT(2,j) >= UNI(1,j)
           CC0(i,j)=2;
        elseif (1-PT(3,j)) < UNI(1,j)
           CC0(i,j)=3;
        end
    end
    
end

for i = 1:nIs
    for j = 1:Sim
     
      if Tier1(i,1) == 1  
        if CC0(i,j) == 1 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject(i,j)==1))
           TOTALNUMEZi(1,i) = TOTALNUMEZi(1,i)+1; 
        end
        if CC0(i,j) == 2 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject(i,j)==1))
           TOTALNUMEZi(2,i) = TOTALNUMEZi(2,i)+1; 
        end
        if CC0(i,j) == 3 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject(i,j)==1))
           TOTALNUMEZi(3,i) = TOTALNUMEZi(3,i)+1; 
        end
      end
      
    if (Tier1(i,1) == 2 | Tier1(i,1) ==11 | Tier1(i,1)==12)  
        if CC0(i,j) == 1 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject(i,j)==1))
           TOTALNUMSZi(1,i) = TOTALNUMSZi(1,i)+1; 
        end
        if CC0(i,j) == 2 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject(i,j)==1))
           TOTALNUMSZi(2,i) = TOTALNUMSZi(2,i)+1; 
        end
        if CC0(i,j) == 3 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject(i,j)==1))
           TOTALNUMSZi(3,i) = TOTALNUMSZi(3,i)+1; 
        end
    end    

    if Tier1(i,1) == 6  
        if CC0(i,j) == 1 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject(i,j)==1))
           TOTALNUMCZi(1,i) = TOTALNUMCZi(1,i)+1; 
        end
        if CC0(i,j) == 2 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject(i,j)==1))
           TOTALNUMCZi(2,i) = TOTALNUMCZi(2,i)+1; 
        end
        if CC0(i,j) == 3 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject(i,j)==1))
           TOTALNUMCZi(3,i) = TOTALNUMCZi(3,i)+1; 
        end
    end   
    
   if Tier1(i,1) == 8  
        if CC0(i,j) == 1 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject(i,j)==1))
           TOTALNUMFZi(1,i) = TOTALNUMFZi(1,i)+1; 
        end
        if CC0(i,j) == 2 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject(i,j)==1))
           TOTALNUMFZi(2,i) = TOTALNUMFZi(2,i)+1; 
        end
        if CC0(i,j) == 3 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject(i,j)==1))
           TOTALNUMFZi(3,i) = TOTALNUMFZi(3,i)+1; 
        end
    end  
    end
end

for i = 1:nIs
	if Tier1(i,1)==1
    TOTALNUMEZ = TOTALNUMEZ + TOTALNUMEZi(:,i)*(Sim/sum(TOTALNUMEZi(:,i)));    
    end
    if (Tier1(i,1) == 2 | Tier1(i,1) ==11 | Tier1(i,1)==12)  
    TOTALNUMSZ = TOTALNUMSZ + TOTALNUMSZi(:,i)*(Sim/sum(TOTALNUMSZi(:,i)));
    end
    if Tier1(i,1)==6
    TOTALNUMCZ = TOTALNUMCZ + TOTALNUMCZi(:,i)*(Sim/sum(TOTALNUMCZi(:,i)));
    end
    if Tier1(i,1)==8
    TOTALNUMFZ = TOTALNUMFZ + TOTALNUMFZi(:,i)*(Sim/sum(TOTALNUMFZi(:,i)));
    end 
end


for i = 1:nIs  
    if Tier1(i,1) == 1
    TOTALCOSTEZ(1,1) = TOTALCOSTEZ(1,1) + (1/K)*sum(reshape(PlanPaid1(i,:,1),1,K),2)*TOTALNUMEZi(1,i)*(Sim/sum(TOTALNUMEZi(:,i)));
    TOTALCOSTEZ(2,1) = TOTALCOSTEZ(2,1) + (1/K)*sum(reshape(PlanPaid1(i,:,2),1,K),2)*TOTALNUMEZi(2,i)*(Sim/sum(TOTALNUMEZi(:,i)));
    TOTALCOSTEZ(3,1) = TOTALCOSTEZ(3,1) + (1/K)*sum(reshape(PlanPaid1(i,:,3),1,K),2)*TOTALNUMEZi(3,i)*(Sim/sum(TOTALNUMEZi(:,i)));
    end
   
    if (Tier1(i,1) == 2 | Tier1(i,1) ==11 | Tier1(i,1)==12)  
    TOTALCOSTSZ(1,1) = TOTALCOSTSZ(1,1) + (1/K)*sum(reshape(PlanPaid1(i,:,1),1,K),2)*TOTALNUMSZi(1,i)*(Sim/sum(TOTALNUMSZi(:,i)));
    TOTALCOSTSZ(2,1) = TOTALCOSTSZ(2,1) + (1/K)*sum(reshape(PlanPaid1(i,:,2),1,K),2)*TOTALNUMSZi(2,i)*(Sim/sum(TOTALNUMSZi(:,i)));
    TOTALCOSTSZ(3,1) = TOTALCOSTSZ(3,1) + (1/K)*sum(reshape(PlanPaid1(i,:,3),1,K),2)*TOTALNUMSZi(3,i)*(Sim/sum(TOTALNUMSZi(:,i)));
    end
    
    if Tier1(i,1)==6
    TOTALCOSTCZ(1,1) = TOTALCOSTCZ(1,1) + (1/K)*sum(reshape(PlanPaid1(i,:,1),1,K),2)*TOTALNUMCZi(1,i)*(Sim/sum(TOTALNUMCZi(:,i)));
    TOTALCOSTCZ(2,1) = TOTALCOSTCZ(2,1) + (1/K)*sum(reshape(PlanPaid1(i,:,2),1,K),2)*TOTALNUMCZi(2,i)*(Sim/sum(TOTALNUMCZi(:,i)));
    TOTALCOSTCZ(3,1) = TOTALCOSTCZ(3,1) + (1/K)*sum(reshape(PlanPaid1(i,:,3),1,K),2)*TOTALNUMCZi(3,i)*(Sim/sum(TOTALNUMCZi(:,i)));
    end    
       
   if Tier1(i,1)==8
   TOTALCOSTFZ(1,1) = TOTALCOSTFZ(1,1) + (1/K)*sum(reshape(PlanPaid1(i,:,1),1,K),2)*TOTALNUMFZi(1,i)*(Sim/sum(TOTALNUMFZi(:,i)));
   TOTALCOSTFZ(2,1) = TOTALCOSTFZ(2,1) + (1/K)*sum(reshape(PlanPaid1(i,:,2),1,K),2)*TOTALNUMFZi(2,i)*(Sim/sum(TOTALNUMFZi(:,i)));
   TOTALCOSTFZ(3,1) = TOTALCOSTFZ(3,1) + (1/K)*sum(reshape(PlanPaid1(i,:,3),1,K),2)*TOTALNUMFZi(3,i)*(Sim/sum(TOTALNUMFZi(:,i)));
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AVGCOSTEZ = zeros(3,1);
AVGCOSTSZ = zeros(3,1);
AVGCOSTCZ = zeros(3,1);
AVGCOSTFZ = zeros(3,1);
 
AVGCOSTEZ(1,1) = TOTALCOSTEZ(1,1)/TOTALNUMEZ(1,1);
AVGCOSTEZ(2,1) = TOTALCOSTEZ(2,1)/TOTALNUMEZ(2,1);
AVGCOSTEZ(3,1) = TOTALCOSTEZ(3,1)/TOTALNUMEZ(3,1);

AVGCOSTSZ(1,1) = TOTALCOSTSZ(1,1)/TOTALNUMSZ(1,1);
AVGCOSTSZ(2,1) = TOTALCOSTSZ(2,1)/TOTALNUMSZ(2,1);
AVGCOSTSZ(3,1) = TOTALCOSTSZ(3,1)/TOTALNUMSZ(3,1);

AVGCOSTCZ(1,1) = TOTALCOSTCZ(1,1)/TOTALNUMCZ(1,1);
AVGCOSTCZ(2,1) = TOTALCOSTCZ(2,1)/TOTALNUMCZ(2,1);
AVGCOSTCZ(3,1) = TOTALCOSTCZ(3,1)/TOTALNUMCZ(3,1);

AVGCOSTFZ(1,1) = TOTALCOSTFZ(1,1)/TOTALNUMFZ(1,1);
AVGCOSTFZ(2,1) = TOTALCOSTFZ(2,1)/TOTALNUMFZ(2,1);
AVGCOSTFZ(3,1) = TOTALCOSTFZ(3,1)/TOTALNUMFZ(3,1);

MktShare0 = zeros(3,1);

MktShare0(1,1) = (TOTALNUMFZ(1,1) + TOTALNUMCZ(1,1) + TOTALNUMSZ(1,1) + TOTALNUMEZ(1,1))/Sim;
MktShare0(2,1) = (TOTALNUMFZ(2,1) + TOTALNUMCZ(2,1) + TOTALNUMSZ(2,1) + TOTALNUMEZ(2,1))/Sim;
MktShare0(3,1) = (TOTALNUMFZ(3,1) + TOTALNUMCZ(3,1) + TOTALNUMSZ(3,1) + TOTALNUMEZ(3,1))/Sim;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Flat subsidy factor for PPO1200 as a function of income status, as given at firm
%%%%% see section 6 in paper for more details
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sub1 = 0.97;
Sub2 = 0.93;
Sub3 = 0.83;
Sub4 = 0.71;
Sub5 = 0.64;

%%%%%%% Premium 'pass-through' parameter. This is 1 in general in our environment, but can shift this
%%%%%%% to incorporate more or less marginal pass through as desired. Equals 1 in environment in general
%%%%%%% which is what is used for main analysis. Z can also be > 1 if you want to consider is an 
%%%%%%% administrative loading factor that is incorporated into analysis. Here, we let this be 1.10 (so administrative = 10% of total costs) 
%%%%%%% but it can be easily varied to examine different assumptions on this. In actual analysis set = 1.10, 10% administrative costs.  

Z = 1.1;
Z1 = 1.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Now I want to construct prices for year 2 given year 1 enrollment %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% See section 6 for more details. Prices depend on average costs of past enrollment 
%%%% conditional on family status. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Premiums multipled by 1-Marginal Tax Rate since premiums are tax deductible / tax free 
%%%%%%%%%%% This is number at end, differs by income tier. 

for i = 1:nIs

   if Inc1(i,1) == 1
      if Tier1(i,1) == 1
          P1(2,i,:)= (AVGCOSTEZ(1,1)-AVGCOSTEZ(3,1)*Sub1)*0.82*Z;
          P2(2,i,:)= (AVGCOSTEZ(2,1)-AVGCOSTEZ(3,1)*Sub1)*0.82*Z;
          P3(2,i,:)= AVGCOSTEZ(3,1)*0.82*(1-Sub1)*Z;          
      end
      if (Tier1(i,1) == 2 | Tier1(i,1) ==11 | Tier1(i,1)==12)
          P1(2,i,:)= (AVGCOSTSZ(1,1)-AVGCOSTSZ(3,1)*Sub1)*0.82*Z;
          P2(2,i,:)= (AVGCOSTSZ(2,1)-AVGCOSTSZ(3,1)*Sub1)*0.82*Z;
          P3(2,i,:)= AVGCOSTSZ(3,1)*0.82*(1-Sub1)*Z;          
      end
      if Tier1(i,1) == 6
          P1(2,i,:)= (AVGCOSTCZ(1,1)-AVGCOSTCZ(3,1)*Sub1)*0.82*Z1;
          P2(2,i,:)= (AVGCOSTCZ(2,1)-AVGCOSTCZ(3,1)*Sub1)*0.82*Z1;
          P3(2,i,:)= AVGCOSTCZ(3,1)*0.82*(1-Sub1)*Z1;          
      end      
      if Tier1(i,1) == 8
          P1(2,i,:)= (AVGCOSTFZ(1,1)-AVGCOSTFZ(3,1)*Sub1)*0.82*Z;
          P2(2,i,:)= (AVGCOSTFZ(2,1)-AVGCOSTFZ(3,1)*Sub1)*0.82*Z;
          P3(2,i,:)= AVGCOSTFZ(3,1)*0.82*(1-Sub1)*Z;          
      end
   end
   
   if Inc1(i,1) == 2
      if Tier1(i,1) == 1
          P1(2,i,:)= (AVGCOSTEZ(1,1)-AVGCOSTEZ(3,1)*Sub2)*0.73*Z;
          P2(2,i,:)= (AVGCOSTEZ(2,1)-AVGCOSTEZ(3,1)*Sub2)*0.73*Z;
          P3(2,i,:)= AVGCOSTEZ(3,1)*0.73*(1-Sub2)*Z;          
      end
      if (Tier1(i,1) == 2 | Tier1(i,1) ==11 | Tier1(i,1)==12)
          P1(2,i,:)= (AVGCOSTSZ(1,1)-AVGCOSTSZ(3,1)*Sub2)*0.80*Z;
          P2(2,i,:)= (AVGCOSTSZ(2,1)-AVGCOSTSZ(3,1)*Sub2)*0.80*Z;
          P3(2,i,:)= AVGCOSTSZ(3,1)*0.80*(1-Sub2)*Z;          
      end
      if Tier1(i,1) == 6
          P1(2,i,:)= (AVGCOSTCZ(1,1)-AVGCOSTCZ(3,1)*Sub2)*0.80*Z1;
          P2(2,i,:)= (AVGCOSTCZ(2,1)-AVGCOSTCZ(3,1)*Sub2)*0.80*Z1;
          P3(2,i,:)= AVGCOSTCZ(3,1)*0.80*(1-Sub2)*Z1;          
      end      
      if Tier1(i,1) == 8
          P1(2,i,:)= (AVGCOSTFZ(1,1)-AVGCOSTFZ(3,1)*Sub2)*0.80*Z;
          P2(2,i,:)= (AVGCOSTFZ(2,1)-AVGCOSTFZ(3,1)*Sub2)*0.80*Z;
          P3(2,i,:)= AVGCOSTFZ(3,1)*0.80*(1-Sub2)*Z;          
      end
   end
   if Inc1(i,1) == 3
      if Tier1(i,1) == 1
          P1(2,i,:)= (AVGCOSTEZ(1,1)-AVGCOSTEZ(3,1)*Sub3)*0.69*Z;
          P2(2,i,:)= (AVGCOSTEZ(2,1)-AVGCOSTEZ(3,1)*Sub3)*0.69*Z;
          P3(2,i,:)= AVGCOSTEZ(3,1)*0.69*(1-Sub3)*Z;          
      end
      if (Tier1(i,1) == 2 | Tier1(i,1) ==11 | Tier1(i,1)==12)
          P1(2,i,:)= (AVGCOSTSZ(1,1)-AVGCOSTSZ(3,1)*Sub3)*0.69*Z;
          P2(2,i,:)= (AVGCOSTSZ(2,1)-AVGCOSTSZ(3,1)*Sub3)*0.69*Z;
          P3(2,i,:)= AVGCOSTSZ(3,1)*0.69*(1-Sub3)*Z;          
      end
      if Tier1(i,1) == 6
          P1(2,i,:)= (AVGCOSTCZ(1,1)-AVGCOSTCZ(3,1)*Sub3)*0.69*Z1;
          P2(2,i,:)= (AVGCOSTCZ(2,1)-AVGCOSTCZ(3,1)*Sub3)*0.69*Z1;
          P3(2,i,:)= AVGCOSTCZ(3,1)*0.69*(1-Sub3)*Z1;          
      end      
      if Tier1(i,1) == 8
          P1(2,i,:)= (AVGCOSTFZ(1,1)-AVGCOSTFZ(3,1)*Sub3)*0.69*Z;
          P2(2,i,:)= (AVGCOSTFZ(2,1)-AVGCOSTFZ(3,1)*Sub3)*0.69*Z;
          P3(2,i,:)= AVGCOSTFZ(3,1)*0.69*(1-Sub3)*Z;          
      end
   end
   
   if Inc1(i,1) == 4
      if Tier1(i,1) == 1
          P1(2,i,:)= (AVGCOSTEZ(1,1)-AVGCOSTEZ(3,1)*Sub4)*0.66*Z;
          P2(2,i,:)= (AVGCOSTEZ(2,1)-AVGCOSTEZ(3,1)*Sub4)*0.66*Z;
          P3(2,i,:)= AVGCOSTEZ(3,1)*0.66*(1-Sub4)*Z;          
      end
      if (Tier1(i,1) == 2 | Tier1(i,1) ==11 | Tier1(i,1)==12)
          P1(2,i,:)= (AVGCOSTSZ(1,1)-AVGCOSTSZ(3,1)*Sub4)*0.66*Z;
          P2(2,i,:)= (AVGCOSTSZ(2,1)-AVGCOSTSZ(3,1)*Sub4)*0.66*Z;
          P3(2,i,:)= AVGCOSTSZ(3,1)*0.66*(1-Sub4)*Z;          
      end
      if Tier1(i,1) == 6
          P1(2,i,:)= (AVGCOSTCZ(1,1)-AVGCOSTCZ(3,1)*Sub4)*0.66*Z1;
          P2(2,i,:)= (AVGCOSTCZ(2,1)-AVGCOSTCZ(3,1)*Sub4)*0.66*Z1;
          P3(2,i,:)= AVGCOSTCZ(3,1)*0.66*(1-Sub4)*Z1;          
      end      
      if Tier1(i,1) == 8
          P1(2,i,:)= (AVGCOSTFZ(1,1)-AVGCOSTFZ(3,1)*Sub4)*0.66*Z;
          P2(2,i,:)= (AVGCOSTFZ(2,1)-AVGCOSTFZ(3,1)*Sub4)*0.66*Z;
          P3(2,i,:)= AVGCOSTFZ(3,1)*0.66*(1-Sub4)*Z;          
      end
   end   
   
   if Inc1(i,1) == 5
      if Tier1(i,1) == 1
          P1(2,i,:)= (AVGCOSTEZ(1,1)-AVGCOSTEZ(3,1)*Sub5)*0.61*Z;
          P2(2,i,:)= (AVGCOSTEZ(2,1)-AVGCOSTEZ(3,1)*Sub5)*0.61*Z;
          P3(2,i,:)= AVGCOSTEZ(3,1)*0.61*(1-Sub5)*Z;          
      end
      if (Tier1(i,1) == 2 | Tier1(i,1) ==11 | Tier1(i,1)==12)
          P1(2,i,:)= (AVGCOSTSZ(1,1)-AVGCOSTSZ(3,1)*Sub5)*0.61*Z;
          P2(2,i,:)= (AVGCOSTSZ(2,1)-AVGCOSTSZ(3,1)*Sub5)*0.61*Z;
          P3(2,i,:)= AVGCOSTSZ(3,1)*0.61*(1-Sub5)*Z;          
      end
      if Tier1(i,1) == 6
          P1(2,i,:)= (AVGCOSTCZ(1,1)-AVGCOSTCZ(3,1)*Sub5)*0.61*Z1;
          P2(2,i,:)= (AVGCOSTCZ(2,1)-AVGCOSTCZ(3,1)*Sub5)*0.61*Z1;
          P3(2,i,:)= AVGCOSTCZ(3,1)*0.61*(1-Sub5)*Z1;          
      end      
      if Tier1(i,1) == 8
          P1(2,i,:)= (AVGCOSTFZ(1,1)-AVGCOSTFZ(3,1)*Sub5)*0.61*Z;
          P2(2,i,:)= (AVGCOSTFZ(2,1)-AVGCOSTFZ(3,1)*Sub5)*0.61*Z;
          P3(2,i,:)= AVGCOSTFZ(3,1)*0.61*(1-Sub5)*Z;          
      end
   end    
 
end

%%%%%%% Set premiums equal to maximum of value computed above and 150, which is minimum observed in data (e.g. they never go below this level) 

P1(2,:,:) = max(reshape(P1(2,:,:),nIs,Sim),150);
P2(2,:,:) = max(reshape(P2(2,:,:),nIs,Sim),150);
P3(2,:,:) = max(reshape(P3(2,:,:),nIs,Sim),150);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Prices to display at end, calculation same as above, except
%%%%%%%%%%% above is filled in for people in large matrices to be used in calculations
%%%%%%%%%%% here calculated for final display. Also, tax benefit not taken into account, 
%%%%%%%%%%% these are actual pre-tax premiums charged. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    P3EZ = zeros(3,5);
    P3SZ = zeros(3,5);
    P3CZ = zeros(3,5);
    P3FZ = zeros(3,5);
    
    P3EZ(1,1) = (AVGCOSTEZ(1,1)-AVGCOSTEZ(3,1)*Sub1)*Z;
    P3EZ(2,1) = (AVGCOSTEZ(2,1)-AVGCOSTEZ(3,1)*Sub1)*Z;
    P3EZ(3,1) =  AVGCOSTEZ(3,1)*(1-Sub1)*Z;          
      
    P3SZ(1,1) = (AVGCOSTSZ(1,1)-AVGCOSTSZ(3,1)*Sub1)*Z;
    P3SZ(2,1) = (AVGCOSTSZ(2,1)-AVGCOSTSZ(3,1)*Sub1)*Z;
    P3SZ(3,1) =  AVGCOSTSZ(3,1)*(1-Sub1)*Z;          
   
    P3CZ(1,1) = (AVGCOSTCZ(1,1)-AVGCOSTCZ(3,1)*Sub1)*Z1;
    P3CZ(2,1) = (AVGCOSTCZ(2,1)-AVGCOSTCZ(3,1)*Sub1)*Z1;
    P3CZ(3,1) =  AVGCOSTCZ(3,1)*(1-Sub1)*Z1;          
    
    P3FZ(1,1) = (AVGCOSTFZ(1,1)-AVGCOSTFZ(3,1)*Sub1)*Z;
    P3FZ(2,1) = (AVGCOSTFZ(2,1)-AVGCOSTFZ(3,1)*Sub1)*Z;
    P3FZ(3,1) =  AVGCOSTFZ(3,1)*(1-Sub1)*Z;          
  
    P3EZ(1,2) = (AVGCOSTEZ(1,1)-AVGCOSTEZ(3,1)*Sub2)*Z;
    P3EZ(2,2) = (AVGCOSTEZ(2,1)-AVGCOSTEZ(3,1)*Sub2)*Z;
    P3EZ(3,2) =  AVGCOSTEZ(3,1)*(1-Sub2)*Z;          
      
    P3SZ(1,2) = (AVGCOSTSZ(1,1)-AVGCOSTSZ(3,1)*Sub2)*Z;
    P3SZ(2,2) = (AVGCOSTSZ(2,1)-AVGCOSTSZ(3,1)*Sub2)*Z;
    P3SZ(3,2) = AVGCOSTSZ(3,1)*(1-Sub2)*Z;           
   
    P3CZ(1,2) = (AVGCOSTCZ(1,1)-AVGCOSTCZ(3,1)*Sub2)*Z1;
    P3CZ(2,2) = (AVGCOSTCZ(2,1)-AVGCOSTCZ(3,1)*Sub2)*Z1;
    P3CZ(3,2) =  AVGCOSTCZ(3,1)*(1-Sub2)*Z1;         
    
    P3FZ(1,2) = (AVGCOSTFZ(1,1)-AVGCOSTFZ(3,1)*Sub2)*Z;
    P3FZ(2,2) = (AVGCOSTFZ(2,1)-AVGCOSTFZ(3,1)*Sub2)*Z;
    P3FZ(3,2) =  AVGCOSTFZ(3,1)*(1-Sub2)*Z;      

    P3EZ(1,3) = (AVGCOSTEZ(1,1)-AVGCOSTEZ(3,1)*Sub3)*Z;
    P3EZ(2,3) = (AVGCOSTEZ(2,1)-AVGCOSTEZ(3,1)*Sub3)*Z;
    P3EZ(3,3) =  AVGCOSTEZ(3,1)*(1-Sub3)*Z;         
      
    P3SZ(1,3) = (AVGCOSTSZ(1,1)-AVGCOSTSZ(3,1)*Sub3)*Z;
    P3SZ(2,3) = (AVGCOSTSZ(2,1)-AVGCOSTSZ(3,1)*Sub3)*Z;
    P3SZ(3,3) =  AVGCOSTSZ(3,1)*(1-Sub3)*Z;          
   
    P3CZ(1,3) = (AVGCOSTCZ(1,1)-AVGCOSTCZ(3,1)*Sub3)*Z1;
    P3CZ(2,3) = (AVGCOSTCZ(2,1)-AVGCOSTCZ(3,1)*Sub3)*Z1;
    P3CZ(3,3) = AVGCOSTCZ(3,1)*(1-Sub3)*Z1;        
    
    P3FZ(1,3) = (AVGCOSTFZ(1,1)-AVGCOSTFZ(3,1)*Sub3)*Z;
    P3FZ(2,3) = (AVGCOSTFZ(2,1)-AVGCOSTFZ(3,1)*Sub3)*Z;
    P3FZ(3,3) = AVGCOSTFZ(3,1)*(1-Sub3)*Z; 
            
    P3EZ(1,4) = (AVGCOSTEZ(1,1)-AVGCOSTEZ(3,1)*Sub4)*Z;
    P3EZ(2,4) = (AVGCOSTEZ(2,1)-AVGCOSTEZ(3,1)*Sub4)*Z;
    P3EZ(3,4) = AVGCOSTEZ(3,1)*(1-Sub4)*Z;         
      
    P3SZ(1,4) = (AVGCOSTSZ(1,1)-AVGCOSTSZ(3,1)*Sub4)*Z;
    P3SZ(2,4) = (AVGCOSTSZ(2,1)-AVGCOSTSZ(3,1)*Sub4)*Z;
    P3SZ(3,4) =  AVGCOSTSZ(3,1)*(1-Sub4)*Z;           
   
    P3CZ(1,4) = (AVGCOSTCZ(1,1)-AVGCOSTCZ(3,1)*Sub4)*Z1;
    P3CZ(2,4) = (AVGCOSTCZ(2,1)-AVGCOSTCZ(3,1)*Sub4)*Z1;
    P3CZ(3,4) =  AVGCOSTCZ(3,1)*(1-Sub4)*Z1;        
    
    P3FZ(1,4) = (AVGCOSTFZ(1,1)-AVGCOSTFZ(3,1)*Sub4)*Z;
    P3FZ(2,4) = (AVGCOSTFZ(2,1)-AVGCOSTFZ(3,1)*Sub4)*Z;
    P3FZ(3,4) = AVGCOSTFZ(3,1)*(1-Sub4)*Z;

    P3EZ(1,5) = (AVGCOSTEZ(1,1)-AVGCOSTEZ(3,1)*Sub5)*Z;
    P3EZ(2,5) = (AVGCOSTEZ(2,1)-AVGCOSTEZ(3,1)*Sub5)*Z;
    P3EZ(3,5) =  AVGCOSTEZ(3,1)*(1-Sub5)*Z;        
      
    P3SZ(1,5) = (AVGCOSTSZ(1,1)-AVGCOSTSZ(3,1)*Sub5)*Z;
    P3SZ(2,5) = (AVGCOSTSZ(2,1)-AVGCOSTSZ(3,1)*Sub5)*Z;
    P3SZ(3,5) = AVGCOSTSZ(3,1)*(1-Sub5)*Z;             
   
    P3CZ(1,5) = (AVGCOSTCZ(1,1)-AVGCOSTCZ(3,1)*Sub5)*Z1;
    P3CZ(2,5) = (AVGCOSTCZ(2,1)-AVGCOSTCZ(3,1)*Sub5)*Z1;
    P3CZ(3,5) =  AVGCOSTCZ(3,1)*(1-Sub5)*Z1;       
    
    P3FZ(1,5) = (AVGCOSTFZ(1,1)-AVGCOSTFZ(3,1)*Sub5)*Z;
    P3FZ(2,5) = (AVGCOSTFZ(2,1)-AVGCOSTFZ(3,1)*Sub5)*Z;
    P3FZ(3,5) = AVGCOSTFZ(3,1)*(1-Sub5)*Z;

    
P3FZ = max(P3FZ,150);
P3SZ = max(P3SZ,150);
P3CZ = max(P3CZ,150);
P3EZ = max(P3EZ,150);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Develop incumbent plan indicators for next year for simulated choices
%%%%%% incorporate with reduced inertia value so this impact utility / choices
%%%%%% as prices / choices change over time in counterfactual environment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1:nIs
    for s = 1:Sim
       if CC0(i,s)==1
        choice12(i,s) = 1;
       end
       if CC0(i,s)==2
        choice22(i,s) = 1;
       end
       if CC0(i,s)==3
        choice32(i,s) = 1;
       end
    end
end

for k = 1:K
    SwC12(:,k,:) = SCC2.*choice12;
    SwC22(:,k,:) = SCC2.*choice22;
    SwC32(:,k,:) = SCC2.*choice32;
end

%%%% New expected utility baselines calculated because of new premiums and inertia 

for i = 1:nIs
    for k = 1:K
        EUvector12(i,k,:) = max(0.1,(W - reshape(PPO250OP2(i,k,:),1,Sim) - reshape(P1(2,i,:),1,Sim)));
        EUvector22(i,k,:) = max(0.1,(W - reshape(PPO500OP2(i,k,:),1,Sim) - reshape(P2(2,i,:),1,Sim)));
        EUvector32(i,k,:) = max(0.1,(W - reshape(PPO1200OP2(i,k,:),1,Sim) - reshape(P3(2,i,:),1,Sim)));
    end
end

EUvector12 = exp(-((EUvector12 +SwC12).*max(RA,0.0000001)));
EUvector22 = exp(-((EUvector22 + SwC22 + CHTC + Epsilon5002).*max(RA,0.0000001)));
EUvector32 = exp(-((EUvector32 + CDHP + CHTC + SwC32 + Epsilon12002).*max(RA, 0.0000001)));

for n = 1:Sim
    EUvectorAx(:,1,n) = (1/K)*sum(reshape(EUvector12(:,:,n),nIs,K),2);
    EUvectorAx(:,2,n) = (1/K)*sum(reshape(EUvector22(:,:,n),nIs,K),2);
    EUvectorAx(:,3,n) = (1/K)*sum(reshape(EUvector32(:,:,n),nIs,K),2);
end

%%%%%%%%%%% Simulate year 2 choices given new prices and choices
%%%%%%%%%%% Process for this point forward same as for year 1, so see year 1 above for 
%%%%%%%%%%% commented code. 

CC1 = zeros(nIs,Sim);
MktShare1 = zeros(nPlans,1);

TOTALCOSTE = zeros(3,1);
TOTALNUME = zeros(3,1);
TOTALNUMEi = zeros(3,nIs);

TOTALCOSTS = zeros(3,1);
TOTALNUMS = zeros(3,1);
TOTALNUMSi = zeros(3,nIs);

TOTALCOSTC = zeros(3,1);
TOTALNUMC = zeros(3,1);
TOTALNUMCi = zeros(3,nIs);

TOTALCOSTF = zeros(3,1);
TOTALNUMF = zeros(3,1);
TOTALNUMFi = zeros(3,nIs);


for i = 1:nIs 
    PT2 = ones(nPlans,Sim);
    PT2(1:nPlans-1,:) =  (reshape(EUvectorAx(i,3,:),1,Sim)'*ones(1,nPlans-1))'./reshape(EUvectorAx(i,[1:nPlans-1],:),nPlans-1,Sim);
    PTAgg2 = zeros(1,Sim);
    PTAgg2 = sum(PT2,1);
    
    PT2(1,:) = PT2(1,:)./PTAgg2;
    PT2(2,:) = PT2(2,:)./PTAgg2;
    PT2(3,:) = PT2(3,:)./PTAgg2;

    PT2 = ((PT2'.^3)./(sum((PT2'.^3),2)*ones(1,nPlans)))';
    
    UNI = unifrnd(0,1,1,Sim);
    
    for j = 1:Sim
        if PT2(1,j) >= UNI(1,j)
           CC1(i,j)=1;
        elseif PT2(1,j)+PT2(2,j) >= UNI(1,j)
           CC1(i,j)=2;
        elseif (1-PT2(3,j)) < UNI(1,j)
           CC1(i,j)=3;
        end
    end
    
end

for i = 1:nIs
    for j = 1:Sim      
      if Tier2(i,1) == 1  
        if CC1(i,j) == 1 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject2(i,j)==1))
           TOTALNUMEi(1,i) = TOTALNUMEi(1,i)+1; 
        end
        if CC1(i,j) == 2 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject2(i,j)==1))
           TOTALNUMEi(2,i) = TOTALNUMEi(2,i)+1; 
        end
        if CC1(i,j) == 3 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject2(i,j)==1))
           TOTALNUMEi(3,i) = TOTALNUMEi(3,i)+1; 
        end
      end
      
    if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)  
        if CC1(i,j) == 1 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject2(i,j)==1))
           TOTALNUMSi(1,i) = TOTALNUMSi(1,i)+1; 
        end
        if CC1(i,j) == 2 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject2(i,j)==1))
           TOTALNUMSi(2,i) = TOTALNUMSi(2,i)+1; 
        end
        if CC1(i,j) == 3 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject2(i,j)==1))
           TOTALNUMSi(3,i) = TOTALNUMSi(3,i)+1; 
        end
    end    

    if Tier2(i,1) == 6  
        if CC1(i,j) == 1 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject2(i,j)==1))
           TOTALNUMCi(1,i) = TOTALNUMCi(1,i)+1; 
        end
        if CC1(i,j) == 2 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject2(i,j)==1))
           TOTALNUMCi(2,i) = TOTALNUMCi(2,i)+1; 
        end
        if CC1(i,j) == 3 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject2(i,j)==1))
           TOTALNUMCi(3,i) = TOTALNUMCi(3,i)+1; 
        end
    end   
    
   if Tier2(i,1) == 8  
        if CC1(i,j) == 1 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject2(i,j)==1))
           TOTALNUMFi(1,i) = TOTALNUMFi(1,i)+1; 
        end
        if CC1(i,j) == 2 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject2(i,j)==1))
           TOTALNUMFi(2,i) = TOTALNUMFi(2,i)+1; 
        end
        if CC1(i,j) == 3 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject2(i,j)==1))
           TOTALNUMFi(3,i) = TOTALNUMFi(3,i)+1; 
        end
    end  
    end
end

for i = 1:nIs 
    if Tier2(i,1)==1
    TOTALNUME = TOTALNUME + TOTALNUMEi(:,i)*(Sim/sum(TOTALNUMEi(:,i)));    
    end
    if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)  
    TOTALNUMS = TOTALNUMS + TOTALNUMSi(:,i)*(Sim/sum(TOTALNUMSi(:,i)));
    end
    if Tier2(i,1)==6
    TOTALNUMC = TOTALNUMC + TOTALNUMCi(:,i)*(Sim/sum(TOTALNUMCi(:,i)));
    end
    if Tier2(i,1)==8
    TOTALNUMF = TOTALNUMF + TOTALNUMFi(:,i)*(Sim/sum(TOTALNUMFi(:,i)));
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nIs
   if Tier2(i,1) == 1
    TOTALCOSTE(1,1) = TOTALCOSTE(1,1) + (1/K)*sum(reshape(PlanPaid2(i,:,1),1,K),2)*TOTALNUMEi(1,i)*(Sim/sum(TOTALNUMEi(:,i)));
    TOTALCOSTE(2,1) = TOTALCOSTE(2,1) + (1/K)*sum(reshape(PlanPaid2(i,:,2),1,K),2)*TOTALNUMEi(2,i)*(Sim/sum(TOTALNUMEi(:,i)));
    TOTALCOSTE(3,1) = TOTALCOSTE(3,1) + (1/K)*sum(reshape(PlanPaid2(i,:,3),1,K),2)*TOTALNUMEi(3,i)*(Sim/sum(TOTALNUMEi(:,i)));
    end
   
    if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)  
    TOTALCOSTS(1,1) = TOTALCOSTS(1,1) + (1/K)*sum(reshape(PlanPaid2(i,:,1),1,K),2)*TOTALNUMSi(1,i)*(Sim/sum(TOTALNUMSi(:,i)));
    TOTALCOSTS(2,1) = TOTALCOSTS(2,1) + (1/K)*sum(reshape(PlanPaid2(i,:,2),1,K),2)*TOTALNUMSi(2,i)*(Sim/sum(TOTALNUMSi(:,i)));
    TOTALCOSTS(3,1) = TOTALCOSTS(3,1) + (1/K)*sum(reshape(PlanPaid2(i,:,3),1,K),2)*TOTALNUMSi(3,i)*(Sim/sum(TOTALNUMSi(:,i)));
    end
    
    if Tier2(i,1)==6
    TOTALCOSTC(1,1) = TOTALCOSTC(1,1) + (1/K)*sum(reshape(PlanPaid2(i,:,1),1,K),2)*TOTALNUMCi(1,i)*(Sim/sum(TOTALNUMCi(:,i)));
    TOTALCOSTC(2,1) = TOTALCOSTC(2,1) + (1/K)*sum(reshape(PlanPaid2(i,:,2),1,K),2)*TOTALNUMCi(2,i)*(Sim/sum(TOTALNUMCi(:,i)));
    TOTALCOSTC(3,1) = TOTALCOSTC(3,1) + (1/K)*sum(reshape(PlanPaid2(i,:,3),1,K),2)*TOTALNUMCi(3,i)*(Sim/sum(TOTALNUMCi(:,i)));
    end    
       
   if Tier2(i,1)==8
   TOTALCOSTF(1,1) = TOTALCOSTF(1,1) + (1/K)*sum(reshape(PlanPaid2(i,:,1),1,K),2)*TOTALNUMFi(1,i)*(Sim/sum(TOTALNUMFi(:,i)));
   TOTALCOSTF(2,1) = TOTALCOSTF(2,1) + (1/K)*sum(reshape(PlanPaid2(i,:,2),1,K),2)*TOTALNUMFi(2,i)*(Sim/sum(TOTALNUMFi(:,i)));
   TOTALCOSTF(3,1) = TOTALCOSTF(3,1) + (1/K)*sum(reshape(PlanPaid2(i,:,3),1,K),2)*TOTALNUMFi(3,i)*(Sim/sum(TOTALNUMFi(:,i)));
   end
end

AVGCOSTE = zeros(3,1);
AVGCOSTS = zeros(3,1);
AVGCOSTC = zeros(3,1);
AVGCOSTF = zeros(3,1);

AVGCOSTE(1,1) = TOTALCOSTE(1,1)/TOTALNUME(1,1);
AVGCOSTE(2,1) = TOTALCOSTE(2,1)/TOTALNUME(2,1);
AVGCOSTE(3,1) = TOTALCOSTE(3,1)/TOTALNUME(3,1);

AVGCOSTS(1,1) = TOTALCOSTS(1,1)/TOTALNUMS(1,1);
AVGCOSTS(2,1) = TOTALCOSTS(2,1)/TOTALNUMS(2,1);
AVGCOSTS(3,1) = TOTALCOSTS(3,1)/TOTALNUMS(3,1);

AVGCOSTC(1,1) = TOTALCOSTC(1,1)/TOTALNUMC(1,1);
AVGCOSTC(2,1) = TOTALCOSTC(2,1)/TOTALNUMC(2,1);
AVGCOSTC(3,1) = TOTALCOSTC(3,1)/TOTALNUMC(3,1);

AVGCOSTF(1,1) = TOTALCOSTF(1,1)/TOTALNUMF(1,1);
AVGCOSTF(2,1) = TOTALCOSTF(2,1)/TOTALNUMF(2,1);
AVGCOSTF(3,1) = TOTALCOSTF(3,1)/TOTALNUMF(3,1);

MktShare1 = zeros(3,1);

MktShare1(1,1) = (TOTALNUMF(1,1) + TOTALNUMC(1,1) + TOTALNUMS(1,1) + TOTALNUME(1,1))/Sim;
MktShare1(2,1) = (TOTALNUMF(2,1) + TOTALNUMC(2,1) + TOTALNUMS(2,1) + TOTALNUME(2,1))/Sim;
MktShare1(3,1) = (TOTALNUMF(3,1) + TOTALNUMC(3,1) + TOTALNUMS(3,1) + TOTALNUME(3,1))/Sim;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NOTE: Subsidies consistent from year to year! %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sub1 = 0.97;
Sub2 = 0.93;
Sub3 = 0.83;
Sub4 = 0.71;
Sub5 = 0.64;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Now I want to construct prices for year 3          %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nIs

   if Inc2(i,1) == 1
      if Tier2(i,1) == 1
          P1(3,i,:)= (AVGCOSTE(1,1)-AVGCOSTE(3,1)*Sub1)*0.82*Z;
          P2(3,i,:)= (AVGCOSTE(2,1)-AVGCOSTE(3,1)*Sub1)*0.82*Z;
          P3(3,i,:)= AVGCOSTE(3,1)*0.82*(1-Sub1)*Z;          
      end
      if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)
          P1(3,i,:)= (AVGCOSTS(1,1)*Z-AVGCOSTS(3,1)*Sub1)*0.82*Z;
          P2(3,i,:)= (AVGCOSTS(2,1)*Z-AVGCOSTS(3,1)*Sub1)*0.82*Z;
          P3(3,i,:)= AVGCOSTS(3,1)*0.82*(1-Sub1)*Z;          
      end
      if Tier2(i,1) == 6
          P1(3,i,:)= (AVGCOSTC(1,1)*Z-AVGCOSTC(3,1)*Sub1)*0.82*Z1;
          P2(3,i,:)= (AVGCOSTC(2,1)*Z-AVGCOSTC(3,1)*Sub1)*0.82*Z1;
          P3(3,i,:)= AVGCOSTC(3,1)*Z*0.82*(1-Sub1)*Z1;          
      end      
      if Tier2(i,1) == 8
          P1(3,i,:)= (AVGCOSTF(1,1)-AVGCOSTF(3,1)*Sub1)*0.82*Z;
          P2(3,i,:)= (AVGCOSTF(2,1)-AVGCOSTF(3,1)*Sub1)*0.82*Z;
          P3(3,i,:)= AVGCOSTF(3,1)*0.82*(1-Sub1)*Z;          
      end
   end
   
   if Inc2(i,1) == 2
      if Tier2(i,1) == 1
          P1(3,i,:)= (AVGCOSTE(1,1)-AVGCOSTE(3,1)*Sub2)*0.73*Z;
          P2(3,i,:)= (AVGCOSTE(2,1)-AVGCOSTE(3,1)*Sub2)*0.73*Z;
          P3(3,i,:)= AVGCOSTE(3,1)*0.73*(1-Sub2)*Z;          
      end
      if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)
          P1(3,i,:)= (AVGCOSTS(1,1)-AVGCOSTS(3,1)*Sub2)*0.80*Z;
          P2(3,i,:)= (AVGCOSTS(2,1)-AVGCOSTS(3,1)*Sub2)*0.80*Z;
          P3(3,i,:)= AVGCOSTS(3,1)*0.80*(1-Sub2)*Z;          
      end
      if Tier2(i,1) == 6
          P1(3,i,:)= (AVGCOSTC(1,1)-AVGCOSTC(3,1)*Sub2)*0.80*Z1;
          P2(3,i,:)= (AVGCOSTC(2,1)-AVGCOSTC(3,1)*Sub2)*0.80*Z1;
          P3(3,i,:)= AVGCOSTC(3,1)*0.80*(1-Sub2)*Z1;          
      end      
      if Tier2(i,1) == 8
          P1(3,i,:)= (AVGCOSTF(1,1)-AVGCOSTF(3,1)*Sub2)*0.80*Z;
          P2(3,i,:)= (AVGCOSTF(2,1)-AVGCOSTF(3,1)*Sub2)*0.80*Z;
          P3(3,i,:)= AVGCOSTF(3,1)*0.80*(1-Sub2)*Z;          
      end
   end
   if Inc2(i,1) == 3
      if Tier2(i,1) == 1
          P1(3,i,:)= (AVGCOSTE(1,1)-AVGCOSTE(3,1)*Sub3)*0.69*Z;
          P2(3,i,:)= (AVGCOSTE(2,1)-AVGCOSTE(3,1)*Sub3)*0.69*Z;
          P3(3,i,:)= AVGCOSTE(3,1)*0.69*(1-Sub3)*Z;          
      end
      if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)
          P1(3,i,:)= (AVGCOSTS(1,1)-AVGCOSTS(3,1)*Sub3)*0.69*Z;
          P2(3,i,:)= (AVGCOSTS(2,1)-AVGCOSTS(3,1)*Sub3)*0.69*Z;
          P3(3,i,:)= AVGCOSTS(3,1)*0.69*(1-Sub3)*Z;          
      end
      if Tier2(i,1) == 6
          P1(3,i,:)= (AVGCOSTC(1,1)-AVGCOSTC(3,1)*Sub3)*0.69*Z1;
          P2(3,i,:)= (AVGCOSTC(2,1)-AVGCOSTC(3,1)*Sub3)*0.69*Z1;
          P3(3,i,:)= AVGCOSTC(3,1)*0.69*(1-Sub3)*Z1;          
      end      
      if Tier2(i,1) == 8
          P1(3,i,:)= (AVGCOSTF(1,1)-AVGCOSTF(3,1)*Sub3)*0.69*Z;
          P2(3,i,:)= (AVGCOSTF(2,1)-AVGCOSTF(3,1)*Sub3)*0.69*Z;
          P3(3,i,:)= AVGCOSTF(3,1)*0.69*(1-Sub3)*Z;          
      end
   end
   
   if Inc2(i,1) == 4
      if Tier2(i,1) == 1
          P1(3,i,:)= (AVGCOSTE(1,1)-AVGCOSTE(3,1)*Sub4)*0.66*Z;
          P2(3,i,:)= (AVGCOSTE(2,1)-AVGCOSTE(3,1)*Sub4)*0.66*Z;
          P3(3,i,:)= AVGCOSTE(3,1)*0.66*(1-Sub4)*Z;          
      end
      if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)
          P1(3,i,:)= (AVGCOSTS(1,1)-AVGCOSTS(3,1)*Sub4)*0.66*Z;
          P2(3,i,:)= (AVGCOSTS(2,1)-AVGCOSTS(3,1)*Sub4)*0.66*Z;
          P3(3,i,:)= AVGCOSTS(3,1)*0.66*(1-Sub4)*Z;          
      end
      if Tier2(i,1) == 6
          P1(3,i,:)= (AVGCOSTC(1,1)-AVGCOSTC(3,1)*Sub4)*0.66*Z1;
          P2(3,i,:)= (AVGCOSTC(2,1)-AVGCOSTC(3,1)*Sub4)*0.66*Z1;
          P3(3,i,:)= AVGCOSTC(3,1)*0.66*(1-Sub4)*Z1;          
      end      
      if Tier2(i,1) == 8
          P1(3,i,:)= (AVGCOSTF(1,1)-AVGCOSTF(3,1)*Sub4)*0.66*Z;
          P2(3,i,:)= (AVGCOSTF(2,1)-AVGCOSTF(3,1)*Sub4)*0.66*Z;
          P3(3,i,:)= AVGCOSTF(3,1)*0.66*(1-Sub4)*Z;          
      end
   end   
   
   if Inc2(i,1) == 5
      if Tier2(i,1) == 1
          P1(3,i,:)= (AVGCOSTE(1,1)-AVGCOSTE(3,1)*Sub5)*0.61*Z;
          P2(3,i,:)= (AVGCOSTE(2,1)-AVGCOSTE(3,1)*Sub5)*0.61*Z;
          P3(3,i,:)= AVGCOSTE(3,1)*0.61*(1-Sub5)*Z;          
      end
      if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)
          P1(3,i,:)= (AVGCOSTS(1,1)-AVGCOSTS(3,1)*Sub5)*0.61*Z;
          P2(3,i,:)= (AVGCOSTS(2,1)-AVGCOSTS(3,1)*Sub5)*0.61*Z;
          P3(3,i,:)= AVGCOSTS(3,1)*0.61*(1-Sub5)*Z;          
      end
      if Tier2(i,1) == 6
          P1(3,i,:)= (AVGCOSTC(1,1)-AVGCOSTC(3,1)*Sub5)*0.61*Z1;
          P2(3,i,:)= (AVGCOSTC(2,1)-AVGCOSTC(3,1)*Sub5)*0.61*Z1;
          P3(3,i,:)= AVGCOSTC(3,1)*0.61*(1-Sub5)*Z1;          
      end      
      if Tier2(i,1) == 8
          P1(3,i,:)= (AVGCOSTF(1,1)-AVGCOSTF(3,1)*Sub5)*0.61*Z;
          P2(3,i,:)= (AVGCOSTF(2,1)-AVGCOSTF(3,1)*Sub5)*0.61*Z;
          P3(3,i,:)= AVGCOSTF(3,1)*0.61*(1-Sub5)*Z;          
      end
   end    
 
end

P1(3,:,:) = max(reshape(P1(3,:,:),nIs,Sim),150);
P2(3,:,:) = max(reshape(P2(3,:,:),nIs,Sim),150);
P3(3,:,:) = max(reshape(P3(3,:,:),nIs,Sim),150);

    P3E = zeros(3,5);
    P3S = zeros(3,5);
    P3C = zeros(3,5);
    P3F = zeros(3,5);
    
    P3E(1,1) = (AVGCOSTE(1,1)-AVGCOSTE(3,1)*Sub1)*Z;
    P3E(2,1) = (AVGCOSTE(2,1)-AVGCOSTE(3,1)*Sub1)*Z;
    P3E(3,1) =  AVGCOSTE(3,1)*(1-Sub1)*Z;          
      
    P3S(1,1) = (AVGCOSTS(1,1)-AVGCOSTS(3,1)*Sub1)*Z;
    P3S(2,1) = (AVGCOSTS(2,1)-AVGCOSTS(3,1)*Sub1)*Z;
    P3S(3,1) =  AVGCOSTS(3,1)*(1-Sub1)*Z;          
   
    P3C(1,1) = (AVGCOSTC(1,1)-AVGCOSTC(3,1)*Sub1)*Z1;
    P3C(2,1) = (AVGCOSTC(2,1)-AVGCOSTC(3,1)*Sub1)*Z1;
    P3C(3,1) =  AVGCOSTC(3,1)*(1-Sub1)*Z1;          
    
    P3F(1,1) = (AVGCOSTF(1,1)-AVGCOSTF(3,1)*Sub1)*Z;
    P3F(2,1) = (AVGCOSTF(2,1)-AVGCOSTF(3,1)*Sub1)*Z;
    P3F(3,1) =  AVGCOSTF(3,1)*(1-Sub1)*Z;          
  
    P3E(1,2) = (AVGCOSTE(1,1)-AVGCOSTE(3,1)*Sub2)*Z;
    P3E(2,2) = (AVGCOSTE(2,1)-AVGCOSTE(3,1)*Sub2)*Z;
    P3E(3,2) =  AVGCOSTE(3,1)*(1-Sub2)*Z;          
      
    P3S(1,2) = (AVGCOSTS(1,1)-AVGCOSTS(3,1)*Sub2)*Z;
    P3S(2,2) = (AVGCOSTS(2,1)-AVGCOSTS(3,1)*Sub2)*Z;
    P3S(3,2) = AVGCOSTS(3,1)*(1-Sub2)*Z;           
   
    P3C(1,2) = (AVGCOSTC(1,1)-AVGCOSTC(3,1)*Sub2)*Z1;
    P3C(2,2) = (AVGCOSTC(2,1)-AVGCOSTC(3,1)*Sub2)*Z1;
    P3C(3,2) =  AVGCOSTC(3,1)*(1-Sub2)*Z1;         
    
    P3F(1,2) = (AVGCOSTF(1,1)-AVGCOSTF(3,1)*Sub2)*Z;
    P3F(2,2) = (AVGCOSTF(2,1)-AVGCOSTF(3,1)*Sub2)*Z;
    P3F(3,2) =  AVGCOSTF(3,1)*(1-Sub2)*Z;      

    P3E(1,3) = (AVGCOSTE(1,1)-AVGCOSTE(3,1)*Sub3)*Z;
    P3E(2,3) = (AVGCOSTE(2,1)-AVGCOSTE(3,1)*Sub3)*Z;
    P3E(3,3) =  AVGCOSTE(3,1)*(1-Sub3)*Z;         
      
    P3S(1,3) = (AVGCOSTS(1,1)-AVGCOSTS(3,1)*Sub3)*Z;
    P3S(2,3) = (AVGCOSTS(2,1)-AVGCOSTS(3,1)*Sub3)*Z;
    P3S(3,3) =  AVGCOSTS(3,1)*(1-Sub3)*Z;          
   
    P3C(1,3) = (AVGCOSTC(1,1)-AVGCOSTC(3,1)*Sub3)*Z1;
    P3C(2,3) = (AVGCOSTC(2,1)-AVGCOSTC(3,1)*Sub3)*Z1;
    P3C(3,3) = AVGCOSTC(3,1)*(1-Sub3)*Z1;        
    
    P3F(1,3) = (AVGCOSTF(1,1)-AVGCOSTF(3,1)*Sub3)*Z;
    P3F(2,3) = (AVGCOSTF(2,1)-AVGCOSTF(3,1)*Sub3)*Z;
    P3F(3,3) = AVGCOSTF(3,1)*(1-Sub3)*Z; 
            
    P3E(1,4) = (AVGCOSTE(1,1)-AVGCOSTE(3,1)*Sub4)*Z;
    P3E(2,4) = (AVGCOSTE(2,1)-AVGCOSTE(3,1)*Sub4)*Z;
    P3E(3,4) = AVGCOSTE(3,1)*(1-Sub4)*Z;         
      
    P3S(1,4) = (AVGCOSTS(1,1)-AVGCOSTS(3,1)*Sub4)*Z;
    P3S(2,4) = (AVGCOSTS(2,1)-AVGCOSTS(3,1)*Sub4)*Z;
    P3S(3,4) =  AVGCOSTS(3,1)*(1-Sub4)*Z;           
   
    P3C(1,4) = (AVGCOSTC(1,1)-AVGCOSTC(3,1)*Sub4)*Z1;
    P3C(2,4) = (AVGCOSTC(2,1)-AVGCOSTC(3,1)*Sub4)*Z1;
    P3C(3,4) =  AVGCOSTC(3,1)*(1-Sub4)*Z1;        
    
    P3F(1,4) = (AVGCOSTF(1,1)-AVGCOSTF(3,1)*Sub4)*Z;
    P3F(2,4) = (AVGCOSTF(2,1)-AVGCOSTF(3,1)*Sub4)*Z;
    P3F(3,4) = AVGCOSTF(3,1)*(1-Sub4)*Z;

    P3E(1,5) = (AVGCOSTE(1,1)-AVGCOSTE(3,1)*Sub5)*Z;
    P3E(2,5) = (AVGCOSTE(2,1)-AVGCOSTE(3,1)*Sub5)*Z;
    P3E(3,5) =  AVGCOSTE(3,1)*(1-Sub5)*Z;        
      
    P3S(1,5) = (AVGCOSTS(1,1)-AVGCOSTS(3,1)*Sub5)*Z;
    P3S(2,5) = (AVGCOSTS(2,1)-AVGCOSTS(3,1)*Sub5)*Z;
    P3S(3,5) = AVGCOSTS(3,1)*(1-Sub5)*Z;             
   
    P3C(1,5) = (AVGCOSTC(1,1)-AVGCOSTC(3,1)*Sub5)*Z1;
    P3C(2,5) = (AVGCOSTC(2,1)-AVGCOSTC(3,1)*Sub5)*Z1;
    P3C(3,5) =  AVGCOSTC(3,1)*(1-Sub5)*Z1;       
    
    P3F(1,5) = (AVGCOSTF(1,1)-AVGCOSTF(3,1)*Sub5)*Z;
    P3F(2,5) = (AVGCOSTF(2,1)-AVGCOSTF(3,1)*Sub5)*Z;
    P3F(3,5) = AVGCOSTF(3,1)*(1-Sub5)*Z;

    
P3F = max(P3F,150);
P3S = max(P3S,150);
P3C = max(P3C,150);
P3E = max(P3E,150);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Compute New Year 3 Choices  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1:nIs
    for s = 1:Sim
       if CC1(i,s)==1
            choice13(i,s) = 1;
       end
       if CC1(i,s)==2
            choice23(i,s) = 1;
       end
       if CC1(i,s)==3
            choice33(i,s) = 1;
       end
    end
end

for k = 1:K
    SwC13(:,k,:) = SCC3.*choice13;
    SwC23(:,k,:) = SCC3.*choice23;
    SwC33(:,k,:) = SCC3.*choice33;
end

%%%% New expected utility baselines calculated because of new premiums and inertia 

for i = 1:nIs
    for k = 1:K
        EUvector13(i,k,:) = max(0.1,(W - reshape(PPO250OP3(i,k,:),1,Sim) - reshape(P1(3,i,:),1,Sim)));
        EUvector23(i,k,:) = max(0.1,(W - reshape(PPO500OP3(i,k,:),1,Sim) - reshape(P2(3,i,:),1,Sim)));
        EUvector33(i,k,:) = max(0.1,(W - reshape(PPO1200OP3(i,k,:),1,Sim) - reshape(P3(3,i,:),1,Sim)));
    end
end

EUvector13 = exp(-((EUvector13 +SwC13).*max(RA,0.0000001)));
EUvector23 = exp(-((EUvector23 + SwC23 + CHTC + Epsilon5003).*max(RA,0.0000001)));
EUvector33 = exp(-((EUvector33 + CDHP + CHTC + SwC33 + Epsilon12003).*max(RA, 0.0000001)));

for n = 1:Sim
    EUvectorBx(:,1,n) = (1/K)*sum(reshape(EUvector13(:,:,n),nIs,K),2);
    EUvectorBx(:,2,n) = (1/K)*sum(reshape(EUvector23(:,:,n),nIs,K),2);
    EUvectorBx(:,3,n) = (1/K)*sum(reshape(EUvector33(:,:,n),nIs,K),2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Year 3 Choices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CC2 = zeros(nIs,Sim);
MktShare2 = zeros(nPlans,1);

TOTALCOSTE2 = zeros(3,1);
TOTALNUME2 = zeros(3,1);
TOTALNUME2i = zeros(3,nIs);

TOTALCOSTS2 = zeros(3,1);
TOTALNUMS2 = zeros(3,1);
TOTALNUMS2i = zeros(3,nIs);

TOTALCOSTC2 = zeros(3,1);
TOTALNUMC2 = zeros(3,1);
TOTALNUMC2i = zeros(3,nIs);

TOTALCOSTF2 = zeros(3,1);
TOTALNUMF2 = zeros(3,1);
TOTALNUMF2i = zeros(3,nIs);

for i = 1:nIs 

    PT3 = ones(nPlans,Sim);
    PT3(1:nPlans-1,:) =  (reshape(EUvectorBx(i,3,:),1,Sim)'*ones(1,nPlans-1))'./reshape(EUvectorBx(i,[1:nPlans-1],:),nPlans-1,Sim);
    PTAgg3 = zeros(1,Sim);
    PTAgg3 = sum(PT3,1);
    
    PT3(1,:) = PT3(1,:)./PTAgg3;
    PT3(2,:) = PT3(2,:)./PTAgg3;
    PT3(3,:) = PT3(3,:)./PTAgg3;
    
    PT3 = ((PT3'.^3)./(sum((PT3'.^3),2)*ones(1,nPlans)))';
    
    UNI = unifrnd(0,1,1,Sim);
    
    for j = 1:Sim
        if PT3(1,j) >= UNI(1,j)
           CC2(i,j)=1;
        elseif PT3(1,j)+PT3(2,j) >= UNI(1,j)
           CC2(i,j)=2;
        elseif (1-PT3(3,j)) < UNI(1,j)
           CC2(i,j)=3;
        end
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nIs
    for j = 1:Sim
      if Tier2(i,1) == 1  
        if CC2(i,j) == 1 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUME2i(1,i) = TOTALNUME2i(1,i)+1; 
        end
        if CC2(i,j) == 2 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUME2i(2,i) = TOTALNUME2i(2,i)+1; 
        end
        if CC2(i,j) == 3 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUME2i(3,i) = TOTALNUME2i(3,i)+1; 
        end
      end
      
    if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)  
        if CC2(i,j) == 1 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMS2i(1,i) = TOTALNUMS2i(1,i)+1; 
        end
        if CC2(i,j) == 2 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMS2i(2,i) = TOTALNUMS2i(2,i)+1; 
        end
        if CC2(i,j) == 3 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMS2i(3,i) = TOTALNUMS2i(3,i)+1; 
        end
    end    

    if Tier2(i,1) == 6  
        if CC2(i,j) == 1 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMC2i(1,i) = TOTALNUMC2i(1,i)+1; 
        end
        if CC2(i,j) == 2 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMC2i(2,i) = TOTALNUMC2i(2,i)+1; 
        end
        if CC2(i,j) == 3 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMC2i(3,i) = TOTALNUMC2i(3,i)+1; 
        end
    end   
    
   if Tier2(i,1) == 8  
        if CC2(i,j) == 1 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMF2i(1,i) = TOTALNUMF2i(1,i)+1; 
        end
        if CC2(i,j) == 2 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMF2i(2,i) = TOTALNUMF2i(2,i)+1; 
        end
        if CC2(i,j) == 3 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMF2i(3,i) = TOTALNUMF2i(3,i)+1; 
        end
    end  
    end
end

for i = 1:nIs 
    if Tier2(i,1)==1
    TOTALNUME2 = TOTALNUME2 + TOTALNUME2i(:,i)*(Sim/sum(TOTALNUME2i(:,i)));    
    end
    if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)  
    TOTALNUMS2 = TOTALNUMS2 + TOTALNUMS2i(:,i)*(Sim/sum(TOTALNUMS2i(:,i)));
    end
    if Tier2(i,1)==6
    TOTALNUMC2 = TOTALNUMC2 + TOTALNUMC2i(:,i)*(Sim/sum(TOTALNUMC2i(:,i)));
    end
    if Tier2(i,1)==8
    TOTALNUMF2 = TOTALNUMF2 + TOTALNUMF2i(:,i)*(Sim/sum(TOTALNUMF2i(:,i)));
    end 
end

for i = 1:nIs 
    if Tier2(i,1) == 1
    TOTALCOSTE2(1,1) = TOTALCOSTE2(1,1) + (1/K)*sum(reshape(PlanPaid3(i,:,1),1,K),2)*TOTALNUME2i(1,i)*(Sim/sum(TOTALNUME2i(:,i)));
    TOTALCOSTE2(2,1) = TOTALCOSTE2(2,1) + (1/K)*sum(reshape(PlanPaid3(i,:,2),1,K),2)*TOTALNUME2i(2,i)*(Sim/sum(TOTALNUME2i(:,i)));
    TOTALCOSTE2(3,1) = TOTALCOSTE2(3,1) + (1/K)*sum(reshape(PlanPaid3(i,:,3),1,K),2)*TOTALNUME2i(3,i)*(Sim/sum(TOTALNUME2i(:,i)));
    end
   
    if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)  
    TOTALCOSTS2(1,1) = TOTALCOSTS2(1,1) + (1/K)*sum(reshape(PlanPaid3(i,:,1),1,K),2)*TOTALNUMS2i(1,i)*(Sim/sum(TOTALNUMS2i(:,i)));
    TOTALCOSTS2(2,1) = TOTALCOSTS2(2,1) + (1/K)*sum(reshape(PlanPaid3(i,:,2),1,K),2)*TOTALNUMS2i(2,i)*(Sim/sum(TOTALNUMS2i(:,i)));
    TOTALCOSTS2(3,1) = TOTALCOSTS2(3,1) + (1/K)*sum(reshape(PlanPaid3(i,:,3),1,K),2)*TOTALNUMS2i(3,i)*(Sim/sum(TOTALNUMS2i(:,i)));
    end
    
    if Tier2(i,1)==6
    TOTALCOSTC2(1,1) = TOTALCOSTC2(1,1) + (1/K)*sum(reshape(PlanPaid3(i,:,1),1,K),2)*TOTALNUMC2i(1,i)*(Sim/sum(TOTALNUMC2i(:,i)));
    TOTALCOSTC2(2,1) = TOTALCOSTC2(2,1) + (1/K)*sum(reshape(PlanPaid3(i,:,2),1,K),2)*TOTALNUMC2i(2,i)*(Sim/sum(TOTALNUMC2i(:,i)));
    TOTALCOSTC2(3,1) = TOTALCOSTC2(3,1) + (1/K)*sum(reshape(PlanPaid3(i,:,3),1,K),2)*TOTALNUMC2i(3,i)*(Sim/sum(TOTALNUMC2i(:,i)));
    end    
       
   if Tier2(i,1)==8
   TOTALCOSTF2(1,1) = TOTALCOSTF2(1,1) + (1/K)*sum(reshape(PlanPaid3(i,:,1),1,K),2)*TOTALNUMF2i(1,i)*(Sim/sum(TOTALNUMF2i(:,i)));
   TOTALCOSTF2(2,1) = TOTALCOSTF2(2,1) + (1/K)*sum(reshape(PlanPaid3(i,:,2),1,K),2)*TOTALNUMF2i(2,i)*(Sim/sum(TOTALNUMF2i(:,i)));
   TOTALCOSTF2(3,1) = TOTALCOSTF2(3,1) + (1/K)*sum(reshape(PlanPaid3(i,:,3),1,K),2)*TOTALNUMF2i(3,i)*(Sim/sum(TOTALNUMF2i(:,i)));
   end
end

AVGCOSTE2 = zeros(3,1);
AVGCOSTS2 = zeros(3,1);
AVGCOSTC2 = zeros(3,1);
AVGCOSTF2 = zeros(3,1);

AVGCOSTE2(1,1) = TOTALCOSTE2(1,1)/TOTALNUME2(1,1);
AVGCOSTE2(2,1) = TOTALCOSTE2(2,1)/TOTALNUME2(2,1);
AVGCOSTE2(3,1) = TOTALCOSTE2(3,1)/TOTALNUME2(3,1);

AVGCOSTS2(1,1) = TOTALCOSTS2(1,1)/TOTALNUMS2(1,1);
AVGCOSTS2(2,1) = TOTALCOSTS2(2,1)/TOTALNUMS2(2,1);
AVGCOSTS2(3,1) = TOTALCOSTS2(3,1)/TOTALNUMS2(3,1);

AVGCOSTC2(1,1) = TOTALCOSTC2(1,1)/TOTALNUMC2(1,1);
AVGCOSTC2(2,1) = TOTALCOSTC2(2,1)/TOTALNUMC2(2,1);
AVGCOSTC2(3,1) = TOTALCOSTC2(3,1)/TOTALNUMC2(3,1);

AVGCOSTF2(1,1) = TOTALCOSTF2(1,1)/TOTALNUMF2(1,1);
AVGCOSTF2(2,1) = TOTALCOSTF2(2,1)/TOTALNUMF2(2,1);
AVGCOSTF2(3,1) = TOTALCOSTF2(3,1)/TOTALNUMF2(3,1);

MktShare2 = zeros(3,1);

MktShare2(1,1) = (TOTALNUMF2(1,1) + TOTALNUMC2(1,1) + TOTALNUMS2(1,1) + TOTALNUME2(1,1))/Sim;
MktShare2(2,1) = (TOTALNUMF2(2,1) + TOTALNUMC2(2,1) + TOTALNUMS2(2,1) + TOTALNUME2(2,1))/Sim;
MktShare2(3,1) = (TOTALNUMF2(3,1) + TOTALNUMC2(3,1) + TOTALNUMS2(3,1) + TOTALNUME2(3,1))/Sim;

for i = 1:nIs

   if Inc2(i,1) == 1
      if Tier2(i,1) == 1
          P1(4,i,:)= (AVGCOSTE2(1,1)-AVGCOSTE2(3,1)*Sub1)*0.82*Z;
          P2(4,i,:)= (AVGCOSTE2(2,1)-AVGCOSTE2(3,1)*Sub1)*0.82*Z;
          P3(4,i,:)= AVGCOSTE2(3,1)*0.82*(1-Sub1)*Z;          
      end
      if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)
          P1(4,i,:)= (AVGCOSTS2(1,1)-AVGCOSTS2(3,1)*Sub1)*0.82*Z;
          P2(4,i,:)= (AVGCOSTS2(2,1)-AVGCOSTS2(3,1)*Sub1)*0.82*Z;
          P3(4,i,:)= AVGCOSTS2(3,1)*0.82*(1-Sub1)*Z;          
      end
      if Tier2(i,1) == 6
          P1(4,i,:)= (AVGCOSTC2(1,1)-AVGCOSTC2(3,1)*Sub1)*0.82*Z1;
          P2(4,i,:)= (AVGCOSTC2(2,1)-AVGCOSTC2(3,1)*Sub1)*0.82*Z1;
          P3(4,i,:)= AVGCOSTC2(3,1)*0.82*(1-Sub1)*Z1;          
      end      
      if Tier2(i,1) == 8
          P1(4,i,:)= (AVGCOSTF2(1,1)-AVGCOSTF2(3,1)*Sub1)*0.82*Z;
          P2(4,i,:)= (AVGCOSTF2(2,1)-AVGCOSTF2(3,1)*Sub1)*0.82*Z;
          P3(4,i,:)= AVGCOSTF2(3,1)*0.82*(1-Sub1)*Z;          
      end
   end
   
   if Inc2(i,1) == 2
      if Tier2(i,1) == 1
          P1(4,i,:)= (AVGCOSTE2(1,1)-AVGCOSTE2(3,1)*Sub2)*0.73*Z;
          P2(4,i,:)= (AVGCOSTE2(2,1)-AVGCOSTE2(3,1)*Sub2)*0.73*Z;
          P3(4,i,:)= AVGCOSTE2(3,1)*0.73*(1-Sub2)*Z;          
      end
      if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)
          P1(4,i,:)= (AVGCOSTS2(1,1)-AVGCOSTS2(3,1)*Sub2)*0.80*Z;
          P2(4,i,:)= (AVGCOSTS2(2,1)-AVGCOSTS2(3,1)*Sub2)*0.80*Z;
          P3(4,i,:)= AVGCOSTS2(3,1)*0.80*(1-Sub2)*Z;          
      end
      if Tier2(i,1) == 6
          P1(4,i,:)= (AVGCOSTC2(1,1)-AVGCOSTC2(3,1)*Sub2)*0.80*Z1;
          P2(4,i,:)= (AVGCOSTC2(2,1)-AVGCOSTC2(3,1)*Sub2)*0.80*Z1;
          P3(4,i,:)= AVGCOSTC2(3,1)*0.80*(1-Sub2)*Z1;          
      end      
      if Tier2(i,1) == 8
          P1(4,i,:)= (AVGCOSTF2(1,1)-AVGCOSTF2(3,1)*Sub2)*0.80*Z;
          P2(4,i,:)= (AVGCOSTF2(2,1)-AVGCOSTF2(3,1)*Sub2)*0.80*Z;
          P3(4,i,:)= AVGCOSTF2(3,1)*0.80*(1-Sub2)*Z;          
      end
   end
   if Inc2(i,1) == 3
      if Tier2(i,1) == 1
          P1(4,i,:)= (AVGCOSTE2(1,1)-AVGCOSTE2(3,1)*Sub3)*0.69*Z;
          P2(4,i,:)= (AVGCOSTE2(2,1)-AVGCOSTE2(3,1)*Sub3)*0.69*Z;
          P3(4,i,:)= AVGCOSTE2(3,1)*0.69*(1-Sub3)*Z;          
      end
      if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)
          P1(4,i,:)= (AVGCOSTS2(1,1)-AVGCOSTS2(3,1)*Sub3)*0.69*Z;
          P2(4,i,:)= (AVGCOSTS2(2,1)-AVGCOSTS2(3,1)*Sub3)*0.69*Z;
          P3(4,i,:)= AVGCOSTS2(3,1)*0.69*(1-Sub3)*Z;          
      end
      if Tier2(i,1) == 6
          P1(4,i,:)= (AVGCOSTC2(1,1)-AVGCOSTC2(3,1)*Sub3)*0.69*Z1;
          P2(4,i,:)= (AVGCOSTC2(2,1)-AVGCOSTC2(3,1)*Sub3)*0.69*Z1;
          P3(4,i,:)= AVGCOSTC2(3,1)*0.69*(1-Sub3)*Z1;          
      end      
      if Tier2(i,1) == 8
          P1(4,i,:)= (AVGCOSTF2(1,1)-AVGCOSTF2(3,1)*Sub3)*0.69*Z;
          P2(4,i,:)= (AVGCOSTF2(2,1)-AVGCOSTF2(3,1)*Sub3)*0.69*Z;
          P3(4,i,:)= AVGCOSTF2(3,1)*0.69*(1-Sub3)*Z;          
      end
   end
   
   if Inc2(i,1) == 4
      if Tier2(i,1) == 1
          P1(4,i,:)= (AVGCOSTE2(1,1)-AVGCOSTE2(3,1)*Sub4)*0.66*Z;
          P2(4,i,:)= (AVGCOSTE2(2,1)-AVGCOSTE2(3,1)*Sub4)*0.66*Z;
          P3(4,i,:)= AVGCOSTE2(3,1)*0.66*(1-Sub4)*Z;          
      end
      if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)
          P1(4,i,:)= (AVGCOSTS2(1,1)-AVGCOSTS2(3,1)*Sub4)*0.66*Z;
          P2(4,i,:)= (AVGCOSTS2(2,1)-AVGCOSTS2(3,1)*Sub4)*0.66*Z;
          P3(4,i,:)= AVGCOSTS2(3,1)*0.66*(1-Sub4)*Z;          
      end
      if Tier2(i,1) == 6
          P1(4,i,:)= (AVGCOSTC2(1,1)-AVGCOSTC2(3,1)*Sub4)*0.66*Z1;
          P2(4,i,:)= (AVGCOSTC2(2,1)-AVGCOSTC2(3,1)*Sub4)*0.66*Z1;
          P3(4,i,:)= AVGCOSTC2(3,1)*0.66*(1-Sub4)*Z1;          
      end      
      if Tier2(i,1) == 8
          P1(4,i,:)= (AVGCOSTF2(1,1)-AVGCOSTF2(3,1)*Sub4)*0.66*Z;
          P2(4,i,:)= (AVGCOSTF2(2,1)-AVGCOSTF2(3,1)*Sub4)*0.66*Z;
          P3(4,i,:)= AVGCOSTF2(3,1)*0.66*(1-Sub4)*Z;          
      end
   end   
   
   if Inc2(i,1) == 5
      if Tier2(i,1) == 1
          P1(4,i,:)= (AVGCOSTE2(1,1)-AVGCOSTE2(3,1)*Sub5)*0.61*Z;
          P2(4,i,:)= (AVGCOSTE2(2,1)-AVGCOSTE2(3,1)*Sub5)*0.61*Z;
          P3(4,i,:)= AVGCOSTE2(3,1)*0.61*(1-Sub5)*Z;          
      end
      if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)
          P1(4,i,:)= (AVGCOSTS2(1,1)-AVGCOSTS2(3,1)*Sub5)*0.61*Z;
          P2(4,i,:)= (AVGCOSTS2(2,1)-AVGCOSTS2(3,1)*Sub5)*0.61*Z;
          P3(4,i,:)= AVGCOSTS2(3,1)*0.61*(1-Sub5)*Z;          
      end
      if Tier2(i,1) == 6
          P1(4,i,:)= (AVGCOSTC2(1,1)-AVGCOSTC2(3,1)*Sub5)*0.61*Z1;
          P2(4,i,:)= (AVGCOSTC2(2,1)-AVGCOSTC2(3,1)*Sub5)*0.61*Z1;
          P3(4,i,:)= AVGCOSTC2(3,1)*0.61*(1-Sub5)*Z1;          
      end      
      if Tier2(i,1) == 8
          P1(4,i,:)= (AVGCOSTF2(1,1)-AVGCOSTF2(3,1)*Sub5)*0.61*Z;
          P2(4,i,:)= (AVGCOSTF2(2,1)-AVGCOSTF2(3,1)*Sub5)*0.61*Z;
          P3(4,i,:)= AVGCOSTF2(3,1)*0.61*(1-Sub5)*Z;          
      end
   end    
 
end

P1(4,:,:) = max(reshape(P1(4,:,:),nIs,Sim),150);
P2(4,:,:) = max(reshape(P2(4,:,:),nIs,Sim),150);
P3(4,:,:) = max(reshape(P3(4,:,:),nIs,Sim),150);

    P3E3 = zeros(3,5);
    P3S3 = zeros(3,5);
    P3C3 = zeros(3,5);
    P3F3 = zeros(3,5);
    
    P3E3(1,1) = (AVGCOSTE2(1,1)-AVGCOSTE2(3,1)*Sub1)*Z;
    P3E3(2,1) = (AVGCOSTE2(2,1)-AVGCOSTE2(3,1)*Sub1)*Z;
    P3E3(3,1) =  AVGCOSTE2(3,1)*(1-Sub1)*Z;          
      
    P3S3(1,1) = (AVGCOSTS2(1,1)-AVGCOSTS2(3,1)*Sub1)*Z;
    P3S3(2,1) = (AVGCOSTS2(2,1)-AVGCOSTS2(3,1)*Sub1)*Z;
    P3S3(3,1) =  AVGCOSTS2(3,1)*(1-Sub1)*Z;          
   
    P3C3(1,1) = (AVGCOSTC2(1,1)-AVGCOSTC2(3,1)*Sub1)*Z1;
    P3C3(2,1) = (AVGCOSTC2(2,1)-AVGCOSTC2(3,1)*Sub1)*Z1;
    P3C3(3,1) =  AVGCOSTC2(3,1)*(1-Sub1)*Z1;          
    
    P3F3(1,1) = (AVGCOSTF2(1,1)-AVGCOSTF2(3,1)*Sub1)*Z;
    P3F3(2,1) = (AVGCOSTF2(2,1)-AVGCOSTF2(3,1)*Sub1)*Z;
    P3F3(3,1) =  AVGCOSTF2(3,1)*(1-Sub1)*Z;          
  
    P3E3(1,2) = (AVGCOSTE2(1,1)-AVGCOSTE2(3,1)*Sub2)*Z;
    P3E3(2,2) = (AVGCOSTE2(2,1)-AVGCOSTE2(3,1)*Sub2)*Z;
    P3E3(3,2) =  AVGCOSTE2(3,1)*(1-Sub2)*Z;          
      
    P3S3(1,2) = (AVGCOSTS2(1,1)-AVGCOSTS2(3,1)*Sub2)*Z;
    P3S3(2,2) = (AVGCOSTS2(2,1)-AVGCOSTS2(3,1)*Sub2)*Z;
    P3S3(3,2) = AVGCOSTS2(3,1)*(1-Sub2)*Z;           
   
    P3C3(1,2) = (AVGCOSTC2(1,1)-AVGCOSTC2(3,1)*Sub2)*Z1;
    P3C3(2,2) = (AVGCOSTC2(2,1)-AVGCOSTC2(3,1)*Sub2)*Z1;
    P3C3(3,2) =  AVGCOSTC2(3,1)*(1-Sub2)*Z1;         
    
    P3F3(1,2) = (AVGCOSTF2(1,1)-AVGCOSTF2(3,1)*Sub2)*Z;
    P3F3(2,2) = (AVGCOSTF2(2,1)-AVGCOSTF2(3,1)*Sub2)*Z;
    P3F3(3,2) =  AVGCOSTF2(3,1)*(1-Sub2)*Z;      

    P3E3(1,3) = (AVGCOSTE2(1,1)-AVGCOSTE2(3,1)*Sub3)*Z;
    P3E3(2,3) = (AVGCOSTE2(2,1)-AVGCOSTE2(3,1)*Sub3)*Z;
    P3E3(3,3) =  AVGCOSTE2(3,1)*(1-Sub3)*Z;         
      
    P3S3(1,3) = (AVGCOSTS2(1,1)-AVGCOSTS2(3,1)*Sub3)*Z;
    P3S3(2,3) = (AVGCOSTS2(2,1)-AVGCOSTS2(3,1)*Sub3)*Z;
    P3S3(3,3) =  AVGCOSTS2(3,1)*(1-Sub3)*Z;          
   
    P3C3(1,3) = (AVGCOSTC2(1,1)-AVGCOSTC2(3,1)*Sub3)*Z1;
    P3C3(2,3) = (AVGCOSTC2(2,1)-AVGCOSTC2(3,1)*Sub3)*Z1;
    P3C3(3,3) = AVGCOSTC2(3,1)*(1-Sub3)*Z1;        
    
    P3F3(1,3) = (AVGCOSTF2(1,1)-AVGCOSTF2(3,1)*Sub3)*Z;
    P3F3(2,3) = (AVGCOSTF2(2,1)-AVGCOSTF2(3,1)*Sub3)*Z;
    P3F3(3,3) = AVGCOSTF2(3,1)*(1-Sub3)*Z; 
            
    P3E3(1,4) = (AVGCOSTE2(1,1)-AVGCOSTE2(3,1)*Sub4)*Z;
    P3E3(2,4) = (AVGCOSTE2(2,1)-AVGCOSTE2(3,1)*Sub4)*Z;
    P3E3(3,4) = AVGCOSTE2(3,1)*(1-Sub4)*Z;         
      
    P3S3(1,4) = (AVGCOSTS2(1,1)-AVGCOSTS2(3,1)*Sub4)*Z;
    P3S3(2,4) = (AVGCOSTS2(2,1)-AVGCOSTS2(3,1)*Sub4)*Z;
    P3S3(3,4) =  AVGCOSTS2(3,1)*(1-Sub4)*Z;           
   
    P3C3(1,4) = (AVGCOSTC2(1,1)-AVGCOSTC2(3,1)*Sub4)*Z1;
    P3C3(2,4) = (AVGCOSTC2(2,1)-AVGCOSTC2(3,1)*Sub4)*Z1;
    P3C3(3,4) =  AVGCOSTC2(3,1)*(1-Sub4)*Z1;        
    
    P3F3(1,4) = (AVGCOSTF2(1,1)-AVGCOSTF2(3,1)*Sub4)*Z;
    P3F3(2,4) = (AVGCOSTF2(2,1)-AVGCOSTF2(3,1)*Sub4)*Z;
    P3F3(3,4) = AVGCOSTF2(3,1)*(1-Sub4)*Z;

    P3E3(1,5) = (AVGCOSTE2(1,1)-AVGCOSTE2(3,1)*Sub5)*Z;
    P3E3(2,5) = (AVGCOSTE2(2,1)-AVGCOSTE2(3,1)*Sub5)*Z;
    P3E3(3,5) =  AVGCOSTE2(3,1)*(1-Sub5)*Z;        
      
    P3S3(1,5) = (AVGCOSTS2(1,1)-AVGCOSTS2(3,1)*Sub5)*Z;
    P3S3(2,5) = (AVGCOSTS2(2,1)-AVGCOSTS2(3,1)*Sub5)*Z;
    P3S3(3,5) = AVGCOSTS2(3,1)*(1-Sub5)*Z;             
   
    P3C3(1,5) = (AVGCOSTC2(1,1)-AVGCOSTC2(3,1)*Sub5)*Z1;
    P3C3(2,5) = (AVGCOSTC2(2,1)-AVGCOSTC2(3,1)*Sub5)*Z1;
    P3C3(3,5) =  AVGCOSTC2(3,1)*(1-Sub5)*Z1;       
    
    P3F3(1,5) = (AVGCOSTF2(1,1)-AVGCOSTF2(3,1)*Sub5)*Z;
    P3F3(2,5) = (AVGCOSTF2(2,1)-AVGCOSTF2(3,1)*Sub5)*Z;
    P3F3(3,5) = AVGCOSTF2(3,1)*(1-Sub5)*Z;

P3F3 = max(P3F3,150);
P3S3 = max(P3S3,150);
P3C3 = max(P3C3,150);
P3E3 = max(P3E3,150);

%%%%%%%%%%%%%% Incumbent Plans for next year 

for i = 1:nIs
    for s = 1:Sim
       
       if CC2(i,s)==1
            choice14(i,s) = 1;
       end
       if CC2(i,s)==2
            choice24(i,s) = 1;
       end
       if CC2(i,s)==3
            choice34(i,s) = 1;
       end
    
    end
end

for k = 1:K
    SwC14(:,k,:) = SCC3.*choice14;
    SwC24(:,k,:) = SCC3.*choice24;
    SwC34(:,k,:) = SCC3.*choice34;
end

%%%% New expected utility baselines calculated because of new premiums and inertia 

for i = 1:nIs
    for k = 1:K
        EUvector14(i,k,:) = max(0.1,(W - reshape(PPO250OP3(i,k,:),1,Sim) - reshape(P1(4,i,:),1,Sim)));
        EUvector24(i,k,:) = max(0.1,(W - reshape(PPO500OP3(i,k,:),1,Sim) - reshape(P2(4,i,:),1,Sim)));
        EUvector34(i,k,:) = max(0.1,(W - reshape(PPO1200OP3(i,k,:),1,Sim) - reshape(P3(4,i,:),1,Sim)));
    end
end

EUvector14 = exp(-((EUvector14 +SwC14).*max(RA,0.0000001)));
EUvector24 = exp(-((EUvector24 + SwC24 + CHTC + Epsilon5004).*max(RA,0.0000001)));
EUvector34 = exp(-((EUvector34 + CDHP + CHTC + SwC34 + Epsilon12004).*max(RA, 0.0000001)));

%%%%%%%%%%%%%%%%%%% Now sum CARA VNM to get expected utility %%%%%%%%%%%%

for n = 1:Sim
    EUvectorC(:,1,n) = (1/K)*sum(reshape(EUvector14(:,:,n),nIs,K),2);
    EUvectorC(:,2,n) = (1/K)*sum(reshape(EUvector24(:,:,n),nIs,K),2);
    EUvectorC(:,3,n) = (1/K)*sum(reshape(EUvector34(:,:,n),nIs,K),2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Year 4 Choices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CC3 = zeros(nIs,Sim);
MktShare3 = zeros(nPlans,1);

TOTALCOSTE3 = zeros(3,1);
TOTALNUME3 = zeros(3,1);
TOTALNUME3i = zeros(3,nIs);

TOTALCOSTS3 = zeros(3,1);
TOTALNUMS3 = zeros(3,1);
TOTALNUMS3i = zeros(3,nIs);

TOTALCOSTC3 = zeros(3,1);
TOTALNUMC3 = zeros(3,1);
TOTALNUMC3i = zeros(3,nIs);

TOTALCOSTF3 = zeros(3,1);
TOTALNUMF3 = zeros(3,1);
TOTALNUMF3i = zeros(3,nIs);

for i = 1:nIs 
    
    PT4 = ones(nPlans,Sim);
    PT4(1:nPlans-1,:) =  (reshape(EUvectorC(i,3,:),1,Sim)'*ones(1,nPlans-1))'./reshape(EUvectorC(i,[1:nPlans-1],:),nPlans-1,Sim);
    PTAgg4 = zeros(1,Sim);
    PTAgg4 = sum(PT4,1);
    
    PT4(1,:) = PT4(1,:)./PTAgg4;
    PT4(2,:) = PT4(2,:)./PTAgg4;
    PT4(3,:) = PT4(3,:)./PTAgg4;

    PT4 = ((PT4'.^3)./(sum((PT4'.^3),2)*ones(1,nPlans)))';
    
    UNI = unifrnd(0,1,1,Sim);
    
    for j = 1:Sim
        if PT4(1,j) >= UNI(1,j)
           CC3(i,j)=1;
        elseif PT4(1,j)+PT4(2,j) >= UNI(1,j)
           CC3(i,j)=2;
        elseif (1-PT4(3,j)) < UNI(1,j)
           CC3(i,j)=3;
        end
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nIs
    for j = 1:Sim
      if Tier2(i,1) == 1  
        if CC3(i,j) == 1 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUME3i(1,i) = TOTALNUME3i(1,i)+1; 
        end
        if CC3(i,j) == 2 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUME3i(2,i) = TOTALNUME3i(2,i)+1; 
        end
        if CC3(i,j) == 3 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUME3i(3,i) = TOTALNUME3i(3,i)+1; 
        end
      end
      
    if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)  
        if CC3(i,j) == 1 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMS3i(1,i) = TOTALNUMS3i(1,i)+1; 
        end
        if CC3(i,j) == 2 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMS3i(2,i) = TOTALNUMS3i(2,i)+1; 
        end
        if CC3(i,j) == 3 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMS3i(3,i) = TOTALNUMS3i(3,i)+1; 
        end
    end    

    if Tier2(i,1) == 6  
        if CC3(i,j) == 1 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMC3i(1,i) = TOTALNUMC3i(1,i)+1; 
        end
        if CC3(i,j) == 2 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMC3i(2,i) = TOTALNUMC3i(2,i)+1; 
        end
        if CC3(i,j) == 3 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMC3i(3,i) = TOTALNUMC3i(3,i)+1; 
        end
    end   
    
   if Tier2(i,1) == 8  
        if CC3(i,j) == 1 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMF3i(1,i) = TOTALNUMF3i(1,i)+1; 
        end
        if CC3(i,j) == 2 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMF3i(2,i) = TOTALNUMF3i(2,i)+1; 
        end
        if CC3(i,j) == 3 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMF3i(3,i) = TOTALNUMF3i(3,i)+1; 
        end
    end  
    end
end

for i = 1:nIs
    if Tier2(i,1)==1
    TOTALNUME3 = TOTALNUME3 + TOTALNUME3i(:,i)*(Sim/sum(TOTALNUME3i(:,i)));    
    end
    if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)  
    TOTALNUMS3 = TOTALNUMS3 + TOTALNUMS3i(:,i)*(Sim/sum(TOTALNUMS3i(:,i)));
    end
    if Tier2(i,1)==6
    TOTALNUMC3 = TOTALNUMC3 + TOTALNUMC3i(:,i)*(Sim/sum(TOTALNUMC3i(:,i)));
    end
    if Tier2(i,1)==8
    TOTALNUMF3 = TOTALNUMF3 + TOTALNUMF3i(:,i)*(Sim/sum(TOTALNUMF3i(:,i)));
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nIs
      
    if Tier2(i,1) == 1
    TOTALCOSTE3(1,1) = TOTALCOSTE3(1,1) + (1/K)*sum(reshape(PlanPaid3(i,:,1),1,K),2)*TOTALNUME3i(1,i)*(Sim/sum(TOTALNUME3i(:,i)));
    TOTALCOSTE3(2,1) = TOTALCOSTE3(2,1) + (1/K)*sum(reshape(PlanPaid3(i,:,2),1,K),2)*TOTALNUME3i(2,i)*(Sim/sum(TOTALNUME3i(:,i)));
    TOTALCOSTE3(3,1) = TOTALCOSTE3(3,1) + (1/K)*sum(reshape(PlanPaid3(i,:,3),1,K),2)*TOTALNUME3i(3,i)*(Sim/sum(TOTALNUME3i(:,i)));
    end
   
    if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)  
    TOTALCOSTS3(1,1) = TOTALCOSTS3(1,1) + (1/K)*sum(reshape(PlanPaid3(i,:,1),1,K),2)*TOTALNUMS3i(1,i)*(Sim/sum(TOTALNUMS3i(:,i)));
    TOTALCOSTS3(2,1) = TOTALCOSTS3(2,1) + (1/K)*sum(reshape(PlanPaid3(i,:,2),1,K),2)*TOTALNUMS3i(2,i)*(Sim/sum(TOTALNUMS3i(:,i)));
    TOTALCOSTS3(3,1) = TOTALCOSTS3(3,1) + (1/K)*sum(reshape(PlanPaid3(i,:,3),1,K),2)*TOTALNUMS3i(3,i)*(Sim/sum(TOTALNUMS3i(:,i)));
    end
    
    if Tier2(i,1)==6
    TOTALCOSTC3(1,1) = TOTALCOSTC3(1,1) + (1/K)*sum(reshape(PlanPaid3(i,:,1),1,K),2)*TOTALNUMC3i(1,i)*(Sim/sum(TOTALNUMC3i(:,i)));
    TOTALCOSTC3(2,1) = TOTALCOSTC3(2,1) + (1/K)*sum(reshape(PlanPaid3(i,:,2),1,K),2)*TOTALNUMC3i(2,i)*(Sim/sum(TOTALNUMC3i(:,i)));
    TOTALCOSTC3(3,1) = TOTALCOSTC3(3,1) + (1/K)*sum(reshape(PlanPaid3(i,:,3),1,K),2)*TOTALNUMC3i(3,i)*(Sim/sum(TOTALNUMC3i(:,i)));
    end    
       
   if Tier2(i,1)==8
   TOTALCOSTF3(1,1) = TOTALCOSTF3(1,1) + (1/K)*sum(reshape(PlanPaid3(i,:,1),1,K),2)*TOTALNUMF3i(1,i)*(Sim/sum(TOTALNUMF3i(:,i)));
   TOTALCOSTF3(2,1) = TOTALCOSTF3(2,1) + (1/K)*sum(reshape(PlanPaid3(i,:,2),1,K),2)*TOTALNUMF3i(2,i)*(Sim/sum(TOTALNUMF3i(:,i)));
   TOTALCOSTF3(3,1) = TOTALCOSTF3(3,1) + (1/K)*sum(reshape(PlanPaid3(i,:,3),1,K),2)*TOTALNUMF3i(3,i)*(Sim/sum(TOTALNUMF3i(:,i)));
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AVGCOSTE3 = zeros(3,1);
AVGCOSTS3 = zeros(3,1);
AVGCOSTC3 = zeros(3,1);
AVGCOSTF3 = zeros(3,1);

%%%%%%%% Incorporate Adminsitrative Loadinig Factors and Risk Buffer Funds

AVGCOSTE3(1,1) = TOTALCOSTE3(1,1)/TOTALNUME3(1,1);
AVGCOSTE3(2,1) = TOTALCOSTE3(2,1)/TOTALNUME3(2,1);
AVGCOSTE3(3,1) = TOTALCOSTE3(3,1)/TOTALNUME3(3,1);

AVGCOSTS3(1,1) = TOTALCOSTS3(1,1)/TOTALNUMS3(1,1);
AVGCOSTS3(2,1) = TOTALCOSTS3(2,1)/TOTALNUMS3(2,1);
AVGCOSTS3(3,1) = TOTALCOSTS3(3,1)/TOTALNUMS3(3,1);

AVGCOSTC3(1,1) = TOTALCOSTC3(1,1)/TOTALNUMC3(1,1);
AVGCOSTC3(2,1) = TOTALCOSTC3(2,1)/TOTALNUMC3(2,1);
AVGCOSTC3(3,1) = TOTALCOSTC3(3,1)/TOTALNUMC3(3,1);

AVGCOSTF3(1,1) = TOTALCOSTF3(1,1)/TOTALNUMF3(1,1);
AVGCOSTF3(2,1) = TOTALCOSTF3(2,1)/TOTALNUMF3(2,1);
AVGCOSTF3(3,1) = TOTALCOSTF3(3,1)/TOTALNUMF3(3,1);

MktShare3 = zeros(3,1);

MktShare3(1,1) = (TOTALNUMF3(1,1) + TOTALNUMC3(1,1) + TOTALNUMS3(1,1) + TOTALNUME3(1,1))/Sim;
MktShare3(2,1) = (TOTALNUMF3(2,1) + TOTALNUMC3(2,1) + TOTALNUMS3(2,1) + TOTALNUME3(2,1))/Sim;
MktShare3(3,1) = (TOTALNUMF3(3,1) + TOTALNUMC3(3,1) + TOTALNUMS3(3,1) + TOTALNUME3(3,1))/Sim;

%%%%%%%%%%%% Prices for Next Year

for i = 1:nIs

   if Inc2(i,1) == 1
      if Tier2(i,1) == 1
          P1(5,i,:)= (AVGCOSTE3(1,1)-AVGCOSTE3(3,1)*Sub1)*0.82*Z;
          P2(5,i,:)= (AVGCOSTE3(2,1)-AVGCOSTE3(3,1)*Sub1)*0.82*Z;
          P3(5,i,:)= AVGCOSTE3(3,1)*0.82*(1-Sub1)*Z;          
      end
      if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)
          P1(5,i,:)= (AVGCOSTS3(1,1)-AVGCOSTS3(3,1)*Sub1)*0.82*Z;
          P2(5,i,:)= (AVGCOSTS3(2,1)-AVGCOSTS3(3,1)*Sub1)*0.82*Z;
          P3(5,i,:)= AVGCOSTS3(3,1)*0.82*(1-Sub1)*Z;          
      end
      if Tier2(i,1) == 6
          P1(5,i,:)= (AVGCOSTC3(1,1)-AVGCOSTC3(3,1)*Sub1)*0.82*Z1;
          P2(5,i,:)= (AVGCOSTC3(2,1)-AVGCOSTC3(3,1)*Sub1)*0.82*Z1;
          P3(5,i,:)= AVGCOSTC3(3,1)*0.82*(1-Sub1)*Z1;          
      end      
      if Tier2(i,1) == 8
          P1(5,i,:)= (AVGCOSTF3(1,1)-AVGCOSTF3(3,1)*Sub1)*0.82*Z;
          P2(5,i,:)= (AVGCOSTF3(2,1)-AVGCOSTF3(3,1)*Sub1)*0.82*Z;
          P3(5,i,:)= AVGCOSTF3(3,1)*0.82*(1-Sub1)*Z;          
      end
   end
   
   if Inc2(i,1) == 2
      if Tier2(i,1) == 1
          P1(5,i,:)= (AVGCOSTE3(1,1)-AVGCOSTE3(3,1)*Sub2)*0.73*Z;
          P2(5,i,:)= (AVGCOSTE3(2,1)-AVGCOSTE3(3,1)*Sub2)*0.73*Z;
          P3(5,i,:)= AVGCOSTE3(3,1)*0.73*(1-Sub2)*Z;          
      end
      if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)
          P1(5,i,:)= (AVGCOSTS3(1,1)-AVGCOSTS3(3,1)*Sub2)*0.80*Z;
          P2(5,i,:)= (AVGCOSTS3(2,1)-AVGCOSTS3(3,1)*Sub2)*0.80*Z;
          P3(5,i,:)= AVGCOSTS3(3,1)*0.80*(1-Sub2)*Z;          
      end
      if Tier2(i,1) == 6
          P1(5,i,:)= (AVGCOSTC3(1,1)-AVGCOSTC3(3,1)*Sub2)*0.80*Z1;
          P2(5,i,:)= (AVGCOSTC3(2,1)-AVGCOSTC3(3,1)*Sub2)*0.80*Z1;
          P3(5,i,:)= AVGCOSTC3(3,1)*0.80*(1-Sub2)*Z1;          
      end      
      if Tier2(i,1) == 8
          P1(5,i,:)= (AVGCOSTF3(1,1)-AVGCOSTF3(3,1)*Sub2)*0.80*Z;
          P2(5,i,:)= (AVGCOSTF3(2,1)-AVGCOSTF3(3,1)*Sub2)*0.80*Z;
          P3(5,i,:)= AVGCOSTF3(3,1)*0.80*(1-Sub2)*Z;          
      end
   end
   if Inc2(i,1) == 3
      if Tier2(i,1) == 1
          P1(5,i,:)= (AVGCOSTE3(1,1)-AVGCOSTE3(3,1)*Sub3)*0.69*Z;
          P2(5,i,:)= (AVGCOSTE3(2,1)-AVGCOSTE3(3,1)*Sub3)*0.69*Z;
          P3(5,i,:)= AVGCOSTE3(3,1)*0.69*(1-Sub3)*Z;          
      end
      if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)
          P1(5,i,:)= (AVGCOSTS3(1,1)-AVGCOSTS3(3,1)*Sub3)*0.69*Z;
          P2(5,i,:)= (AVGCOSTS3(2,1)-AVGCOSTS3(3,1)*Sub3)*0.69*Z;
          P3(5,i,:)= AVGCOSTS3(3,1)*0.69*(1-Sub3)*Z;          
      end
      if Tier2(i,1) == 6
          P1(5,i,:)= (AVGCOSTC3(1,1)-AVGCOSTC3(3,1)*Sub3)*0.69*Z1;
          P2(5,i,:)= (AVGCOSTC3(2,1)-AVGCOSTC3(3,1)*Sub3)*0.69*Z1;
          P3(5,i,:)= AVGCOSTC3(3,1)*0.69*(1-Sub3)*Z1;          
      end      
      if Tier2(i,1) == 8
          P1(5,i,:)= (AVGCOSTF3(1,1)-AVGCOSTF3(3,1)*Sub3)*0.69*Z;
          P2(5,i,:)= (AVGCOSTF3(2,1)-AVGCOSTF3(3,1)*Sub3)*0.69*Z;
          P3(5,i,:)= AVGCOSTF3(3,1)*0.69*(1-Sub3)*Z;          
      end
   end
   
   if Inc2(i,1) == 4
      if Tier2(i,1) == 1
          P1(5,i,:)= (AVGCOSTE3(1,1)-AVGCOSTE3(3,1)*Sub4)*0.66*Z;
          P2(5,i,:)= (AVGCOSTE3(2,1)-AVGCOSTE3(3,1)*Sub4)*0.66*Z;
          P3(5,i,:)= AVGCOSTE3(3,1)*0.66*(1-Sub4)*Z;          
      end
      if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)
          P1(5,i,:)= (AVGCOSTS3(1,1)-AVGCOSTS3(3,1)*Sub4)*0.66*Z;
          P2(5,i,:)= (AVGCOSTS3(2,1)-AVGCOSTS3(3,1)*Sub4)*0.66*Z;
          P3(5,i,:)= AVGCOSTS3(3,1)*0.66*(1-Sub4)*Z;          
      end
      if Tier2(i,1) == 6
          P1(5,i,:)= (AVGCOSTC3(1,1)-AVGCOSTC3(3,1)*Sub4)*0.66*Z1;
          P2(5,i,:)= (AVGCOSTC3(2,1)-AVGCOSTC3(3,1)*Sub4)*0.66*Z1;
          P3(5,i,:)= AVGCOSTC3(3,1)*0.66*(1-Sub4)*Z1;          
      end      
      if Tier2(i,1) == 8
          P1(5,i,:)= (AVGCOSTF3(1,1)-AVGCOSTF3(3,1)*Sub4)*0.66*Z;
          P2(5,i,:)= (AVGCOSTF3(2,1)-AVGCOSTF3(3,1)*Sub4)*0.66*Z;
          P3(5,i,:)= AVGCOSTF3(3,1)*0.66*(1-Sub4)*Z;          
      end
   end   
   
   if Inc2(i,1) == 5
      if Tier2(i,1) == 1
          P1(5,i,:)= (AVGCOSTE3(1,1)-AVGCOSTE3(3,1)*Sub5)*0.61*Z;
          P2(5,i,:)= (AVGCOSTE3(2,1)-AVGCOSTE3(3,1)*Sub5)*0.61*Z;
          P3(5,i,:)= AVGCOSTE3(3,1)*0.61*(1-Sub5)*Z;          
      end
      if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)
          P1(5,i,:)= (AVGCOSTS3(1,1)-AVGCOSTS3(3,1)*Sub5)*0.61*Z;
          P2(5,i,:)= (AVGCOSTS3(2,1)-AVGCOSTS3(3,1)*Sub5)*0.61*Z;
          P3(5,i,:)= AVGCOSTS3(3,1)*0.61*(1-Sub5)*Z;          
      end
      if Tier2(i,1) == 6
          P1(5,i,:)= (AVGCOSTC3(1,1)-AVGCOSTC3(3,1)*Sub5)*0.61*Z1;
          P2(5,i,:)= (AVGCOSTC3(2,1)-AVGCOSTC3(3,1)*Sub5)*0.61*Z1;
          P3(5,i,:)= AVGCOSTC3(3,1)*0.61*(1-Sub5)*Z1;          
      end      
      if Tier2(i,1) == 8
          P1(5,i,:)= (AVGCOSTF3(1,1)-AVGCOSTF3(3,1)*Sub5)*0.61*Z;
          P2(5,i,:)= (AVGCOSTF3(2,1)-AVGCOSTF3(3,1)*Sub5)*0.61*Z;
          P3(5,i,:)= AVGCOSTF3(3,1)*0.61*(1-Sub5)*Z;          
      end
   end    
 
end


P1(5,:,:) = max(reshape(P1(5,:,:),nIs,Sim),150);
P2(5,:,:) = max(reshape(P2(5,:,:),nIs,Sim),150);
P3(5,:,:) = max(reshape(P3(5,:,:),nIs,Sim),150);

    P3E4 = zeros(3,5);
    P3S4 = zeros(3,5);
    P3C4 = zeros(3,5);
    P3F4 = zeros(3,5);
    
    P3E4(1,1) = (AVGCOSTE3(1,1)-AVGCOSTE3(3,1)*Sub1)*Z;
    P3E4(2,1) = (AVGCOSTE3(2,1)-AVGCOSTE3(3,1)*Sub1)*Z;
    P3E4(3,1) =  AVGCOSTE3(3,1)*(1-Sub1)*Z;          
      
    P3S4(1,1) = (AVGCOSTS3(1,1)-AVGCOSTS3(3,1)*Sub1)*Z;
    P3S4(2,1) = (AVGCOSTS3(2,1)-AVGCOSTS3(3,1)*Sub1)*Z;
    P3S4(3,1) =  AVGCOSTS3(3,1)*(1-Sub1)*Z;          
   
    P3C4(1,1) = (AVGCOSTC3(1,1)-AVGCOSTC3(3,1)*Sub1)*Z1;
    P3C4(2,1) = (AVGCOSTC3(2,1)-AVGCOSTC3(3,1)*Sub1)*Z1;
    P3C4(3,1) =  AVGCOSTC3(3,1)*(1-Sub1)*Z1;          
    
    P3F4(1,1) = (AVGCOSTF3(1,1)-AVGCOSTF3(3,1)*Sub1)*Z;
    P3F4(2,1) = (AVGCOSTF3(2,1)-AVGCOSTF3(3,1)*Sub1)*Z;
    P3F4(3,1) =  AVGCOSTF3(3,1)*(1-Sub1)*Z;          
  
    P3E4(1,2) = (AVGCOSTE3(1,1)-AVGCOSTE3(3,1)*Sub2)*Z;
    P3E4(2,2) = (AVGCOSTE3(2,1)-AVGCOSTE3(3,1)*Sub2)*Z;
    P3E4(3,2) =  AVGCOSTE3(3,1)*(1-Sub2)*Z;          
      
    P3S4(1,2) = (AVGCOSTS3(1,1)-AVGCOSTS3(3,1)*Sub2)*Z;
    P3S4(2,2) = (AVGCOSTS3(2,1)-AVGCOSTS3(3,1)*Sub2)*Z;
    P3S4(3,2) = AVGCOSTS3(3,1)*(1-Sub2)*Z;           
   
    P3C4(1,2) = (AVGCOSTC3(1,1)-AVGCOSTC3(3,1)*Sub2)*Z1;
    P3C4(2,2) = (AVGCOSTC3(2,1)-AVGCOSTC3(3,1)*Sub2)*Z1;
    P3C4(3,2) =  AVGCOSTC3(3,1)*(1-Sub2)*Z1;         
    
    P3F4(1,2) = (AVGCOSTF3(1,1)-AVGCOSTF3(3,1)*Sub2)*Z;
    P3F4(2,2) = (AVGCOSTF3(2,1)-AVGCOSTF3(3,1)*Sub2)*Z;
    P3F4(3,2) =  AVGCOSTF3(3,1)*(1-Sub2)*Z;      

    P3E4(1,3) = (AVGCOSTE3(1,1)-AVGCOSTE3(3,1)*Sub3)*Z;
    P3E4(2,3) = (AVGCOSTE3(2,1)-AVGCOSTE3(3,1)*Sub3)*Z;
    P3E4(3,3) =  AVGCOSTE3(3,1)*(1-Sub3)*Z;         
      
    P3S4(1,3) = (AVGCOSTS3(1,1)-AVGCOSTS3(3,1)*Sub3)*Z;
    P3S4(2,3) = (AVGCOSTS3(2,1)-AVGCOSTS3(3,1)*Sub3)*Z;
    P3S4(3,3) =  AVGCOSTS3(3,1)*(1-Sub3)*Z;          
   
    P3C4(1,3) = (AVGCOSTC3(1,1)-AVGCOSTC3(3,1)*Sub3)*Z1;
    P3C4(2,3) = (AVGCOSTC3(2,1)-AVGCOSTC3(3,1)*Sub3)*Z1;
    P3C4(3,3) = AVGCOSTC3(3,1)*(1-Sub3)*Z1;        
    
    P3F4(1,3) = (AVGCOSTF3(1,1)-AVGCOSTF3(3,1)*Sub3)*Z;
    P3F4(2,3) = (AVGCOSTF3(2,1)-AVGCOSTF3(3,1)*Sub3)*Z;
    P3F4(3,3) = AVGCOSTF3(3,1)*(1-Sub3)*Z; 
            
    P3E4(1,4) = (AVGCOSTE3(1,1)-AVGCOSTE3(3,1)*Sub4)*Z;
    P3E4(2,4) = (AVGCOSTE3(2,1)-AVGCOSTE3(3,1)*Sub4)*Z;
    P3E4(3,4) = AVGCOSTE3(3,1)*(1-Sub4)*Z;         
      
    P3S4(1,4) = (AVGCOSTS3(1,1)-AVGCOSTS3(3,1)*Sub4)*Z;
    P3S4(2,4) = (AVGCOSTS3(2,1)-AVGCOSTS3(3,1)*Sub4)*Z;
    P3S4(3,4) =  AVGCOSTS3(3,1)*(1-Sub4)*Z;           
   
    P3C4(1,4) = (AVGCOSTC3(1,1)-AVGCOSTC3(3,1)*Sub4)*Z1;
    P3C4(2,4) = (AVGCOSTC3(2,1)-AVGCOSTC3(3,1)*Sub4)*Z1;
    P3C4(3,4) =  AVGCOSTC3(3,1)*(1-Sub4)*Z1;        
    
    P3F4(1,4) = (AVGCOSTF3(1,1)-AVGCOSTF3(3,1)*Sub4)*Z;
    P3F4(2,4) = (AVGCOSTF3(2,1)-AVGCOSTF3(3,1)*Sub4)*Z;
    P3F4(3,4) = AVGCOSTF3(3,1)*(1-Sub4)*Z;

    P3E4(1,5) = (AVGCOSTE3(1,1)-AVGCOSTE3(3,1)*Sub5)*Z;
    P3E4(2,5) = (AVGCOSTE3(2,1)-AVGCOSTE3(3,1)*Sub5)*Z;
    P3E4(3,5) =  AVGCOSTE3(3,1)*(1-Sub5)*Z;        
      
    P3S4(1,5) = (AVGCOSTS3(1,1)-AVGCOSTS3(3,1)*Sub5)*Z;
    P3S4(2,5) = (AVGCOSTS3(2,1)-AVGCOSTS3(3,1)*Sub5)*Z;
    P3S4(3,5) = AVGCOSTS3(3,1)*(1-Sub5)*Z;             
   
    P3C4(1,5) = (AVGCOSTC3(1,1)-AVGCOSTC3(3,1)*Sub5)*Z1;
    P3C4(2,5) = (AVGCOSTC3(2,1)-AVGCOSTC3(3,1)*Sub5)*Z1;
    P3C4(3,5) =  AVGCOSTC3(3,1)*(1-Sub5)*Z1;       
    
    P3F4(1,5) = (AVGCOSTF3(1,1)-AVGCOSTF3(3,1)*Sub5)*Z;
    P3F4(2,5) = (AVGCOSTF3(2,1)-AVGCOSTF3(3,1)*Sub5)*Z;
    P3F4(3,5) = AVGCOSTF3(3,1)*(1-Sub5)*Z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P3F4 = max(P3F4,150);
P3S4 = max(P3S4,150);
P3C4 = max(P3C4,150);
P3E4 = max(P3E4,150);

%%%%%%%%% Set Incumbent Plans for Next Year

for i = 1:nIs
    for s = 1:Sim
       
       if CC3(i,s)==1
            choice15(i,s) = 1;
       end
       if CC3(i,s)==2
            choice25(i,s) = 1;
       end
       if CC3(i,s)==3
            choice35(i,s) = 1;
       end
    
    end
end

for k = 1:K
    SwC15(:,k,:) = SCC3.*choice15;
    SwC25(:,k,:) = SCC3.*choice25;
    SwC35(:,k,:) = SCC3.*choice35;
end


%%%% New expected utility baselines calculated because of new premiums and inertia 

for i = 1:nIs
    for k = 1:K
        EUvector15(i,k,:) = max(0.1,(W - reshape(PPO250OP3(i,k,:),1,Sim) - reshape(P1(5,i,:),1,Sim)));
        EUvector25(i,k,:) = max(0.1,(W - reshape(PPO500OP3(i,k,:),1,Sim) - reshape(P2(5,i,:),1,Sim)));
        EUvector35(i,k,:) = max(0.1,(W - reshape(PPO1200OP3(i,k,:),1,Sim) - reshape(P3(5,i,:),1,Sim)));
    end
end

EUvector15 = exp(-((EUvector15 +SwC15).*max(RA,0.0000001)));
EUvector25 = exp(-((EUvector25 + SwC25 + CHTC + Epsilon5005).*max(RA,0.0000001)));
EUvector35 = exp(-((EUvector35 + CDHP + CHTC + SwC35 + Epsilon12005).*max(RA, 0.0000001)));

for n = 1:Sim
    EUvectorD(:,1,n) = (1/K)*sum(reshape(EUvector15(:,:,n),nIs,K),2);
    EUvectorD(:,2,n) = (1/K)*sum(reshape(EUvector25(:,:,n),nIs,K),2);
    EUvectorD(:,3,n) = (1/K)*sum(reshape(EUvector35(:,:,n),nIs,K),2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Year 5 Chocies %%%%%%%%%%%%%%%%%%%%%%%%%%%

CC4 = zeros(nIs,Sim);
MktShare4 = zeros(nPlans,1);

TOTALCOSTE4 = zeros(3,1);
TOTALNUME4 = zeros(3,1);
TOTALNUME4i = zeros(3,nIs);

TOTALCOSTS4 = zeros(3,1);
TOTALNUMS4 = zeros(3,1);
TOTALNUMS4i = zeros(3,nIs);

TOTALCOSTC4 = zeros(3,1);
TOTALNUMC4 = zeros(3,1);
TOTALNUMC4i = zeros(3,nIs);

TOTALCOSTF4 = zeros(3,1);
TOTALNUMF4 = zeros(3,1);
TOTALNUMF4i = zeros(3,nIs);

for i = 1:nIs 
    
    PT5 = ones(nPlans,Sim);
    PT5(1:nPlans-1,:) =  (reshape(EUvectorD(i,3,:),1,Sim)'*ones(1,nPlans-1))'./reshape(EUvectorD(i,[1:nPlans-1],:),nPlans-1,Sim);
    PTAgg5 = zeros(1,Sim);
    PTAgg5 = sum(PT5,1);
    
    PT5(1,:) = PT5(1,:)./PTAgg5;
    PT5(2,:) = PT5(2,:)./PTAgg5;
    PT5(3,:) = PT5(3,:)./PTAgg5;
    
    PT5 = ((PT5'.^3)./(sum((PT5'.^3),2)*ones(1,nPlans)))';
    
    UNI = unifrnd(0,1,1,Sim);
    
    for j = 1:Sim
        if PT5(1,j) >= UNI(1,j)
           CC4(i,j)=1;
        elseif PT5(1,j)+PT5(2,j) >= UNI(1,j)
           CC4(i,j)=2;
        elseif (1-PT5(3,j)) < UNI(1,j)
           CC4(i,j)=3;
        end
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nIs
    for j = 1:Sim
      if Tier2(i,1) == 1  
        if CC4(i,j) == 1 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUME4i(1,i) = TOTALNUME4i(1,i)+1; 
        end
        if CC4(i,j) == 2 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUME4i(2,i) = TOTALNUME4i(2,i)+1; 
        end
        if CC4(i,j) == 3 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUME4i(3,i) = TOTALNUME4i(3,i)+1; 
        end
      end
      
    if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)  
        if CC4(i,j) == 1 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMS4i(1,i) = TOTALNUMS4i(1,i)+1; 
        end
        if CC4(i,j) == 2 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMS4i(2,i) = TOTALNUMS4i(2,i)+1; 
        end
        if CC4(i,j) == 3 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMS4i(3,i) = TOTALNUMS4i(3,i)+1; 
        end
    end    

    if Tier2(i,1) == 6  
        if CC4(i,j) == 1 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMC4i(1,i) = TOTALNUMC4i(1,i)+1; 
        end
        if CC4(i,j) == 2 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMC4i(2,i) = TOTALNUMC4i(2,i)+1; 
        end
        if CC4(i,j) == 3 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMC4i(3,i) = TOTALNUMC4i(3,i)+1; 
        end
    end   
    
   if Tier2(i,1) == 8  
        if CC4(i,j) == 1 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMF4i(1,i) = TOTALNUMF4i(1,i)+1; 
        end
        if CC4(i,j) == 2 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMF4i(2,i) = TOTALNUMF4i(2,i)+1; 
        end
        if CC4(i,j) == 3 & (AcceptRejectF(i,j)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,j)==1))
           TOTALNUMF4i(3,i) = TOTALNUMF4i(3,i)+1; 
        end
    end  
    end
end

for i = 1:nIs
    if Tier2(i,1)==1
    TOTALNUME4 = TOTALNUME4 + TOTALNUME4i(:,i)*(Sim/sum(TOTALNUME4i(:,i)));    
    end
    if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)  
    TOTALNUMS4 = TOTALNUMS4 + TOTALNUMS4i(:,i)*(Sim/sum(TOTALNUMS4i(:,i)));
    end
    if Tier2(i,1)==6
    TOTALNUMC4 = TOTALNUMC4 + TOTALNUMC4i(:,i)*(Sim/sum(TOTALNUMC4i(:,i)));
    end
    if Tier2(i,1)==8
    TOTALNUMF4 = TOTALNUMF4 + TOTALNUMF4i(:,i)*(Sim/sum(TOTALNUMF4i(:,i)));
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nIs
    if Tier2(i,1) == 1
    TOTALCOSTE4(1,1) = TOTALCOSTE4(1,1) + (1/K)*sum(reshape(PlanPaid3(i,:,1),1,K),2)*TOTALNUME4i(1,i)*(Sim/sum(TOTALNUME4i(:,i)));
    TOTALCOSTE4(2,1) = TOTALCOSTE4(2,1) + (1/K)*sum(reshape(PlanPaid3(i,:,2),1,K),2)*TOTALNUME4i(2,i)*(Sim/sum(TOTALNUME4i(:,i)));
    TOTALCOSTE4(3,1) = TOTALCOSTE4(3,1) + (1/K)*sum(reshape(PlanPaid3(i,:,3),1,K),2)*TOTALNUME4i(3,i)*(Sim/sum(TOTALNUME4i(:,i)));
    end
   
    if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)  
    TOTALCOSTS4(1,1) = TOTALCOSTS4(1,1) + (1/K)*sum(reshape(PlanPaid3(i,:,1),1,K),2)*TOTALNUMS4i(1,i)*(Sim/sum(TOTALNUMS4i(:,i)));
    TOTALCOSTS4(2,1) = TOTALCOSTS4(2,1) + (1/K)*sum(reshape(PlanPaid3(i,:,2),1,K),2)*TOTALNUMS4i(2,i)*(Sim/sum(TOTALNUMS4i(:,i)));
    TOTALCOSTS4(3,1) = TOTALCOSTS4(3,1) + (1/K)*sum(reshape(PlanPaid3(i,:,3),1,K),2)*TOTALNUMS4i(3,i)*(Sim/sum(TOTALNUMS4i(:,i)));
    end
    
    if Tier2(i,1)==6
    TOTALCOSTC4(1,1) = TOTALCOSTC4(1,1) + (1/K)*sum(reshape(PlanPaid3(i,:,1),1,K),2)*TOTALNUMC4i(1,i)*(Sim/sum(TOTALNUMC4i(:,i)));
    TOTALCOSTC4(2,1) = TOTALCOSTC4(2,1) + (1/K)*sum(reshape(PlanPaid3(i,:,2),1,K),2)*TOTALNUMC4i(2,i)*(Sim/sum(TOTALNUMC4i(:,i)));
    TOTALCOSTC4(3,1) = TOTALCOSTC4(3,1) + (1/K)*sum(reshape(PlanPaid3(i,:,3),1,K),2)*TOTALNUMC4i(3,i)*(Sim/sum(TOTALNUMC4i(:,i)));
    end    
       
   if Tier2(i,1)==8
   TOTALCOSTF4(1,1) = TOTALCOSTF4(1,1) + (1/K)*sum(reshape(PlanPaid3(i,:,1),1,K),2)*TOTALNUMF4i(1,i)*(Sim/sum(TOTALNUMF4i(:,i)));
   TOTALCOSTF4(2,1) = TOTALCOSTF4(2,1) + (1/K)*sum(reshape(PlanPaid3(i,:,2),1,K),2)*TOTALNUMF4i(2,i)*(Sim/sum(TOTALNUMF4i(:,i)));
   TOTALCOSTF4(3,1) = TOTALCOSTF4(3,1) + (1/K)*sum(reshape(PlanPaid3(i,:,3),1,K),2)*TOTALNUMF4i(3,i)*(Sim/sum(TOTALNUMF4i(:,i)));
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AVGCOSTE4 = zeros(3,1);
AVGCOSTS4 = zeros(3,1);
AVGCOSTC4 = zeros(3,1);
AVGCOSTF4 = zeros(3,1);

%%%%%%%% Incorporate Adminsitrative Loadinig Factors and Risk Buffer Funds

AVGCOSTE4(1,1) = TOTALCOSTE4(1,1)/TOTALNUME4(1,1);
AVGCOSTE4(2,1) = TOTALCOSTE4(2,1)/TOTALNUME4(2,1);
AVGCOSTE4(3,1) = TOTALCOSTE4(3,1)/TOTALNUME4(3,1);

AVGCOSTS4(1,1) = TOTALCOSTS4(1,1)/TOTALNUMS4(1,1);
AVGCOSTS4(2,1) = TOTALCOSTS4(2,1)/TOTALNUMS4(2,1);
AVGCOSTS4(3,1) = TOTALCOSTS4(3,1)/TOTALNUMS4(3,1);

AVGCOSTC4(1,1) = TOTALCOSTC4(1,1)/TOTALNUMC4(1,1);
AVGCOSTC4(2,1) = TOTALCOSTC4(2,1)/TOTALNUMC4(2,1);
AVGCOSTC4(3,1) = TOTALCOSTC4(3,1)/TOTALNUMC4(3,1);

AVGCOSTF4(1,1) = TOTALCOSTF4(1,1)/TOTALNUMF4(1,1);
AVGCOSTF4(2,1) = TOTALCOSTF4(2,1)/TOTALNUMF4(2,1);
AVGCOSTF4(3,1) = TOTALCOSTF4(3,1)/TOTALNUMF4(3,1);

MktShare4 = zeros(3,1);

MktShare4(1,1) = (TOTALNUMF4(1,1) + TOTALNUMC4(1,1) + TOTALNUMS4(1,1) + TOTALNUME4(1,1))/Sim;
MktShare4(2,1) = (TOTALNUMF4(2,1) + TOTALNUMC4(2,1) + TOTALNUMS4(2,1) + TOTALNUME4(2,1))/Sim;
MktShare4(3,1) = (TOTALNUMF4(3,1) + TOTALNUMC4(3,1) + TOTALNUMS4(3,1) + TOTALNUME4(3,1))/Sim;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nIs

   if Inc2(i,1) == 1
      if Tier2(i,1) == 1
          P1(6,i,:)= (AVGCOSTE4(1,1)-AVGCOSTE4(3,1)*Sub1)*0.82*Z;
          P2(6,i,:)= (AVGCOSTE4(2,1)-AVGCOSTE4(3,1)*Sub1)*0.82*Z;
          P3(6,i,:)= AVGCOSTE4(3,1)*0.82*(1-Sub1)*Z;          
      end
      if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)
          P1(6,i,:)= (AVGCOSTS4(1,1)-AVGCOSTS4(3,1)*Sub1)*0.82*Z;
          P2(6,i,:)= (AVGCOSTS4(2,1)-AVGCOSTS4(3,1)*Sub1)*0.82*Z;
          P3(6,i,:)= AVGCOSTS4(3,1)*0.82*(1-Sub1)*Z;          
      end
      if Tier2(i,1) == 6
          P1(6,i,:)= (AVGCOSTC4(1,1)-AVGCOSTC4(3,1)*Sub1)*0.82*Z1;
          P2(6,i,:)= (AVGCOSTC4(2,1)-AVGCOSTC4(3,1)*Sub1)*0.82*Z1;
          P3(6,i,:)= AVGCOSTC4(3,1)*0.82*(1-Sub1)*Z1;          
      end      
      if Tier2(i,1) == 8
          P1(6,i,:)= (AVGCOSTF4(1,1)-AVGCOSTF4(3,1)*Sub1)*0.82*Z;
          P2(6,i,:)= (AVGCOSTF4(2,1)-AVGCOSTF4(3,1)*Sub1)*0.82*Z;
          P3(6,i,:)= AVGCOSTF4(3,1)*0.82*(1-Sub1)*Z;          
      end
   end
   
   if Inc2(i,1) == 2
      if Tier2(i,1) == 1
          P1(6,i,:)= (AVGCOSTE4(1,1)-AVGCOSTE4(3,1)*Sub2)*0.73*Z;
          P2(6,i,:)= (AVGCOSTE4(2,1)-AVGCOSTE4(3,1)*Sub2)*0.73*Z;
          P3(6,i,:)= AVGCOSTE4(3,1)*0.73*(1-Sub2)*Z;          
      end
      if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)
          P1(6,i,:)= (AVGCOSTS4(1,1)-AVGCOSTS4(3,1)*Sub2)*0.80*Z;
          P2(6,i,:)= (AVGCOSTS4(2,1)-AVGCOSTS4(3,1)*Sub2)*0.80*Z;
          P3(6,i,:)= AVGCOSTS4(3,1)*0.80*(1-Sub2)*Z;          
      end
      if Tier2(i,1) == 6
          P1(6,i,:)= (AVGCOSTC4(1,1)-AVGCOSTC4(3,1)*Sub2)*0.80*Z1;
          P2(6,i,:)= (AVGCOSTC4(2,1)-AVGCOSTC4(3,1)*Sub2)*0.80*Z1;
          P3(6,i,:)= AVGCOSTC4(3,1)*0.80*(1-Sub2)*Z1;          
      end      
      if Tier2(i,1) == 8
          P1(6,i,:)= (AVGCOSTF4(1,1)-AVGCOSTF4(3,1)*Sub2)*0.80*Z;
          P2(6,i,:)= (AVGCOSTF4(2,1)-AVGCOSTF4(3,1)*Sub2)*0.80*Z;
          P3(6,i,:)= AVGCOSTF4(3,1)*0.80*(1-Sub2)*Z;          
      end
   end
   if Inc2(i,1) == 3
      if Tier2(i,1) == 1
          P1(6,i,:)= (AVGCOSTE4(1,1)-AVGCOSTE4(3,1)*Sub3)*0.69*Z;
          P2(6,i,:)= (AVGCOSTE4(2,1)-AVGCOSTE4(3,1)*Sub3)*0.69*Z;
          P3(6,i,:)= AVGCOSTE4(3,1)*0.69*(1-Sub3)*Z;          
      end
      if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)
          P1(6,i,:)= (AVGCOSTS4(1,1)-AVGCOSTS4(3,1)*Sub3)*0.69*Z;
          P2(6,i,:)= (AVGCOSTS4(2,1)-AVGCOSTS4(3,1)*Sub3)*0.69*Z;
          P3(6,i,:)= AVGCOSTS4(3,1)*0.69*(1-Sub3)*Z;          
      end
      if Tier2(i,1) == 6
          P1(6,i,:)= (AVGCOSTC4(1,1)-AVGCOSTC4(3,1)*Sub3)*0.69*Z1;
          P2(6,i,:)= (AVGCOSTC4(2,1)-AVGCOSTC4(3,1)*Sub3)*0.69*Z1;
          P3(6,i,:)= AVGCOSTC4(3,1)*0.69*(1-Sub3)*Z1;          
      end      
      if Tier2(i,1) == 8
          P1(6,i,:)= (AVGCOSTF4(1,1)-AVGCOSTF4(3,1)*Sub3)*0.69*Z;
          P2(6,i,:)= (AVGCOSTF4(2,1)-AVGCOSTF4(3,1)*Sub3)*0.69*Z;
          P3(6,i,:)= AVGCOSTF4(3,1)*0.69*(1-Sub3)*Z;          
      end
   end
   
   if Inc2(i,1) == 4
      if Tier2(i,1) == 1
          P1(6,i,:)= (AVGCOSTE4(1,1)-AVGCOSTE4(3,1)*Sub4)*0.66*Z;
          P2(6,i,:)= (AVGCOSTE4(2,1)-AVGCOSTE4(3,1)*Sub4)*0.66*Z;
          P3(6,i,:)= AVGCOSTE4(3,1)*0.66*(1-Sub4)*Z;          
      end
      if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)
          P1(6,i,:)= (AVGCOSTS4(1,1)-AVGCOSTS4(3,1)*Sub4)*0.66*Z;
          P2(6,i,:)= (AVGCOSTS4(2,1)-AVGCOSTS4(3,1)*Sub4)*0.66*Z;
          P3(6,i,:)= AVGCOSTS4(3,1)*0.66*(1-Sub4)*Z;          
      end
      if Tier2(i,1) == 6
          P1(6,i,:)= (AVGCOSTC4(1,1)-AVGCOSTC4(3,1)*Sub4)*0.66*Z1;
          P2(6,i,:)= (AVGCOSTC4(2,1)-AVGCOSTC4(3,1)*Sub4)*0.66*Z1;
          P3(6,i,:)= AVGCOSTC4(3,1)*0.66*(1-Sub4)*Z1;          
      end      
      if Tier2(i,1) == 8
          P1(6,i,:)= (AVGCOSTF4(1,1)-AVGCOSTF4(3,1)*Sub4)*0.66*Z;
          P2(6,i,:)= (AVGCOSTF4(2,1)-AVGCOSTF4(3,1)*Sub4)*0.66*Z;
          P3(6,i,:)= AVGCOSTF4(3,1)*0.66*(1-Sub4)*Z;          
      end
   end   
   
   if Inc2(i,1) == 5
      if Tier2(i,1) == 1
          P1(6,i,:)= (AVGCOSTE4(1,1)-AVGCOSTE4(3,1)*Sub5)*0.61*Z;
          P2(6,i,:)= (AVGCOSTE4(2,1)-AVGCOSTE4(3,1)*Sub5)*0.61*Z;
          P3(6,i,:)= AVGCOSTE4(3,1)*0.61*(1-Sub5)*Z;          
      end
      if (Tier2(i,1) == 2 | Tier2(i,1) ==11 | Tier2(i,1)==12)
          P1(6,i,:)= (AVGCOSTS4(1,1)-AVGCOSTS4(3,1)*Sub5)*0.61*Z;
          P2(6,i,:)= (AVGCOSTS4(2,1)-AVGCOSTS4(3,1)*Sub5)*0.61*Z;
          P3(6,i,:)= AVGCOSTS4(3,1)*0.61*(1-Sub5)*Z;          
      end
      if Tier2(i,1) == 6
          P1(6,i,:)= (AVGCOSTC4(1,1)-AVGCOSTC4(3,1)*Sub5)*0.61*Z1;
          P2(6,i,:)= (AVGCOSTC4(2,1)-AVGCOSTC4(3,1)*Sub5)*0.61*Z1;
          P3(6,i,:)= AVGCOSTC4(3,1)*0.61*(1-Sub5)*Z1;          
      end      
      if Tier2(i,1) == 8
          P1(6,i,:)= (AVGCOSTF4(1,1)-AVGCOSTF4(3,1)*Sub5)*0.61*Z;
          P2(6,i,:)= (AVGCOSTF4(2,1)-AVGCOSTF4(3,1)*Sub5)*0.61*Z;
          P3(6,i,:)= AVGCOSTF4(3,1)*0.61*(1-Sub5)*Z;          
      end
   end    
 
end


P1(6,:,:) = max(reshape(P1(6,:,:),nIs,Sim),150);
P2(6,:,:) = max(reshape(P2(6,:,:),nIs,Sim),150);
P3(6,:,:) = max(reshape(P3(6,:,:),nIs,Sim),150);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 P3E5 = zeros(3,5);
    P3S5 = zeros(3,5);
    P3C5 = zeros(3,5);
    P3F5 = zeros(3,5);
    
    P3E5(1,1) = (AVGCOSTE4(1,1)-AVGCOSTE4(3,1)*Sub1)*Z;
    P3E5(2,1) = (AVGCOSTE4(2,1)-AVGCOSTE4(3,1)*Sub1)*Z;
    P3E5(3,1) =  AVGCOSTE4(3,1)*(1-Sub1)*Z;          
      
    P3S5(1,1) = (AVGCOSTS4(1,1)-AVGCOSTS4(3,1)*Sub1)*Z;
    P3S5(2,1) = (AVGCOSTS4(2,1)-AVGCOSTS4(3,1)*Sub1)*Z;
    P3S5(3,1) =  AVGCOSTS4(3,1)*(1-Sub1)*Z;          
   
    P3C5(1,1) = (AVGCOSTC4(1,1)-AVGCOSTC4(3,1)*Sub1)*Z1;
    P3C5(2,1) = (AVGCOSTC4(2,1)-AVGCOSTC4(3,1)*Sub1)*Z1;
    P3C5(3,1) =  AVGCOSTC4(3,1)*(1-Sub1)*Z1;          
    
    P3F5(1,1) = (AVGCOSTF4(1,1)-AVGCOSTF4(3,1)*Sub1)*Z;
    P3F5(2,1) = (AVGCOSTF4(2,1)-AVGCOSTF4(3,1)*Sub1)*Z;
    P3F5(3,1) =  AVGCOSTF4(3,1)*(1-Sub1)*Z;          
  
    P3E5(1,2) = (AVGCOSTE4(1,1)-AVGCOSTE4(3,1)*Sub2)*Z;
    P3E5(2,2) = (AVGCOSTE4(2,1)-AVGCOSTE4(3,1)*Sub2)*Z;
    P3E5(3,2) =  AVGCOSTE4(3,1)*(1-Sub2)*Z;          
      
    P3S5(1,2) = (AVGCOSTS4(1,1)-AVGCOSTS4(3,1)*Sub2)*Z;
    P3S5(2,2) = (AVGCOSTS4(2,1)-AVGCOSTS4(3,1)*Sub2)*Z;
    P3S5(3,2) = AVGCOSTS4(3,1)*(1-Sub2)*Z;           
   
    P3C5(1,2) = (AVGCOSTC4(1,1)-AVGCOSTC4(3,1)*Sub2)*Z1;
    P3C5(2,2) = (AVGCOSTC4(2,1)-AVGCOSTC4(3,1)*Sub2)*Z1;
    P3C5(3,2) =  AVGCOSTC4(3,1)*(1-Sub2)*Z1;         
    
    P3F5(1,2) = (AVGCOSTF4(1,1)-AVGCOSTF4(3,1)*Sub2)*Z;
    P3F5(2,2) = (AVGCOSTF4(2,1)-AVGCOSTF4(3,1)*Sub2)*Z;
    P3F5(3,2) =  AVGCOSTF4(3,1)*(1-Sub2)*Z;      

    P3E5(1,3) = (AVGCOSTE4(1,1)-AVGCOSTE4(3,1)*Sub3)*Z;
    P3E5(2,3) = (AVGCOSTE4(2,1)-AVGCOSTE4(3,1)*Sub3)*Z;
    P3E5(3,3) =  AVGCOSTE4(3,1)*(1-Sub3)*Z;         
      
    P3S5(1,3) = (AVGCOSTS4(1,1)-AVGCOSTS4(3,1)*Sub3)*Z;
    P3S5(2,3) = (AVGCOSTS4(2,1)-AVGCOSTS4(3,1)*Sub3)*Z;
    P3S5(3,3) =  AVGCOSTS4(3,1)*(1-Sub3)*Z;          
   
    P3C5(1,3) = (AVGCOSTC4(1,1)-AVGCOSTC4(3,1)*Sub3)*Z1;
    P3C5(2,3) = (AVGCOSTC4(2,1)-AVGCOSTC4(3,1)*Sub3)*Z1;
    P3C5(3,3) = AVGCOSTC4(3,1)*(1-Sub3)*Z1;        
    
    P3F5(1,3) = (AVGCOSTF4(1,1)-AVGCOSTF4(3,1)*Sub3)*Z;
    P3F5(2,3) = (AVGCOSTF4(2,1)-AVGCOSTF4(3,1)*Sub3)*Z;
    P3F5(3,3) = AVGCOSTF4(3,1)*(1-Sub3)*Z; 
            
    P3E5(1,4) = (AVGCOSTE4(1,1)-AVGCOSTE4(3,1)*Sub4)*Z;
    P3E5(2,4) = (AVGCOSTE4(2,1)-AVGCOSTE4(3,1)*Sub4)*Z;
    P3E5(3,4) = AVGCOSTE4(3,1)*(1-Sub4)*Z;         
      
    P3S5(1,4) = (AVGCOSTS4(1,1)-AVGCOSTS4(3,1)*Sub4)*Z;
    P3S5(2,4) = (AVGCOSTS4(2,1)-AVGCOSTS4(3,1)*Sub4)*Z;
    P3S5(3,4) =  AVGCOSTS4(3,1)*(1-Sub4)*Z;           
   
    P3C5(1,4) = (AVGCOSTC4(1,1)-AVGCOSTC4(3,1)*Sub4)*Z1;
    P3C5(2,4) = (AVGCOSTC4(2,1)-AVGCOSTC4(3,1)*Sub4)*Z1;
    P3C5(3,4) =  AVGCOSTC4(3,1)*(1-Sub4)*Z1;        
    
    P3F5(1,4) = (AVGCOSTF4(1,1)-AVGCOSTF4(3,1)*Sub4)*Z;
    P3F5(2,4) = (AVGCOSTF4(2,1)-AVGCOSTF4(3,1)*Sub4)*Z;
    P3F5(3,4) = AVGCOSTF4(3,1)*(1-Sub4)*Z;

    P3E5(1,5) = (AVGCOSTE4(1,1)-AVGCOSTE4(3,1)*Sub5)*Z;
    P3E5(2,5) = (AVGCOSTE4(2,1)-AVGCOSTE4(3,1)*Sub5)*Z;
    P3E5(3,5) =  AVGCOSTE4(3,1)*(1-Sub5)*Z;        
      
    P3S5(1,5) = (AVGCOSTS4(1,1)-AVGCOSTS4(3,1)*Sub5)*Z;
    P3S5(2,5) = (AVGCOSTS4(2,1)-AVGCOSTS4(3,1)*Sub5)*Z;
    P3S5(3,5) = AVGCOSTS4(3,1)*(1-Sub5)*Z;             
   
    P3C5(1,5) = (AVGCOSTC4(1,1)-AVGCOSTC4(3,1)*Sub5)*Z1;
    P3C5(2,5) = (AVGCOSTC4(2,1)-AVGCOSTC4(3,1)*Sub5)*Z1;
    P3C5(3,5) =  AVGCOSTC4(3,1)*(1-Sub5)*Z1;       
    
    P3F5(1,5) = (AVGCOSTF4(1,1)-AVGCOSTF4(3,1)*Sub5)*Z;
    P3F5(2,5) = (AVGCOSTF4(2,1)-AVGCOSTF4(3,1)*Sub5)*Z;
    P3F5(3,5) = AVGCOSTF4(3,1)*(1-Sub5)*Z;

P3F5 = max(P3F5,150);
P3S5 = max(P3S5,150);
P3C5 = max(P3C5,150);
P3E5 = max(P3E5,150);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%  FIVE YEARS OF CHOICES DONE FOR COUNTERFACTUAL SIMULATIONS. REPORT CHOICE RESULTS NOW
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Since results can change a little with simulation draws, here we just 
  % set up reporting of results for this specific counteractual with simulated data
  % and user can report results themselves given their draws of simulated matrices / random coefficients / epsilons.
  % In actual analysis simulate very large simulation matrix and store for comparitive statics / primary analysis 
  % show in section 6 of main paper. 	

  % NOTE: Results different from actual results, reflecting differences between simulated and actual data. Simulated data 
  % have lower inertia, higher risk preferences, and no demographics correlation structure which makes them quite different 
  % for purposes of this CF. Code / results here to illustrate use of code so it can be used in other settings. 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MktShare Overall for Plans %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MktShare0
MktShare1
MktShare2
MktShare3
MktShare4 
 
%%%%%%%%%%%%% Actual premiums paid each year 
 
%%% Family prices after first year

P3FZ
P3F
P3F3
P3F4
P3F5

%%% With child(ren) prices after first year

P3CZ
P3C
P3C3
P3C4
P3C5

%%% With spouse prices after first year

P3SZ
P3S
P3S3
P3S4
P3S5

%%% Single prices after first year

P3EZ
P3E
P3E3
P3E4
P3E5


%%%%%%%%%%%%%%% Prices Overall = total premium = AVG COST + ADMIN % %%%%%%

%% Single Total Price each year 

AVGCOSTEZ
AVGCOSTE
AVGCOSTE2
AVGCOSTE3
AVGCOSTE4

%% Family Total Price each year 

AVGCOSTFZ
AVGCOSTF
AVGCOSTF2
AVGCOSTF3
AVGCOSTF4

%% with spouse Total Price each year 

AVGCOSTSZ
AVGCOSTS
AVGCOSTS2
AVGCOSTS3
AVGCOSTS4

%% With children Total Price each year 

AVGCOSTCZ
AVGCOSTC
AVGCOSTC2
AVGCOSTC3
AVGCOSTC4

%%%%%%%%%%%%%%% Enrollment by Family Status each year %%%%%%

% Single 

TOTALNUMEZ/400
TOTALNUME/400
TOTALNUME2/400
TOTALNUME3/400
TOTALNUME4/400

% Spouse 

TOTALNUMSZ/400
TOTALNUMS/400
TOTALNUMS2/400
TOTALNUMS3/400
TOTALNUMS4/400

% Children

TOTALNUMCZ/400
TOTALNUMC/400
TOTALNUMC2/400
TOTALNUMC3/400
TOTALNUMC4/400

% Family 

TOTALNUMFZ/400
TOTALNUMF/400
TOTALNUMF2/400
TOTALNUMF3/400
TOTALNUMF4/400

%%%%%%%%%%%%%%%%% For welfare analysis to run appropriately user has to run above code twice. First, they have to run the code for the basic 
%%%%%%%%%%%%%%%%% case where inertia is not reduced from the status quo. They should do this and save all output, especially prices and choices
%%%%%%%%%%%%%%%%% in file called load 'FULLSC CF DATA.mat' which we will reference below commented out. Then, user should run with 75% reduction 
%%%%%%%%%%%%%%%%% and save with file called ''CFSC3FOURTHS.mat'. These will both be referenced in code below and are necessary to do welfare analysis
%%%%%%%%%%%%%%%%% but ARE NOT generated above since this would just require doubling the code to this point, which is already quite lengthy. Thus, 
%%%%%%%%%%%%%%%%% users should create these two files by running above code twice and saving, in order to do actual welfare analysis similar to 
%%%%%%%%%%%%%%%%% that which was done in the paper.Code above is for 75% reduction, to run baseline set SCRED = 1 in code above and save output. 

% IF YOU RAN BASELINE WHERE SCRED = 1.0 then:

% save 'FULLSC CD DATA.mat'

% IF YOU RAN REDUCED INERTIA COUNTERFACTUAL WHERE SCRED = X (ABOVE EXAMPLE = 0.25) then:

 % save 'CFSC3FOURTHSt.mat'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% WELFARE ANALYSIS: Code to do welfare analysis. Compute change in certainty equivalent for counterfactual environment %%%%%%%%%%%%%%%%%%
%%%%%%%% model above relative to simulated environment with preferences and inertia as estimated. As noted above and in READMNE %%%%%%%%%%%%%%%%
%%%%%%%% file this analysis is for one specific counterfactual: 75% inertia reduction and inertia doesn't count for welfare analysis itself %%%%
%%%%%%%% It is easy to vary code above and here to study other counterfactuals shown in section 6 which vary these dimensions %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CLEAR and LOAD Baseline file generated from code above where SCRED = 1

clear
load 'FULLSC CD DATA.mat'

% Generate expected utilities for each plan which will be used to compute CEQ Differences. Note that we 
% exclude inertia and CDHP RC delta from calculations here to keep it simple, though in actual analysis 
% detla is included and we study range of cases where inertia is not counted to where it completely is 

% B in expeced utility final matrices stands for 'baseline' that reduced inertia 'treatment' counterfactual will be compared to 

% NOTE: Welfare change analysis only done for four years after initial choice year since coutnerfactual and actual environment are same in that year. 

for i = 1:nIs
    for k = 1:K
        EUvector12(i,k,:) = W - reshape(PPO250OP2(i,k,:),1,Sim) - reshape(P1(2,i,:),1,Sim);
        EUvector22(i,k,:) = W - reshape(PPO500OP2(i,k,:),1,Sim) - reshape(P2(2,i,:),1,Sim);
        EUvector32(i,k,:) = W - reshape(PPO1200OP2(i,k,:),1,Sim) - reshape(P3(2,i,:),1,Sim);
    end
end

EUvector12 = -(1/max(RA,0.0000001)).*exp(-((EUvector12).*max(RA,0.0000001)));
EUvector22 = -(1/max(RA,0.0000001)).*exp(-((EUvector22).*max(RA,0.0000001)));
EUvector32 = -(1/max(RA,0.0000001)).*exp(-((EUvector32).*max(RA, 0.0000001)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Now sum CARA VNM to get expected utility %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n = 1:Sim
    BEUvectorAx(:,1,n) = (1/K)*sum(reshape(EUvector12(:,:,n),nIs,K),2);
    BEUvectorAx(:,2,n) = (1/K)*sum(reshape(EUvector22(:,:,n),nIs,K),2);
    BEUvectorAx(:,3,n) = (1/K)*sum(reshape(EUvector32(:,:,n),nIs,K),2);
end

for i = 1:nIs
    for k = 1:K
        EUvector13(i,k,:) = W - reshape(PPO250OP3(i,k,:),1,Sim) - reshape(P1(3,i,:),1,Sim);
        EUvector23(i,k,:) = W - reshape(PPO500OP3(i,k,:),1,Sim) - reshape(P2(3,i,:),1,Sim);
        EUvector33(i,k,:) = W - reshape(PPO1200OP3(i,k,:),1,Sim) - reshape(P3(3,i,:),1,Sim);
    end
end

EUvector13 = -(1/max(RA,0.0000001)).*exp(-((EUvector13).*max(RA,0.0000001)));
EUvector23 = -(1/max(RA,0.0000001)).*exp(-((EUvector23).*max(RA,0.0000001)));
EUvector33 = -(1/max(RA,0.0000001)).*exp(-((EUvector33).*max(RA, 0.0000001)));

for n = 1:Sim
    BEUvectorBx(:,1,n) = (1/K)*sum(reshape(EUvector13(:,:,n),nIs,K),2);
    BEUvectorBx(:,2,n) = (1/K)*sum(reshape(EUvector23(:,:,n),nIs,K),2);
    BEUvectorBx(:,3,n) = (1/K)*sum(reshape(EUvector33(:,:,n),nIs,K),2);
end

for i = 1:nIs
    for k = 1:K
        EUvector14(i,k,:) = W - reshape(PPO250OP3(i,k,:),1,Sim) - reshape(P1(4,i,:),1,Sim);
        EUvector24(i,k,:) = W - reshape(PPO500OP3(i,k,:),1,Sim) - reshape(P2(4,i,:),1,Sim);
        EUvector34(i,k,:) = W - reshape(PPO1200OP3(i,k,:),1,Sim) - reshape(P3(4,i,:),1,Sim);
    end
end

EUvector14 = -(1/max(RA,0.0000001)).*exp(-((EUvector14).*max(RA,0.0000001)));
EUvector24 = -(1/max(RA,0.0000001)).*exp(-((EUvector24).*max(RA,0.0000001)));
EUvector34 = -(1/max(RA,0.0000001)).*exp(-((EUvector34).*max(RA, 0.0000001)));

%%%%%%%%%%%%%%%%%%% Now sum CARA VNM to get expected utility %%%%%%%%%%%%

for n = 1:Sim
    BEUvectorCx(:,1,n) = (1/K)*sum(reshape(EUvector14(:,:,n),nIs,K),2);
    BEUvectorCx(:,2,n) = (1/K)*sum(reshape(EUvector24(:,:,n),nIs,K),2);
    BEUvectorCx(:,3,n) = (1/K)*sum(reshape(EUvector34(:,:,n),nIs,K),2);
end

for i = 1:nIs
    for k = 1:K
        EUvector15(i,k,:) = W - reshape(PPO250OP3(i,k,:),1,Sim) - reshape(P1(5,i,:),1,Sim);
        EUvector25(i,k,:) = W - reshape(PPO500OP3(i,k,:),1,Sim) - reshape(P2(5,i,:),1,Sim);
        EUvector35(i,k,:) = W - reshape(PPO1200OP3(i,k,:),1,Sim) - reshape(P3(5,i,:),1,Sim);
    end
end

EUvector15 = -(1/max(RA,0.0000001)).*exp(-((EUvector15).*max(RA,0.0000001)));
EUvector25 = -(1/max(RA,0.0000001)).*exp(-((EUvector25).*max(RA,0.0000001)));
EUvector35 = -(1/max(RA,0.0000001)).*exp(-((EUvector35).*max(RA, 0.0000001)));

%%%%%%%%%%%%%%%%%%% Now sum CARA VNM to get expected utility %%%%%%%%%%%%

for n = 1:Sim
    BEUvectorDx(:,1,n) = (1/K)*sum(reshape(EUvector15(:,:,n),nIs,K),2);
    BEUvectorDx(:,2,n) = (1/K)*sum(reshape(EUvector25(:,:,n),nIs,K),2);
    BEUvectorDx(:,3,n) = (1/K)*sum(reshape(EUvector35(:,:,n),nIs,K),2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Define prices from baseline simulation for each year   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BP1 = zeros(nIs,3);
BP2 = zeros(nIs,3);
BP3 = zeros(nIs,3);
BP4 = zeros(nIs,3);

for i = 1:nIs
   BP1(i,1) = P1(2,i,1);
   BP1(i,2) = P2(2,i,1);
   BP1(i,3) = P3(2,i,1);

   BP2(i,1) = P1(3,i,1);
   BP2(i,2) = P2(3,i,1);
   BP2(i,3) = P3(3,i,1);

   BP3(i,1) = P1(4,i,1);
   BP3(i,2) = P2(4,i,1);
   BP3(i,3) = P3(4,i,1);
   
   BP4(i,1) = P1(5,i,1);
   BP4(i,2) = P2(5,i,1);
   BP4(i,3) = P3(5,i,1);
end

%%%%%%%%%%%%%%%% Store choices made from baseline simulation

CCX1 = CC1;
CCX2 = CC2;
CCX3 = CC3;
CCX4 = CC4;

% clear key variables from baseline file that we'll need to replace with reduced ienrtia values to compute utilities in that scenario and add into same database

clear CC1 CC2 CC3 CC4 P1 P2 P3

%%%%%%%%%%% NOW, load in specific variables from file that saved code run above for 75% reduction in switching cost. 

load 'CFSC3FOURTHSt.mat' CC1 CC2 CC3 CC4 P1 P2 P3

% Store treatment prices 

TP1 = zeros(nIs,3);
TP2 = zeros(nIs,3);
TP3 = zeros(nIs,3);
TP4 = zeros(nIs,3);

for i = 1:nIs
   TP1(i,1) = P1(2,i,1);
   TP1(i,2) = P2(2,i,1);
   TP1(i,3) = P3(2,i,1);

   TP2(i,1) = P1(3,i,1);
   TP2(i,2) = P2(3,i,1);
   TP2(i,3) = P3(3,i,1);

   TP3(i,1) = P1(4,i,1);
   TP3(i,2) = P2(4,i,1);
   TP3(i,3) = P3(4,i,1);
   
   TP4(i,1) = P1(5,i,1);
   TP4(i,2) = P2(5,i,1);
   TP4(i,3) = P3(5,i,1);
end

%%%%%%%%%%%%% Now compute expected utilities for scenario with reduced inertia

for i = 1:nIs
    for k = 1:K
        EUvector12(i,k,:) = W - reshape(PPO250OP2(i,k,:),1,Sim) - reshape(P1(2,i,:),1,Sim);
        EUvector22(i,k,:) = W - reshape(PPO500OP2(i,k,:),1,Sim) - reshape(P2(2,i,:),1,Sim);
        EUvector32(i,k,:) = W - reshape(PPO1200OP2(i,k,:),1,Sim) - reshape(P3(2,i,:),1,Sim);
    end
end

EUvector12 = -(1/max(RA,0.0000001)).*exp(-((EUvector12).*max(RA,0.0000001)));
EUvector22 = -(1/max(RA,0.0000001)).*exp(-((EUvector22).*max(RA,0.0000001)));
EUvector32 = -(1/max(RA,0.0000001)).*exp(-((EUvector32).*max(RA, 0.0000001)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Now sum CARA VNM to get expected utility %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n = 1:Sim
    TEUvectorAx(:,1,n) = (1/K)*sum(reshape(EUvector12(:,:,n),nIs,K),2);
    TEUvectorAx(:,2,n) = (1/K)*sum(reshape(EUvector22(:,:,n),nIs,K),2);
    TEUvectorAx(:,3,n) = (1/K)*sum(reshape(EUvector32(:,:,n),nIs,K),2);
end

for i = 1:nIs
    for k = 1:K
        EUvector13(i,k,:) = W - reshape(PPO250OP3(i,k,:),1,Sim) - reshape(P1(3,i,:),1,Sim);
        EUvector23(i,k,:) = W - reshape(PPO500OP3(i,k,:),1,Sim) - reshape(P2(3,i,:),1,Sim);
        EUvector33(i,k,:) = W - reshape(PPO1200OP3(i,k,:),1,Sim) - reshape(P3(3,i,:),1,Sim);
    end
end

EUvector13 = -(1/max(RA,0.0000001)).*exp(-((EUvector13).*max(RA,0.0000001)));
EUvector23 = -(1/max(RA,0.0000001)).*exp(-((EUvector23).*max(RA,0.0000001)));
EUvector33 = -(1/max(RA,0.0000001)).*exp(-((EUvector33).*max(RA, 0.0000001)));

for n = 1:Sim
    TEUvectorBx(:,1,n) = (1/K)*sum(reshape(EUvector13(:,:,n),nIs,K),2);
    TEUvectorBx(:,2,n) = (1/K)*sum(reshape(EUvector23(:,:,n),nIs,K),2);
    TEUvectorBx(:,3,n) = (1/K)*sum(reshape(EUvector33(:,:,n),nIs,K),2);
end

for i = 1:nIs
    for k = 1:K
        EUvector14(i,k,:) = W - reshape(PPO250OP3(i,k,:),1,Sim) - reshape(P1(4,i,:),1,Sim);
        EUvector24(i,k,:) = W - reshape(PPO500OP3(i,k,:),1,Sim) - reshape(P2(4,i,:),1,Sim);
        EUvector34(i,k,:) = W - reshape(PPO1200OP3(i,k,:),1,Sim) - reshape(P3(4,i,:),1,Sim);
    end
end

EUvector14 = -(1/max(RA,0.0000001)).*exp(-((EUvector14).*max(RA,0.0000001)));
EUvector24 = -(1/max(RA,0.0000001)).*exp(-((EUvector24).*max(RA,0.0000001)));
EUvector34 = -(1/max(RA,0.0000001)).*exp(-((EUvector34).*max(RA, 0.0000001)));

%%%%%%%%%%%%%%%%%%% Now sum CARA VNM to get expected utility %%%%%%%%%%%%

for n = 1:Sim
    TEUvectorCx(:,1,n) = (1/K)*sum(reshape(EUvector14(:,:,n),nIs,K),2);
    TEUvectorCx(:,2,n) = (1/K)*sum(reshape(EUvector24(:,:,n),nIs,K),2);
    TEUvectorCx(:,3,n) = (1/K)*sum(reshape(EUvector34(:,:,n),nIs,K),2);
end

for i = 1:nIs
    for k = 1:K
        EUvector15(i,k,:) = W - reshape(PPO250OP3(i,k,:),1,Sim) - reshape(P1(5,i,:),1,Sim);
        EUvector25(i,k,:) = W - reshape(PPO500OP3(i,k,:),1,Sim) - reshape(P2(5,i,:),1,Sim);
        EUvector35(i,k,:) = W - reshape(PPO1200OP3(i,k,:),1,Sim) - reshape(P3(5,i,:),1,Sim);
    end
end

EUvector15 = -(1/max(RA,0.0000001)).*exp(-((EUvector15).*max(RA,0.0000001)));
EUvector25 = -(1/max(RA,0.0000001)).*exp(-((EUvector25).*max(RA,0.0000001)));
EUvector35 = -(1/max(RA,0.0000001)).*exp(-((EUvector35).*max(RA, 0.0000001)));

%%%%%%%%%%%%%%%%%%% Now sum CARA VNM to get expected utility %%%%%%%%%%%%

for n = 1:Sim
    TEUvectorDx(:,1,n) = (1/K)*sum(reshape(EUvector15(:,:,n),nIs,K),2);
    TEUvectorDx(:,2,n) = (1/K)*sum(reshape(EUvector25(:,:,n),nIs,K),2);
    TEUvectorDx(:,3,n) = (1/K)*sum(reshape(EUvector35(:,:,n),nIs,K),2);
end

%%%% Now save choices from counterfactual environment with reduced inertia 

CCZ1 = CC1;
CCZ2 = CC2;
CCZ3 = CC3;
CCZ4 = CC4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Compute all CEQs for each choice in treatment environment and each choice in baseline environment  %%%%%%%%
%%%%%%%%%%%%%%% These will then be compared to one another in the welfare calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Invert CARA expected utility formula to solve for certainty equivalent as function of total expected utility
%%%% Email the author for a file that shows this derivation, formula given below for computation of CEQ here. 

BCEQ = zeros(nIs,3,Sim);
BCEQ2 = zeros(nIs,3,Sim);
BCEQ3 = zeros(nIs,3,Sim);
BCEQ4 = zeros(nIs,3,Sim);

TCEQ = zeros(nIs,3,Sim);
TCEQ2 = zeros(nIs,3,Sim);
TCEQ3 = zeros(nIs,3,Sim);
TCEQ4 = zeros(nIs,3,Sim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Compute certainty equivalent using formula solved for by inverting CARA function %%%%%%%%%
%%%%%% and setting CEQ from CARA equal to expected utility from enrolling in each plan %%%%%%%%%%
%%%%%% in each counterfactual. Actual analysis used another longer / less-efficient method %%%%%%
%%%%%% to find the same CEQ so this is a more concise way to do this. NOTE: in trials with %%%%%%
%%%%%% this new code with simulated data, runtime for next loop took good amount of time  %%%%%%%
%%%%%% close to one hour, so users may want to experiment with other methods for computing %%%%%%
%%%%%% CEQs / adjusting this new more concise code to run faster. For method used in actual %%%%%
%%%%%% analysis, which leads to same results but is longer code wise, contact the author via %%%%
%%%%%% email at handel@berkeley.edu %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% CEQ are total $ CEQ, wealth assumed 75000
%%% CEQ loss 'at stake' is 75000-CEQ value

RA = max(RA,0.0000001);

BCEQ(:,1,:) = (-log(-RA(:,1,:).*BEUvectorAx(:,1,:)))./(RA(:,1,:));
BCEQ(:,2,:) = (-log(-RA(:,1,:).*BEUvectorAx(:,2,:)))./(RA(:,1,:));
BCEQ(:,3,:) = (-log(-RA(:,1,:).*BEUvectorAx(:,3,:)))./(RA(:,1,:));

BCEQ2(:,1,:) = (-log(-RA(:,1,:).*BEUvectorBx(:,1,:)))./(RA(:,1,:));
BCEQ2(:,2,:) = (-log(-RA(:,1,:).*BEUvectorBx(:,2,:)))./(RA(:,1,:));
BCEQ2(:,3,:) = (-log(-RA(:,1,:).*BEUvectorBx(:,3,:)))./(RA(:,1,:));

BCEQ3(:,1,:) = (-log(-RA(:,1,:).*BEUvectorCx(:,1,:)))./(RA(:,1,:));
BCEQ3(:,2,:) = (-log(-RA(:,1,:).*BEUvectorCx(:,2,:)))./(RA(:,1,:));
BCEQ3(:,3,:) = (-log(-RA(:,1,:).*BEUvectorCx(:,3,:)))./(RA(:,1,:));

BCEQ4(:,1,:) = (-log(-RA(:,1,:).*BEUvectorDx(:,1,:)))./(RA(:,1,:));
BCEQ4(:,2,:) = (-log(-RA(:,1,:).*BEUvectorDx(:,2,:)))./(RA(:,1,:));
BCEQ4(:,3,:) = (-log(-RA(:,1,:).*BEUvectorDx(:,3,:)))./(RA(:,1,:));

TCEQ(:,1,:) = (-log(-RA(:,1,:).*TEUvectorAx(:,1,:)))./(RA(:,1,:));
TCEQ(:,2,:) = (-log(-RA(:,1,:).*TEUvectorAx(:,2,:)))./(RA(:,1,:));
TCEQ(:,3,:) = (-log(-RA(:,1,:).*TEUvectorAx(:,3,:)))./(RA(:,1,:));

TCEQ2(:,1,:) = (-log(-RA(:,1,:).*TEUvectorBx(:,1,:)))./(RA(:,1,:));
TCEQ2(:,2,:) = (-log(-RA(:,1,:).*TEUvectorBx(:,2,:)))./(RA(:,1,:));
TCEQ2(:,3,:) = (-log(-RA(:,1,:).*TEUvectorBx(:,3,:)))./(RA(:,1,:));

TCEQ3(:,1,:) = (-log(-RA(:,1,:).*TEUvectorCx(:,1,:)))./(RA(:,1,:));
TCEQ3(:,2,:) = (-log(-RA(:,1,:).*TEUvectorCx(:,2,:)))./(RA(:,1,:));
TCEQ3(:,3,:) = (-log(-RA(:,1,:).*TEUvectorCx(:,3,:)))./(RA(:,1,:));

TCEQ4(:,1,:) = (-log(-RA(:,1,:).*TEUvectorDx(:,1,:)))./(RA(:,1,:));
TCEQ4(:,2,:) = (-log(-RA(:,1,:).*TEUvectorDx(:,2,:)))./(RA(:,1,:));
TCEQ4(:,3,:) = (-log(-RA(:,1,:).*TEUvectorDx(:,3,:)))./(RA(:,1,:));
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Now Calculate Welfare Change in 4 years resulting from inertia reduction  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ExpectChange1W = zeros(nIs,1);
ExpectChange2W = zeros(nIs,1);
ExpectChange3W = zeros(nIs,1);
ExpectChange4W = zeros(nIs,1);

NUM1W = zeros(nIs,1);
NUM2W = zeros(nIs,1);
NUM3W = zeros(nIs,1);
NUM4W = zeros(nIs,1);

%%%%%%%%%%%%%%%%%% Calculating welfare change at the individual level
%%%%%%%%%%%%%%%%%% Between baseline CEQ from enrollment in baseline environment with inertia 
%%%%%%%%%%%%%%%%%% to CEQ in counterfactual environment with reduced inertia for plan actually chosen. 
%%%%%%%%%%%%%%%%%% Change in plans chosen from baseline to counterfactual AND change in plan premiums 
%%%%%%%%%%%%%%%%%% are both responsible for consumer welfare change from reduced inertia in counterfactual 
%%%%%%%%%%%%%%%%%% environment. 

for i = 1:nIs
    for s=1:Sim
            if (AcceptRejectF(i,s)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject2(i,s)==1))
             ExpectChange1W(i,1) = ExpectChange1W(i,1) + TCEQ(i,CCX1(i,s),s)-BCEQ(i,CCZ1(i,s),s);  
             NUM1W(i,1) = NUM1W(i,1)+1;
            end
            if (AcceptRejectF(i,s)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,s)==1))
             ExpectChange2W(i,1) = ExpectChange2W(i,1) + TCEQ2(i,CCX2(i,s),s)-BCEQ2(i,CCZ2(i,s),s);  
             NUM2W(i,1) = NUM2W(i,1)+1;
            end 
			if (AcceptRejectF(i,s)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,s)==1))
             ExpectChange3W(i,1) = ExpectChange3W(i,1) + TCEQ3(i,CCX3(i,s),s)-BCEQ3(i,CCZ3(i,s),s);  
             NUM3W(i,1) = NUM3W(i,1)+1;
             end
             if (AcceptRejectF(i,s)==1 | (mean(AcceptRejectF(i,:),2)==0 & AcceptReject3(i,s)==1))
             ExpectChange4W(i,1) = ExpectChange4W(i,1) + TCEQ4(i,CCX4(i,s),s)-BCEQ4(i,CCZ4(i,s),s);  
             NUM4W(i,1) = NUM4W(i,1)+1;
             end 
    end 
end

%%%%%% Computing individual level welfare gain / loss from treatment averaging across simulations of random coefficients 

for i = 1:nIs
   ExpectChange1W(i,1) = ExpectChange1W(i,1)/NUM1W(i,1);
   ExpectChange2W(i,1) = ExpectChange2W(i,1)/NUM2W(i,1);
   ExpectChange3W(i,1) = ExpectChange3W(i,1)/NUM3W(i,1);
   ExpectChange4W(i,1) = ExpectChange4W(i,1)/NUM4W(i,1);
end

%%%%%%%% Now begin computation of average change for whole population 
    
GAIN1W = zeros(1,1);
GAIN2W = zeros(1,1);
GAIN3W = zeros(1,1);
GAIN4W = zeros(1,1);

N1W = zeros(1,1);
N2W = zeros(1,1);
N3W = zeros(1,1);
N4W = zeros(1,1);

for i = 1:nIs
        GAIN1W  = GAIN1W + ExpectChange1W(i,1);
        N1W = N1W + 1;
        GAIN2W  = GAIN2W + ExpectChange2W(i,1);
        N2W = N2W + 1;
		GAIN3W  = GAIN3W + ExpectChange3W(i,1);
        N3W = N3W + 1;
        GAIN4W  = GAIN4W + ExpectChange4W(i,1);
        N4W = N4W + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Final per person per year change in CEQ resulting from reduction in inertia

GAIN1W = GAIN1W/N1W;
GAIN2W = GAIN2W/N2W;
GAIN3W = GAIN3W/N3W;
GAIN4W = GAIN4W/N4W;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% END OF WELFARE CODE PROVIDED AS EXAMPLE HERE %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADDITIONAL KEY NOTES TO TRANSLATE EXAMPLE WELFARE CALCULATIONS ABOVE TO OTHER WELFARE CALCULATIONS DONE IN PAPER 

% 1. Above analysis is consumer welfare change per person per year overall i terms of $ certainty equivalent change 
% 2. For % change as done in paper, user should sum up either total premiums at stake in baseline, total premiums plus
%    out of pocket expenses in baseline, or total CEQ losses from Wealth (75000) at stake in baseline (which will be bigger than mean OOP + premiums)
%    See the welfare section of the paper for further discussion of these benchmarks. Computing them shuold be straightforward using 
%    simulated data and code similar to that provided here. 
% 3. Section 6 discusses distributional implications of reduced inertia by conditioning welfare analysis on certain demographics / groups of consumers
%    to do this kind of analysis, just do welfare analysis above but condition on additional information specific to the calculation (e.g. high income) 
% 4. To study counterfactuals with different levels of inertia reduction just save two different files resulting from two different runs of primary counterfactual
%    code supplied before welfare calculation above, save the baseline and the new counterfactual, and re-run the analysis above. 
% 5. To calculate first-best in this environment, and changes from first-best (Table 7 in paper) just compute welfare / CEQ for people all in PPO250 and compare to 
%    specific counterfactual in question (or compare to CEQs realized in observed environment)
% 6. To copmute coutnerfactuals where inertia matters in welfare calculation, just compute total inertial / switching costs $ incurred in counterfactual vs. 
%    baseline. Then, if you want inertia to count as switching cost that fully counts toward welfare analysis, take change in switching costs incurred from 
%    treatment to control and add to change in CEQ given above. If you want X% of inertia to count for welfare, add X% of this difference to welfare CEQ change 
%    given above. See section 6 in paper and table 8 for additional analysis of this point. 
% 7. Note: in actual paper, CEQ is computed using delta / CEQ random coefficient, but including it or not makes minimal difference sicne there is very little 
%   switching back and forth from this plan in the data because of the low mean good size variance of this random coefficient

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%













