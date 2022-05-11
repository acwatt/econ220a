function output = Likelihood(alpha,R,EPS,nIs,nPlans,nTs,HTC,Sim,K,EUvector11,EUvector21,EUvector31,EUvector12,EUvector22,EUvector32,EUvector13,EUvector23,EUvector33,choice,choice12,choice22,choice32,choice13,choice23,choice33,Income,Tier2,IND,FSAY2,FSAY3,QS,managerX,CC1,CC2,CSAL2,CSAL3,MAge)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Likelihood function for simulated maximum likelihood choice model estimation 		%%%%
%%%%%%%%% called in EstimationCode.m file. NOTE: the material here is described in detail	%%%%
%%%%%%%%% also in Online Appendix B in the online materials posted for this paper. 			%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Broadly speaking, the likelihood function sets up the sequence of expected utility %%%
%%%%%%%%% decisions as a function of the candidate parameters and underlying observed 		 %%%
%%%%%%%%% demographic and health risk data. Then, it matches choices predicted under the 	%%%%
%%%%%%%%% candidate parameters to the actual choices made, which determines the 'likelihood' %%%
%%%%%%%%% of observing actual choice sequences. Then, optimizer picks parameters that best 	 %%%
%%%%%%%%% match the simulated choices to the actual choices made							 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%																					%%%%
%%%%%%%%%		Ben Handel																	%%%%
%%%%%%%%%		handel@berkeley.edu															%%%%
%%%%%%%%%																					%%%%
%%%%%%%%%		ADVERSE SELECTION AND INERTIA IN HEALTH INSURANCE MARKETS:					%%%%
%%%%%%%%%		WHEN NUDGING HURTS															%%%%
%%%%%%%%%																					%%%%
%%%%%%%%%		March 2013																	%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%% NOTE: Candidate parameters in alpha are defined after non-linear optimization 
%%%%%%%%% is called in choide model estimation file. That should be referenced if there is
%%%%%%%%% uncertainty about what given alpha(x,y) is representing. 

%%%%%%%%% Set up CDHP preference coefficients as function of underlying parameters

CDHPmean = alpha(2,3);
CDHPmeanE = alpha(2,1);

CDHPmrc = ((CDHPmean')*(ones(1,Sim)))';
CDHPmrcE = ((CDHPmeanE')*(ones(1,Sim)))';

VRC = R(:,1)*alpha(2,4);
VRCE = R(:,1)*alpha(2,2);

CDHPt = zeros(nPlans,Sim);
CDHPtE = zeros(nPlans,Sim);
CDHPt(3,:)= (CDHPmrc+VRC(:,1))';
CDHPtE(3,:)=(CDHPmrcE+VRCE(:,1))';

%%%%%%%%% Set up risk preference coefficients as function of underlying parameters
%%%%%%%%% NOTE: For some non-linear esimation routines, dividing the risk preference coefficient
%%%%%%%%% by 100 or another positive factors helps the routine find an optimum because otherwise 
%%%%%%%%% the CARA coefficients are very small numbers. This normalization IS reflected as well 
%%%%%%%%% in the starting values in EsimationCode_2011_1284.m. The correct / un-normalized 
%%%%%%%%% CARA coefficients are the ones ultimately used / reported. 

RAmean = alpha(1,1)/100;
RAincome = alpha(1,2)/100;
RAage = alpha(1,3)/100;

RArct = RAmean*ones(1,nIs) + RAincome*(Income') + RAage*(MAge');

RArct = (RArct'*ones(1,Sim))';
RAvar = (R(:,2)*(alpha(1,4)/100))*ones(1,nIs);
RArct = (RArct + RAvar)';

%%%%%%%%% Set up PPO250 HTC preferences as function of underlying parameters

CHTC = alpha(4,1)*HTC;

%%%%%%%%% Set up Epsilon and CDHP matrices to fill in as function of individual vs. family

CDHPti = zeros(nIs,Sim);
Epsilon1200t1 = zeros(nIs,Sim);
Epsilon500t1 = zeros(nIs,Sim);
Epsilon1200t2 = zeros(nIs,Sim);
Epsilon500t2 = zeros(nIs,Sim);
Epsilon1200t3 = zeros(nIs,Sim);
Epsilon500t3 = zeros(nIs,Sim);
dummy = ones(nIs,1);

finder0 = find(IND'==0);
CDHPti(finder0,:) = dummy(finder0,:)*CDHPt(3,:);
Epsilon1200t1(finder0,:) = dummy(finder0,:)*(EPS(:,1)*alpha(3,4))';
Epsilon500t1(finder0,:) = dummy(finder0,:)*(EPS(:,2)*alpha(3,3))';
Epsilon1200t2(finder0,:) = dummy(finder0,:)*(EPS(:,3)*alpha(3,4))';
Epsilon500t2(finder0,:) = dummy(finder0,:)*(EPS(:,4)*alpha(3,3))';
Epsilon1200t3(finder0,:) = dummy(finder0,:)*(EPS(:,5)*alpha(3,4))';
Epsilon500t3(finder0,:) = dummy(finder0,:)*(EPS(:,6)*alpha(3,3))';

finder1 = find(IND'==1);
CDHPti(finder1,:) = dummy(finder1,:)*CDHPtE(3,:);
Epsilon1200t1(finder1,:) = dummy(finder1,:)*(EPS(:,1)*alpha(3,2))';
Epsilon500t1(finder1,:) = dummy(finder1,:)*(EPS(:,2)*alpha(3,1))';
Epsilon1200t2(finder1,:) = dummy(finder1,:)*(EPS(:,3)*alpha(3,2))';
Epsilon500t2(finder1,:) = dummy(finder1,:)*(EPS(:,4)*alpha(3,1))';
Epsilon1200t3(finder1,:) = dummy(finder1,:)*(EPS(:,5)*alpha(3,2))';
Epsilon500t3(finder1,:) = dummy(finder1,:)*(EPS(:,6)*alpha(3,1))';

%%%%%%%%%%%%%%%%%%%% Fill in switching costs as function of underlying fundamentals / candidate parameters

SCE = alpha(4,2)'*ones(nIs,Sim);
SC2 = ones(nIs,Sim)*alpha(4,3);
SC3 = ones(nIs,Sim)*alpha(4,3);

SC2(IND'==1,:) = SCE(IND'==1,:);
SC3(IND'==1,:) = SCE(IND'==1,:);

SC2 = SC2 + (FSAY2*ones(1,Sim)*alpha(4,4)) + (Income*alpha(5,1)*ones(1,Sim))+ (QS*ones(1,Sim)*alpha(5,2)) + (managerX*ones(1,Sim)*alpha(5,3)) + (CC1*ones(1,Sim)*alpha(5,4)) + (CSAL2*ones(1,Sim)*alpha(6,1));
SC3 = SC3 + (FSAY3*ones(1,Sim)*alpha(4,4)) + (Income*alpha(5,1)*ones(1,Sim))+ (QS*ones(1,Sim)*alpha(5,2)) + (managerX*ones(1,Sim)*alpha(5,3)) + (CC2*ones(1,Sim)*alpha(5,4)) + (CSAL3*ones(1,Sim)*alpha(6,1));

%%%%%%%%%%%%%%%%%%%%%% Set up 3-D matrices to speed up likelihood computation. Add dimension K of health draws to expedite things with faster matrix operations. 

CDHP = zeros(nIs,K,Sim);

SwC12 = zeros(nIs,K,Sim);
SwC22 = zeros(nIs,K,Sim);
SwC32 = zeros(nIs,K,Sim);
SwC13 = zeros(nIs,K,Sim);
SwC23 = zeros(nIs,K,Sim);
SwC33 = zeros(nIs,K,Sim);

Epsilon12001 = zeros(nIs,K,Sim);
Epsilon5001 = zeros(nIs,K,Sim);
Epsilon12002 = zeros(nIs,K,Sim);
Epsilon5002 = zeros(nIs,K,Sim);
Epsilon12003 = zeros(nIs,K,Sim);
Epsilon5003 = zeros(nIs,K,Sim);

RA = zeros(nIs,K,Sim);

for k =1:K
RA(:,k,:) = RArct;
CDHP(:,k,:) = CDHPti;
SwC12(:,k,:) = SC2.*choice12;
SwC22(:,k,:) = SC2.*choice22;
SwC32(:,k,:) = SC2.*choice32;
SwC13(:,k,:) = SC3.*choice13;
SwC23(:,k,:) = SC3.*choice23;
SwC33(:,k,:) = SC3.*choice33;
Epsilon12001(:,k,:) = Epsilon1200t1;
Epsilon5001(:,k,:) = Epsilon500t1;
Epsilon12002(:,k,:) = Epsilon1200t2;
Epsilon5002(:,k,:) = Epsilon500t2;
Epsilon12003(:,k,:) = Epsilon1200t3;
Epsilon5003(:,k,:) = Epsilon500t3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute expected utillity index by health state for given individual and risk 
% preference simulation draw. Value of utility calculated is negative of true CARA
% preference values, which is accounted for throughout likelihood function. If there 
% is any confusion in the calculation don't hesitate to contact Ben Handel. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RA = max(RA, 0.0000001);
EUvector11 = exp(-((EUvector11).*RA)+300);
EUvector21 = exp(-((EUvector21 + CHTC + Epsilon5001).*RA)+300);
EUvector31 = exp(-((EUvector31 + CDHP + CHTC + Epsilon12001).*RA)+300);
EUvector12 = exp(-((EUvector12 + SwC12).*RA)+300);
EUvector22 = exp(-((EUvector22 + SwC22 + CHTC + Epsilon5002).*RA)+300);
EUvector32 = exp(-((EUvector32 + CDHP + CHTC + SwC32 + Epsilon12002).*RA)+300);
EUvector13 = exp(-((EUvector13 + SwC13).*RA)+300);
EUvector23 = exp(-((EUvector23 + SwC23 + CHTC + Epsilon5003).*RA)+300);
EUvector33 = exp(-((EUvector33 + CDHP + CHTC + SwC33 + Epsilon12003).*RA)+300);

EUvectorA = zeros(nIs,nPlans,Sim);
EUvectorB = zeros(nIs,nPlans,Sim);
EUvectorC = zeros(nIs,nPlans,Sim);

%%%%%%%%%%%%% Compute expected utility for each potential choice for each family in each year (and for each risk preference simulation)

EUvectorA(:,1,:) = mean(EUvector11,2);
EUvectorA(:,2,:) = mean(EUvector21,2);
EUvectorA(:,3,:) = mean(EUvector31,2);
EUvectorB(:,1,:) = mean(EUvector12,2);
EUvectorB(:,2,:) = mean(EUvector22,2);
EUvectorB(:,3,:) = mean(EUvector32,2);
EUvectorC(:,1,:) = mean(EUvector13,2);
EUvectorC(:,2,:) = mean(EUvector23,2);
EUvectorC(:,3,:) = mean(EUvector33,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Finding Simulated Choice Probabilities %%%%%%%%%%%%%%%%%%%%
%%%%%%%% As function of candidate parameters %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Simulated maximum likelihood with %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Smoothed-Accept-Reject Methodology %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% See estimation appendix for further details %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = zeros(nIs,1);
PR = ones(nIs,Sim);
PT = ones(nIs,nTs,nPlans,Sim);
PTAgg = ones(nIs,nTs,1,Sim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%% Define utility relative to PPO1200 to keep things simple 
%%%%%%%%%%%%% Relative utility is PPO1200 utility / PPO X utility, because actual 
%%%%%%%%%%%%% utility calculated above is negative of CARA, so lower values of that function 
%%%%%%%%%%%%% are better (since they would be higher negative numbers)

PT = ones(nIs,nTs,nPlans,Sim);

%Expanding just the PPO1200 utility to a higher dimension so we can do the arithmetic
EUvectorPPO1200 = zeros(nIs,nPlans,Sim,nTs);
for plan = 1:nPlans
	EUvectorPPO1200(:,plan,:,1) = EUvectorA(:,3,:);
	EUvectorPPO1200(:,plan,:,2) = EUvectorB(:,3,:);
	EUvectorPPO1200(:,plan,:,3) = EUvectorC(:,3,:);
end
   

PT(:,1,1:nPlans-1,:) = EUvectorPPO1200(:,1:nPlans-1,:,1) ./ EUvectorA(:,1:nPlans-1,:);
PT(:,2,1:nPlans-1,:) = EUvectorPPO1200(:,1:nPlans-1,:,2) ./ EUvectorB(:,1:nPlans-1,:);
PT(:,3,1:nPlans-1,:) = EUvectorPPO1200(:,1:nPlans-1,:,3) ./ EUvectorC(:,1:nPlans-1,:);
  

%%%%%%%%%% Smoothed Accept-Reject: Smoothing gets rid of discreteness issue with accept-reject estimation
%%%%%%%%%% Smoothing takes 0-1 binary outcome for sequence of choices given utilties, and makes it smooth 
%%%%%%%%%% function of underlying difference in utilities. Here, we raise to power 6 which has effect of 
%%%%%%%%%% having probability of choosing maximum utility sequence approach 1. So, this is a way to have code 
%%%%%%%%%% work smoothly but approximate standard accept-reject with simulated MLE. 
%%%%%%%%%% For more details see Online Appendix B on choice model estimation. 

PT = PT.^6; %Smoothing factor
PTAgg = sum(PT,3);

for plan = 1:nPlans
	PT(:,:,plan,:) = PT(:,:,plan,:) ./ PTAgg;
end


%%%%%%%%%% Compute likelihood of actual sequence of choices made in data, given underlying parameters. This will go into likelihood
%%%%%%%%%% function which is then optimized by ktrlink in choice model estimation code. 
	
for i = 1:nIs
	PR(i,:) = PT(i,1,choice(i,1),:).*PT(i,2,choice(i,2),:).*PT(i,3,choice(i,3),:);
end

P = mean(PR,2);
result = -sum(log(P));

%%%%%%%%%%%%%%%%%%%%%%%% Negative of Log-Likelihood Function since optimizer does minimization 
%%%%%%%%%%%%%%%%%%%%%%%% Built in part to 'keep searching' if function returns Nan

if isnan(result) || isinf(result)
    output = 10^20;
else
    output = result;
end