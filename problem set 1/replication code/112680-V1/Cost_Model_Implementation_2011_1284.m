%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   Cost Model Matlab Implementation %%%%%%%%%%%%
%%%%%%%%%%%%  (PART II - Cost Estimation done   %%%%%%%%%%%%
%%%%%%%%%%%%   in part I in STATA)              %%%%%%%%%%%%
%%%%%%%%%%%%   			                        %%%%%%%%%%%%
%%%%%%%%%%%%   Ben Handel                       %%%%%%%%%%%%
%%%%%%%%%%%%									%%%%%%%%%%%%
%%%%%%%%%%%%   Adverse Selection and Inertia    %%%%%%%%%%%%
%%%%%%%%%%%%   in Health Insurane Markets:      %%%%%%%%%%%%
%%%%%%%%%%%%   When Nudging Hurts               %%%%%%%%%%%%
%%%%%%%%%%%%                                    %%%%%%%%%%%%
%%%%%%%%%%%%   February, 2013                   %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% This file describes the code the translates cost 
%%%%%%%%%%% model estimates from the 'Cost_Model_Estimation'
%%%%%%%%%%% file, which estiamtes the cost model in STATA, 
%%%%%%%%%%% into individual simulations of heatlh risk 
%%%%%%%%%%% distributions in MATLAB. Note that because the 
%%%%%%%%%%% data are proprietary, the data inputs used  
%%%%%%%%%%% are approxiamte representations of the data actually
%%%%%%%%%%% used in the paper. 

load 'ASIN-Handel-SimulatedData_2011_1284.mat'

%%%%%%%%%%% The first part of this file contains code for 
%%%%%%%%%%% Simulating individual health risk distributions
%%%%%%%%%%% from the cost model estimates contained in the 
%%%%%%%%%%% ASIN-Handel-SimulatdData.mat file.

%%%%%%%%%%% This replicates the material in the Online Appendix 
%%%%%%%%%%% On Cost Model Estimation 

%%%%%%%%%%% As noted in the appendix, health risk distributions
%%%%%%%%%%% are first constructed at the individual level, then, 
%%%%%%%%%%% ultimately combined to the family level later in this
%%%%%%%%%%% file to determine plan specific out-of-pocket risk 
%%%%%%%%%%% for each family. 

%%%%%%%%%%% We simulate expenditure risk distributions at the individual
%%%%%%%%%%% level for four categories of medical expenditures: 
%%%%%%%%%%% (i) Mental Health
%%%%%%%%%%% (ii) Pharmacy
%%%%%%%%%%% (iii) Office Visit
%%%%%%%%%%% (iv) All Other (e.g. Hospital) 

%%%%%%%%%%% First we simulate risk distributions for mental health: this 
%%%%%%%%%%% is similar because risk draws here are assumed to be uncorrelated
%%%%%%%%%%% with those from the other three categories based on data analysis
%%%%%%%%%%% discuss in the Cost Model Estimation Appendix

%%%%%%%%%%% We first generate TOTAL EXPENDITURES. We determine OUT OF POCKET EXPENDITURES
%%%%%%%%%%% Based on plan specific designs later in this file 

%%%%%%%%%%% 50 Simulations of health risk draws to describe each individual distribution

Sim = 50;

%%%%%%%%%% Number of individuals in the data

nIs = size(famid,1);

%%%%%%%%%% Matrices to contain Mental Health Total Expenditure Draws over 3 years

MentalHealth1= zeros(nIs,Sim);
MentalHealth2 = zeros(nIs,Sim);
MentalHealth3 = zeros(nIs,Sim);

%%%%%%%%%% Loop generates 50 mental health draws for each individual based on underlying ex ante risk distribution
%%%%%%%%%% Distribution simulated from as the following properties, as estimated in the 'Cost Model Estimation' File
%%%%%%%%%% in STATA:

%%%%%%%%%% Two part distribution: Weibull with shape b and scale c, and second part with discrete mass at 0

%%%%%%%%%% PrZMH = Probability someone has 0 claims in category
%%%%%%%%%% MHb = Intercept for Weibull Shape parameter
%%%%%%%%%% MHGenderCoeffb = Shape parameter difference by gender
%%%%%%%%%% MHAgeCoeffb = Shape parameter slope as function of age
%%%%%%%%%% MHc = Intercept for Weibull Scale parameter
%%%%%%%%%% MHGenderCoeffc = Scale parameter difference by gender
%%%%%%%%%% MHAgeCoeffc = Scale parameter slope as function of age
%%%%%%%%%% Gender = 1 if male

%%%%% NOTE: All marginal distribution parameters are generated in cost estimation file, and based on ex ante cells with past mental health expenditures forming cells

for j = 1:Sim
    X = unifrnd(0,1);
    for i = 1:nIs

%%%%% If uniform draw > Probability of Zero, use Weibull distribution to find draw. Otherwise, draw = 0.
	
     if X >= PrZMH1(i,1)
     MentalHealth1(i,j) = wblrnd(MHb1(i,1)+MHGenderCoeffb1(i,1)*gender(i,1)+MHAgeCoeffb1(i,1)*age(i,1),MHc1(i,1)+MHGenderCoeffc1(i,1)*gender(i,1)+MHAgeCoeffc1(i,1)*age(i,1));
     else
     MentalHealth1(i,j)=0;
     end
     
     if X>= PrZMH2(i,1)
     MentalHealth2(i,j) = wblrnd(MHb2(i,1)+MHGenderCoeffb2(i,1)*gender(i,1)+MHAgeCoeffb2(i,1)*age(i,1),MHc2(i,1)+MHGenderCoeffc2(i,1)*gender(i,1)+MHAgeCoeffc2(i,1)*age(i,1));
     else
     MentalHealth2(i,j)=0;         
     end
     
     if X>= PrZMH3(i,1)
     MentalHealth3(i,j) = wblrnd(MHb3(i,1)+MHGenderCoeffb3(i,1)*gender(i,1)+MHAgeCoeffb3(i,1)*age(i,1),MHc3(i,1)+MHGenderCoeffc3(i,1)*gender(i,1)+MHAgeCoeffc3(i,1)*age(i,1));
     else
     MentalHealth3(i,j)=0;
     end
     
    end
end

%%%%%%%%%%%%% Now simulate draws for other three categories of claims / expenditures 
%%%%%%%%%%%%% More complicated: have to used copula methods to build in correlations 
%%%%%%%%%%%%% between health risk draws across categories. Described in Online Appendix
%%%%%%%%%%%%% on cost estimation. Reason to use copula is to use empirical correlations 
%%%%%%%%%%%%% to join marginal Weibull distribution estimates for expenditures in these 
%%%%%%%%%%%%% categories. 

%%%%%%%%%%%%% Crucial for copula methods: empirical correlations between categories conditional 
%%%%%%%%%%%%% on underlying risk cell are Spearman Rank Correlations. 

%%%% Matrices to hold draws

Pharmacy1 = zeros(nIs,Sim);
Pharmacy2 = zeros(nIs,Sim);
Pharmacy3 = zeros(nIs,Sim);

Hosp1 = zeros(nIs,Sim);
Hosp2 = zeros(nIs,Sim);
Hosp3 = zeros(nIs,Sim);

OfficeVis1 = zeros(nIs,Sim);
OfficeVis2 = zeros(nIs,Sim);
OfficeVis3 = zeros(nIs,Sim);

%%%%% Step one: use normal distribution copula: this represents correlations between 
%%%%% marginal Weibull distributions as correlation structure similar to that 
%%%%% in standard multivariate normal distribution. 

%%%%% Standard tri-variate normal mean and covariance matrices

MU = [0 0 0];
SIGMA1 = [1 0 0; 0 1 0; 0 0 1];
SIGMA2 = [1 0 0; 0 1 0; 0 0 1];
SIGMA3 = [1 0 0; 0 1 0; 0 0 1];

for j = 1:Sim

%%% X will be used to find 0 draws separately from positive draws described by joint Weibull distribution 
%%% (joint Weibull is marginal Weibull combined with copula methods)

    X = unifrnd(0,1);

for i = 1:nIs
   
%%%%%%%% Fill in standard normal variance covariance matrix 
%%%%%%%% with covariances derived from formula that links
%%%%%%%% Spearman rank correlations to covariances 

%%%%%%%% CorrRXH is rank correlation between RX and Hospital / other expenditures (for each health risk bin based on central ACG measure)
%%%%%%%% CorrOVH is rank correlation between Office Visit and Hospital / other expenditures (for each health risk bin based on central ACG measure)
%%%%%%%% CorrRXH is rank correlation between RX and Office Visit / other expenditures (for each health risk bin based on central ACG measure)
 
   SIGMA1(2,1) = sin((CorrRXH1(i,1)/6)*pi)*2;
   SIGMA1(1,2) = sin((CorrRXH1(i,1)/6)*pi)*2;
   SIGMA1(1,3) = sin((CorrOVH1(i,1)/6)*pi)*2; 
   SIGMA1(3,1) = sin((CorrOVH1(i,1)/6)*pi)*2;
   SIGMA1(3,2) = sin((CorrOVRX1(i,1)/6)*pi)*2;
   SIGMA1(2,3) = sin((CorrOVRX1(i,1)/6)*pi)*2;

   SIGMA2(2,1) = sin((CorrRXH2(i,1)/6)*pi)*2;
   SIGMA2(1,2) = sin((CorrRXH2(i,1)/6)*pi)*2;
   SIGMA2(1,3) = sin((CorrOVH2(i,1)/6)*pi)*2; 
   SIGMA2(3,1) = sin((CorrOVH2(i,1)/6)*pi)*2;
   SIGMA2(3,2) = sin((CorrOVRX2(i,1)/6)*pi)*2;
   SIGMA2(2,3) = sin((CorrOVRX2(i,1)/6)*pi)*2;

   SIGMA3(2,1) = sin((CorrRXH3(i,1)/6)*pi)*2;
   SIGMA3(1,2) = sin((CorrRXH3(i,1)/6)*pi)*2;
   SIGMA3(1,3) = sin((CorrOVH3(i,1)/6)*pi)*2; 
   SIGMA3(3,1) = sin((CorrOVH3(i,1)/6)*pi)*2;
   SIGMA3(3,2) = sin((CorrOVRX3(i,1)/6)*pi)*2;
   SIGMA3(2,3) = sin((CorrOVRX3(i,1)/6)*pi)*2;
   
%%%%%%%%% Now, with full standard normal copula, simulate tri-variate draw and convert each 
%%%%%%%%% draw for each category into a marginal cdf value for that category
%%%%%%%%% If correlations are high, then, e.g., trivariate draw will have high cdf values for 
%%%%%%%%% all categories, which we will then take to marginal Weibull distributions   
   
%% R is trivariate standard normal copula draw
%% HP is hospital / all other CDF value 
%% RXP is RX CDF value
%% OVP is office visit CDF value
   
   R1 = mvnrnd(MU, SIGMA1, 1);
   HP1 = normcdf(R1(1,1),0,SIGMA1(1,1));
   RXP1 = normcdf(R1(1,2),0,SIGMA1(2,2));
   OVP1 = normcdf(R1(1,3),0,SIGMA1(3,3));
   
   R2 = mvnrnd(MU, SIGMA2, 1);
   HP2 = normcdf(R2(1,1),0,SIGMA2(1,1));
   RXP2 = normcdf(R2(1,2),0,SIGMA2(2,2));
   OVP2 = normcdf(R2(1,3),0,SIGMA2(3,3));

   R3 = mvnrnd(MU, SIGMA3, 1);
   HP3 = normcdf(R3(1,1),0,SIGMA3(1,1));
   RXP3 = normcdf(R3(1,2),0,SIGMA3(2,2));
   OVP3 = normcdf(R3(1,3),0,SIGMA3(3,3));
   

%%%%%%%%%  Now for each category generate 50 draws given CDF draw from copula that takes into account 
%%%%%%%%%  correlations between categories. This, as in Mental health above, follows two part distribution.
%%%%%%%%%  First part determines if person has zero expendituers with discrete mass, second part simulates
%%%%%%%%%  from Weibull using CDF draw from normal copula above to take into account correlations. 

%% For RX:
%%%%%%%%%% PrZRX = Probability someone has 0 claims in category
%%%%%%%%%% RXb = Intercept for Weibull Shape parameter
%%%%%%%%%% RXGenderCoeffb = Shape parameter difference by gender
%%%%%%%%%% RXAgeCoeffb = Shape parameter slope as function of age
%%%%%%%%%% RXc = Intercept for Weibull Scale parameter
%%%%%%%%%% RXGenderCoeffc = Scale parameter difference by gender
%%%%%%%%%% RXAgeCoeffc = Scale parameter slope as function of age

%%%%% NOTE: All marginal distribution parameters are generated in cost estimation file, and based on ex ante cells with ACG-RX predictive measure 
%%%%% as cell basis for RX, and overall ACG expenditure predictive measure as basis for cells for hospital / other and office visit
   
   if X >= PrZRX1(i,1)
   Pharmacy1(i,j) = wblinv(RXP1,RXb1(i,1) + RXGenderCoeffb1(i,1)*gender(i,1)+RXAgeCoeffb1(i,1)*age(i,1), RXc1(i,1) + RXGenderCoeffc1(i,1)*gender(i,1)+RXAgeCoeffc1(i,1)*age(i,1));
   else 
   Pharmacy1(i,j)=0;    
   end

   if X >= PrZRX2(i,1)
   Pharmacy2(i,j) = wblinv(RXP2,RXb2(i,1) + RXGenderCoeffb2(i,1)*gender(i,1)+RXAgeCoeffb2(i,1)*age(i,1), RXc2(i,1) + RXGenderCoeffc2(i,1)*gender(i,1)+RXAgeCoeffc2(i,1)*age(i,1));
   else 
   Pharmacy2(i,j)=0;    
   end
   
   if X >= PrZRX3(i,1)
   Pharmacy3(i,j) = wblinv(RXP3,RXb3(i,1) + RXGenderCoeffb3(i,1)*gender(i,1)+RXAgeCoeffb3(i,1)*age(i,1), RXc3(i,1) + RXGenderCoeffc3(i,1)*gender(i,1)+RXAgeCoeffc3(i,1)*age(i,1));
   else 
   Pharmacy3(i,j)=0 ;   
   end
   
%% For OV:
%%%%%%%%%% PrZOV = Probability someone has 0 claims in category
%%%%%%%%%% OVb = Intercept for Weibull Shape parameter
%%%%%%%%%% OVc = Intercept for Weibull Scale parameter
     
   if X >= PrZOV1(i,1)
   OfficeVis1(i,j) = wblinv(OVP1,OVb1(i,1), OVc1(i,1));
   else
   OfficeVis1(i,j)=0;
   end

   if X >= PrZOV2(i,1)
   OfficeVis2(i,j) = wblinv(OVP2,OVb2(i,1), OVc2(i,1));
   else
   OfficeVis2(i,j)=0;
   end
   
   if X >= PrZOV3(i,1)
   OfficeVis3(i,j) = wblinv(OVP3,OVb3(i,1), OVc3(i,1));
   else
   OfficeVis3(i,j)=0;
   end 
   
%% For Hospital / All Other:
%%%%%%%%%% PrZHP = Probability someone has 0 claims in category
%%%%%%%%%% HPb = Intercept for Weibull Shape parameter
%%%%%%%%%% HPc = Intercept for Weibull Scale parameter 
 
   if X >= PrZHP1(i,1)
   Hosp1(i,j) = wblinv(HP1,Hb1(i,1), Hc1(i,1));
   else
   Hosp1(i,j)=0;
   end
   
   if X >= PrZHP2(i,1)
   Hosp2(i,j) = wblinv(HP2,Hb2(i,1), Hc2(i,1));
   else
   Hosp2(i,j)=0;
   end
   
   if X >= PrZHP3(i,1)
   Hosp3(i,j) = wblinv(HP3,Hb3(i,1), Hc3(i,1));
   else
   Hosp3(i,j)=0;
   end
  
%%%% Implement discrete probability estimated in cost model for very high hospital / other expenditures
%%%% This is the only category for which this high end discrete mass is relevant. Purpose of including 
%%%% is so Weibull estimates don't take into account very high draws which can substantially skew fit of distribution
%%%% to most of cost data realizations. Process, as described in Appendix, is to estimate discrete mass of individuals 
%%%% with expenditures higher than 40,000, equal to PrHHP. Then, replace with average expenditures for that group of individuals
%%%% in expenditure draws, equal to ~90,000.  
  
   if X >= (1 - PrHHP1(i,1))
   Hosp1(i,j) = 90000;
   end
   
   if X >= (1 - PrHHP2(i,1))
   Hosp2(i,j) = 90000;
   end
   
   if X >= (1 - PrHHP3(i,1))
   Hosp3(i,j) = 90000;
   end
      
   end
end

%%%%%%%%%%%%%% Now, for use in certain parts of analysis, construct sum of medical expenditures 
%%%%%%%%%%%%%% across categories for each simulated draw at the individual level. 

TotalY1 = zeros(nIs,Sim);
TotalY2 = zeros(nIs,Sim);
TotalY3 = zeros(nIs,Sim);

for i = 1:nIs
    for j = 1:Sim
        TotalY1(i,j) = Pharmacy1(i,j) + MentalHealth1(i,j) + OfficeVis1(i,j) + Hosp1(i,j);
        TotalY2(i,j) = Pharmacy2(i,j) + MentalHealth2(i,j) + OfficeVis2(i,j) + Hosp2(i,j);
        TotalY3(i,j) = Pharmacy3(i,j) + MentalHealth3(i,j) + OfficeVis3(i,j) + Hosp3(i,j);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% End of Simulation of Health Risk Distributions									%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Next portion of code maps total expenditures distributions from individuals     %%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Into plan-family-specific distribution of out-of-pocket expenditures            %%%%%%%%%%%%%%
%%%%%%%%%%%%%%% to be used as input into insurance plan choice model                            %%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% Map expenses to out of pocket expenditures using the mapping for each of three plans 
%%%%%%%%%%%%%%% Mapping verified elsewhere in STATA code, figures for mapping verification in Online Appendix. 

%%%%%%%%%%%%%%% Here, OOP mapping presented is for three actual plans studied in paper, appropriately renamed
%%%%%%%%%%%%%%% to protect identity of firm 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% 3 Plan Designs Described in Detail of Table 2 in main text of paper %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Not repeated here to keep code concise, refer to table for details  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matrices to Hold OOP draws for each Individual
% converted into family draws at end of this file
% Note: 1,2,3 in matrix name correspond to three years
% Third index corresponds to plan #

OutofPocket1 = zeros(nIs,Sim,3);
OutofPocket2 = zeros(nIs,Sim,3);
OutofPocket3 = zeros(nIs,Sim,3);

%%%%%%%%%%%%%%%% For first two plans, pharmacy claims and office visit claims and treated with copayments instead
%%%%%%%%%%%%%%%% of with % coinsurance or deductible. Copayments do not apply to the out of pocket maximums

%%%%%%%%%%%%%%%% To model amount paid in copayments, for each person, plan, and year for these two categories
%%%%%%%%%%%%%%%% in simplified manner:
%%%%%%%%%%%%%%%% (i) determine empirically mean % paid conditional on having given level of total category expenditures
%%%%%%%%%%%%%%%% (ii) treat this percentage like implicit coinsurance rate for each person 
%%%%%%%%%%%%%%%% (iii) Compute OOP given total expenditure draw in each of pharmacy and OV category for first two plans

%%%%%%%%%%%%%%%% NOTE: All copayment aspects of plans are subsumed in this reduced form calculation, since most expenses 
%%%%%%%%%%%%%%%% don't apply to this and apply to deductible, coinsurance, OOP max as follows later in code. 

RXoopCP1 = zeros(nIs,Sim);
OVoopCP1 = zeros(nIs,Sim);
RXoopCP2 = zeros(nIs,Sim);
OVoopCP2 = zeros(nIs,Sim);
RXoopCP3 = zeros(nIs,Sim);
OVoopCP3 = zeros(nIs,Sim);

% NOTE: copayment structures are identical for pharmacy and OV for first two plans. Given that, simple model of copayments
% mentioned above and implemented here is applicable to both plans.

for i =1:nIs
    for j=1:Sim
      
        if Pharmacy1(i,j) >=0 && Pharmacy1(i,j)<=200
            RXoopCP1(i,j) = Pharmacy1(i,j)/1.9; 
        end
        if Pharmacy1(i,j) >=200 && Pharmacy1(i,j)<=500
            RXoopCP1(i,j) = Pharmacy1(i,j)/2.2; 
        end        
        if Pharmacy1(i,j) >=500 && Pharmacy1(i,j)<=1000
            RXoopCP1(i,j) = Pharmacy1(i,j)/2.5; 
        end
        if Pharmacy1(i,j) >=1000 && Pharmacy1(i,j)<=2000
            RXoopCP1(i,j) = Pharmacy1(i,j)/3.4; 
        end
        if Pharmacy1(i,j) >=2000 && Pharmacy1(i,j)<=4000
            RXoopCP1(i,j) = Pharmacy1(i,j)/3.9; 
        end 
        if Pharmacy1(i,j) >=4000 && Pharmacy1(i,j)<=8000
            RXoopCP1(i,j) = Pharmacy1(i,j)/5.1; 
        end
        if Pharmacy1(i,j) >=8000
            RXoopCP1(i,j) = 1500 ;
        end

        if Pharmacy2(i,j) >=0 && Pharmacy2(i,j)<=200
            RXoopCP2(i,j) = Pharmacy2(i,j)/1.9; 
        end
        if Pharmacy2(i,j) >=200 && Pharmacy2(i,j)<=500
            RXoopCP2(i,j) = Pharmacy2(i,j)/2.2; 
        end        
        if Pharmacy2(i,j) >=500 && Pharmacy2(i,j)<=1000
            RXoopCP2(i,j) = Pharmacy2(i,j)/2.5; 
        end
        if Pharmacy2(i,j) >=1000 && Pharmacy2(i,j)<=2000
            RXoopCP2(i,j) = Pharmacy2(i,j)/3.4; 
        end
        if Pharmacy2(i,j) >=2000 && Pharmacy2(i,j)<=4000
            RXoopCP2(i,j) = Pharmacy2(i,j)/3.9; 
        end 
        if Pharmacy2(i,j) >=4000 && Pharmacy2(i,j)<=8000
            RXoopCP2(i,j) = Pharmacy2(i,j)/5.1; 
        end
        if Pharmacy2(i,j) >=8000
            RXoopCP2(i,j) = 1500; 
        end
        
		%%%%% NOTE: Year 3 ratios increased below since copayment rates went up for that year by a small amount. 
		
        if Pharmacy3(i,j) >=0 && Pharmacy3(i,j)<=200
            RXoopCP3(i,j) = 1.4*Pharmacy3(i,j)/1.9; 
        end
        if Pharmacy3(i,j) >=200 && Pharmacy3(i,j)<=500
            RXoopCP3(i,j) = 1.4*Pharmacy3(i,j)/2.2; 
        end        
        if Pharmacy3(i,j) >=500 && Pharmacy3(i,j)<=1000
            RXoopCP3(i,j) = 1.4*Pharmacy3(i,j)/2.5; 
        end
        if Pharmacy3(i,j) >=1000 && Pharmacy3(i,j)<=2000
            RXoopCP3(i,j) = 1.4*Pharmacy3(i,j)/3.4; 
        end
        if Pharmacy3(i,j) >=2000 && Pharmacy3(i,j)<=4000
            RXoopCP3(i,j) = 1.4*Pharmacy3(i,j)/3.9; 
        end 
        if Pharmacy3(i,j) >=4000 && Pharmacy3(i,j)<=6500
            RXoopCP3(i,j) = 1.4*Pharmacy3(i,j)/5.1; 
        end
        if Pharmacy3(i,j) >=6500
            RXoopCP3(i,j) = 1500; 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Now Office Visit
        
        if OfficeVis1(i,j)>=0 && OfficeVis1(i,j)<=200
           OVoopCP1(i,j) = OfficeVis1(i,j)/3.3; 
        end
        if OfficeVis1(i,j)>=200 && OfficeVis1(i,j)<=500
           OVoopCP1(i,j) = OfficeVis1(i,j)/3.55; 
        end
        if OfficeVis1(i,j)>=500 && OfficeVis1(i,j)<=1200
           OVoopCP1(i,j) = OfficeVis1(i,j)/4.0 ;
        end
        if OfficeVis1(i,j)>=1500 && OfficeVis1(i,j)<=3000
           OVoopCP1(i,j) = OfficeVis1(i,j)/4.7; 
        end
        if OfficeVis1(i,j)>=3000
           OVoopCP1(i,j) = OfficeVis1(i,j)/5.4 ;
        end

        if OfficeVis2(i,j)>=0 && OfficeVis2(i,j)<=200
           OVoopCP2(i,j) = OfficeVis2(i,j)/3.3 ;
        end
        if OfficeVis2(i,j)>=200 && OfficeVis2(i,j)<=500
           OVoopCP2(i,j) = OfficeVis2(i,j)/3.55 ;
        end
        if OfficeVis2(i,j)>=500 && OfficeVis2(i,j)<=1200
           OVoopCP2(i,j) = OfficeVis2(i,j)/4.0 ;
        end
        if OfficeVis2(i,j)>=1500 && OfficeVis2(i,j)<=3000
           OVoopCP2(i,j) = OfficeVis2(i,j)/4.7 ;
        end
        if OfficeVis2(i,j)>=3000
           OVoopCP2(i,j) = OfficeVis2(i,j)/5.4 ;
        end
        
%%%%% NOTE: Year 3 ratios increased below since copayment rates went up for that year by a small amount. 
		
        if OfficeVis3(i,j)>=0 && OfficeVis3(i,j)<=200
           OVoopCP3(i,j) = 1.05*OfficeVis3(i,j)/3.3; 
        end
        if OfficeVis3(i,j)>=200 && OfficeVis3(i,j)<=500
           OVoopCP3(i,j) = 1.05*OfficeVis3(i,j)/3.55; 
        end
        if OfficeVis3(i,j)>=500 && OfficeVis3(i,j)<=1200
           OVoopCP3(i,j) = 1.05*OfficeVis3(i,j)/4.0 ;
        end
        if OfficeVis3(i,j)>=1500 && OfficeVis3(i,j)<=3000
           OVoopCP3(i,j) = 1.05*OfficeVis3(i,j)/4.7 ;
        end
        if OfficeVis3(i,j)>=3000
           OVoopCP3(i,j) = 1.05*OfficeVis3(i,j)/5.4 ;
        end
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Now, begin portion of code where, for each of three plans, we 
%%%%%%% map total medical expenditures across four categories into one 
%%%%%%% out-of-pocket expenditure amount per draw (of 50 draws). Note that 
%%%%%%% first mapping is done for each individual within family. Then, if family 
%%%%%%% deductible / OOP max constraints bind we correct for this at end of code
%%%%%%% by accordingly adjusting family OOP draws 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% Generate OOP draws for Plan 1 - PPO250 - Year 1

%%%%%%%%%%%%%%% (i) function of income (grouped into 5 tiers) because OOP-max depends on income 
%%%%%%%%%%%%%%% NOTE: I comment calculation details for income tier 1 and plan 1 in year 1. 
%%%%%%%%%%%%%%% Analagous details apply to other plans, income tiers, and year. See table 2 
%%%%%%%%%%%%%%% in main text for additional details on plan designs.

for i = 1:nIs
    for j = 1:Sim

%%%%%%%%%%%%%%% OOP expenditures given total costs, individuals in faimlies with lowest income tier 
	
        if income1(i,1) == 1 

%%%%%%%%%%%%%%% OOP from Hosp / All Other category when no Mental Health Expenditures (treated at different coinsurance but still apply to deductible) 
%%%%%%%%%%%%%%% Scenario 1: Spend Less than OOP-Max but more than deductible

		if Hosp1(i,j) < 7750 && Hosp1(i,j) >=250 && MentalHealth1(i,j) ==0
               OutofPocket1(i,j,1) = 250 + (Hosp1(i,j)-250)*0.1;
        end

%%%%%%%%%%%%%%% Scenario 2: Spend enough that you hit OOP max

		if Hosp1(i,j) >= 7750 && MentalHealth1(i,j) ==0
               OutofPocket1(i,j,1) = 1000;
        end
        
%%%%%%%%%%%%%%% Positive Mental Health Expenditures. Proportional rationing assumed for application to deductible if combined
%%%%%%%%%%%%%%% expenses larger than deductible. MH applies only to deductible, 50% MH coinsurance post-ded does not apply to OOP max. 

		if MentalHealth1(i,j) > 0 && Hosp1(i,j) < 7750 && (Hosp1(i,j) + MentalHealth1(i,j)) >=250
               OutofPocket1(i,j,1) = 250 + (Hosp1(i,j)-250*(1-MentalHealth1(i,j)/(MentalHealth1(i,j) +Hosp1(i,j) +0.001)))*0.1 + (MentalHealth1(i,j)-250*MentalHealth1(i,j)/(MentalHealth1(i,j) +Hosp1(i,j) +0.001))*0.5;  
        end

%%%%%%%%%%%%%%% Case where Hospital / Other expenditures hit OOP max and residual coinsurance from mental health
		
		if MentalHealth1(i,j) > 0 && Hosp1(i,j) >= 7750 
               OutofPocket1(i,j,1) = 1000 + (MentalHealth1(i,j)-250*MentalHealth1(i,j)/(MentalHealth1(i,j) +Hosp1(i,j) +0.001))*0.5;  
        end  
    end

%%%%%%%%%%%%%%%%% Plan 1 - PPO250 - Year 1 Income Tiers 2 and 3	
       
        if (income1(i,1) == 2 | income1(i,1)==3) 
            if Hosp1(i,j) < 17750 && Hosp1(i,j) >=250 && MentalHealth1(i,j) ==0
               OutofPocket1(i,j,1) = 250 + (Hosp1(i,j)-250)*0.1;
            end
            if Hosp1(i,j) >= 17750 && MentalHealth1(i,j) ==0
               OutofPocket1(i,j,1) = 2000;
            end
            if MentalHealth1(i,j) > 0 && Hosp1(i,j) < 17750 && (Hosp1(i,j) + MentalHealth1(i,j)) >=250
               OutofPocket1(i,j,1) = 250 + (Hosp1(i,j)-250*(1-MentalHealth1(i,j)/(MentalHealth1(i,j) +Hosp1(i,j) +0.001)))*0.1 + (MentalHealth1(i,j)-250*MentalHealth1(i,j)/(MentalHealth1(i,j) +Hosp1(i,j) +0.001))*0.5;  
            end
            if MentalHealth1(i,j) > 0 && Hosp1(i,j) >= 17750 
               OutofPocket1(i,j,1) = 2000 + (MentalHealth1(i,j)-250*MentalHealth1(i,j)/(MentalHealth1(i,j) +Hosp1(i,j) +0.001))*0.5;  
            end  
        end        

%%%%%%%%%%%%%%%%% Plan 1 - PPO250 - Year 1 Income Tiers 4 and 5	
		
        if (income1(i,1) == 4 | income1(i,1)==5) 
            if Hosp1(i,j) < 27750 && Hosp1(i,j) >=250 && MentalHealth1(i,j) ==0
               OutofPocket1(i,j,1) = 250 + (Hosp1(i,j)-250)*0.1;
            end
            if Hosp1(i,j) >= 27750 && MentalHealth1(i,j) ==0
               OutofPocket1(i,j,1) = 3000;
            end
            if MentalHealth1(i,j) > 0 && Hosp1(i,j) < 27750 && (Hosp1(i,j) + MentalHealth1(i,j)) >=250
               OutofPocket1(i,j,1) = 250 + (Hosp1(i,j)-250*(1-MentalHealth1(i,j)/(MentalHealth1(i,j) +Hosp1(i,j) +0.001)))*0.1 + (MentalHealth1(i,j)-250*MentalHealth1(i,j)/(MentalHealth1(i,j) +Hosp1(i,j) +0.001))*0.5;  
            end
            if MentalHealth1(i,j) > 0 && Hosp1(i,j) >= 27750 
               OutofPocket1(i,j,1) = 3000 + (MentalHealth1(i,j)-250*MentalHealth1(i,j)/(MentalHealth1(i,j) +Hosp1(i,j) +0.001))*0.5;  
            end  
        end        
    
%%%%%%%%%%%%%%%%% OOP expenditures equal total expenditures for all income tiers if Mental Health and Hospital / All Other < deductible

        if Hosp1(i,j) < 250 && MentalHealth1(i,j)==0
               OutofPocket1(i,j,1) = Hosp1(i,j);
        end 
        
        if (Hosp1(i,j) + MentalHealth1(i,j) <=250) 
               OutofPocket1(i,j,1) = Hosp1(i,j) + MentalHealth1(i,j);
        end

%%%%%%%%%%%%%%%% Finally: this plan has copayments for pharmacy / OOP visits, where implicit coinsurance rate calculated in 
%%%%%%%%%%%%%%%% STATA cost model estimation code and implemented above in this file. Applies to both PPO250 and PPO500
%%%%%%%%%%%%%%%% Since both have same copayment rates for these services 
%%%%%%%%%%%%%%%% OOP projected from these two categories added on for each draw right here, to form total plan OOP for each individual draw and year
		
         OutofPocket1(i,j,1)=OutofPocket1(i,j,1) + OVoopCP1(i,j) + RXoopCP1(i,j);     
         
    end
end

%%%%%%%%%%%%%%% Generate OOP draws for Plan 1 - PPO250 - Year 2

for i = 1:nIs
    for j = 1:Sim
        
        if income2(i,1) == 1 
            if Hosp2(i,j) < 7750 && Hosp2(i,j) >=250 && MentalHealth2(i,j) ==0
               OutofPocket2(i,j,1) = 250 + (Hosp2(i,j)-250)*0.1;
            end
            if Hosp2(i,j) >= 7750 && MentalHealth2(i,j) ==0
               OutofPocket2(i,j,1) = 1000;
            end
            if MentalHealth2(i,j) > 0 && Hosp2(i,j) < 7750 && (Hosp2(i,j) + MentalHealth2(i,j)) >=250
               OutofPocket2(i,j,1) = 250 + (Hosp2(i,j)-250*(1-MentalHealth2(i,j)/(MentalHealth2(i,j) +Hosp2(i,j) +0.001)))*0.1 + (MentalHealth2(i,j)-250*MentalHealth2(i,j)/(MentalHealth2(i,j) +Hosp2(i,j) +0.001))*0.5;  
            end
            if MentalHealth2(i,j) > 0 && Hosp2(i,j) >= 7750 
               OutofPocket2(i,j,1) = 1000 + (MentalHealth2(i,j)-250*MentalHealth2(i,j)/(MentalHealth2(i,j) +Hosp2(i,j) +0.001))*0.5;  
            end  
        end
        
        if (income2(i,1) == 2 | income2(i,1)==3) 
            if Hosp2(i,j) < 17750 && Hosp2(i,j) >=250 && MentalHealth2(i,j) ==0
               OutofPocket2(i,j,1) = 250 + (Hosp2(i,j)-250)*0.1;
            end
            if Hosp2(i,j) >= 17750 && MentalHealth2(i,j) ==0
               OutofPocket2(i,j,1) = 2000;
            end
            if MentalHealth2(i,j) > 0 && Hosp2(i,j) < 17750 && (Hosp2(i,j) + MentalHealth2(i,j)) >=250
               OutofPocket2(i,j,1) = 250 + (Hosp2(i,j)-250*(1-MentalHealth2(i,j)/(MentalHealth2(i,j) +Hosp2(i,j) +0.001)))*0.1 + (MentalHealth2(i,j)-250*MentalHealth2(i,j)/(MentalHealth2(i,j) +Hosp2(i,j) +0.001))*0.5;  
            end
            if MentalHealth2(i,j) > 0 && Hosp2(i,j) >= 17750 
               OutofPocket2(i,j,1) = 2000 + (MentalHealth2(i,j)-250*MentalHealth2(i,j)/(MentalHealth2(i,j) +Hosp2(i,j) +0.001))*0.5;  
            end  
        end        

        if (income2(i,1) == 4 | income2(i,1)==5) 
            if Hosp2(i,j) < 27750 && Hosp2(i,j) >=250 && MentalHealth2(i,j) ==0
               OutofPocket2(i,j,1) = 250 + (Hosp2(i,j)-250)*0.1;
            end
            if Hosp2(i,j) >= 27750 && MentalHealth2(i,j) ==0
               OutofPocket2(i,j,1) = 3000;
            end
            if MentalHealth2(i,j) > 0 && Hosp2(i,j) < 27750 && (Hosp2(i,j) + MentalHealth2(i,j)) >=250
               OutofPocket2(i,j,1) = 250 + (Hosp2(i,j)-250*(1-MentalHealth2(i,j)/(MentalHealth2(i,j) +Hosp2(i,j) +0.001)))*0.1 + (MentalHealth2(i,j)-250*MentalHealth2(i,j)/(MentalHealth2(i,j) +Hosp2(i,j) +0.001))*0.5;  
            end
            if MentalHealth2(i,j) > 0 && Hosp2(i,j) >= 27750 
               OutofPocket2(i,j,1) = 3000 + (MentalHealth2(i,j)-250*MentalHealth2(i,j)/(MentalHealth2(i,j) +Hosp2(i,j) +0.001))*0.5;  
            end  
        end        
        
        if Hosp2(i,j) < 250 && MentalHealth2(i,j)==0
               OutofPocket2(i,j,1) = Hosp2(i,j);
        end 
        
        if (Hosp2(i,j) + MentalHealth2(i,j) <=250)
               OutofPocket2(i,j,1) = Hosp2(i,j) + MentalHealth2(i,j);
        end
    
     OutofPocket2(i,j,1)=OutofPocket2(i,j,1) + OVoopCP2(i,j) + RXoopCP2(i,j); 
        
    end
end

%%%%%%%%%%%%%%% Generate OOP draws for Plan 1 - PPO250 - Year 3

for i = 1:nIs
    for j = 1:Sim
        
        if income2(i,1) == 1 
            if Hosp3(i,j) < 7750 && Hosp3(i,j) >=250 && MentalHealth3(i,j) ==0
               OutofPocket3(i,j,1) = 250 + (Hosp3(i,j)-250)*0.1;
            end
            if Hosp3(i,j) >= 7750 && MentalHealth3(i,j) ==0
               OutofPocket3(i,j,1) = 1000;
            end
            if MentalHealth3(i,j) > 0 && Hosp3(i,j) < 7750 && (Hosp3(i,j) + MentalHealth3(i,j)) >=250
               OutofPocket3(i,j,1) = 250 + (Hosp3(i,j)-250*(1-MentalHealth3(i,j)/(MentalHealth3(i,j) +Hosp3(i,j) +0.001)))*0.1 + (MentalHealth3(i,j)-250*MentalHealth3(i,j)/(MentalHealth3(i,j) +Hosp3(i,j) +0.001))*0.5;  
            end
            if MentalHealth3(i,j) > 0 && Hosp3(i,j) >= 7750 
               OutofPocket3(i,j,1) = 1000 + (MentalHealth3(i,j)-250*MentalHealth3(i,j)/(MentalHealth3(i,j) +Hosp3(i,j) +0.001))*0.5;  
            end  
        end
        
        if (income2(i,1) == 2 | income2(i,1)==3) 
            if Hosp3(i,j) < 17750 && Hosp3(i,j) >=250 && MentalHealth3(i,j) ==0
               OutofPocket3(i,j,1) = 250 + (Hosp3(i,j)-250)*0.1;
            end
            if Hosp3(i,j) >= 17750 && MentalHealth3(i,j) ==0
               OutofPocket3(i,j,1) = 2000;
            end
            if MentalHealth3(i,j) > 0 && Hosp3(i,j) < 17750 && (Hosp3(i,j) + MentalHealth3(i,j)) >=250
               OutofPocket3(i,j,1) = 250 + (Hosp3(i,j)-250*(1-MentalHealth3(i,j)/(MentalHealth3(i,j) +Hosp3(i,j) +0.001)))*0.1 + (MentalHealth3(i,j)-250*MentalHealth3(i,j)/(MentalHealth3(i,j) +Hosp3(i,j) +0.001))*0.5;  
            end
            if MentalHealth3(i,j) > 0 && Hosp3(i,j) >= 17750 
               OutofPocket3(i,j,1) = 2000 + (MentalHealth3(i,j)-250*MentalHealth3(i,j)/(MentalHealth3(i,j) +Hosp3(i,j) +0.001))*0.5;  
            end  
        end        

        if (income2(i,1) == 4 | income2(i,1)==5) 
            if Hosp3(i,j) < 27750 && Hosp3(i,j) >=250 && MentalHealth3(i,j) ==0
               OutofPocket3(i,j,1) = 250 + (Hosp3(i,j)-250)*0.1;
            end
            if Hosp3(i,j) >= 27750 && MentalHealth3(i,j) ==0
               OutofPocket3(i,j,1) = 3000;
            end
            if MentalHealth3(i,j) > 0 && Hosp3(i,j) < 27750 && (Hosp3(i,j) + MentalHealth3(i,j)) >=250
               OutofPocket3(i,j,1) = 250 + (Hosp3(i,j)-250*(1-MentalHealth3(i,j)/(MentalHealth3(i,j) +Hosp3(i,j) +0.001)))*0.1 + (MentalHealth3(i,j)-250*MentalHealth3(i,j)/(MentalHealth3(i,j) +Hosp3(i,j) +0.001))*0.5;  
            end
            if MentalHealth3(i,j) > 0 && Hosp3(i,j) >= 27750 
               OutofPocket3(i,j,1) = 3000 + (MentalHealth3(i,j)-250*MentalHealth3(i,j)/(MentalHealth3(i,j) +Hosp3(i,j) +0.001))*0.5;  
            end  
        end        
        
        if Hosp3(i,j) < 250 && MentalHealth3(i,j)==0
               OutofPocket3(i,j,1) = Hosp3(i,j);
        end 
        
        if (Hosp3(i,j) + MentalHealth3(i,j) <=250)
               OutofPocket3(i,j,1) = Hosp3(i,j) + MentalHealth3(i,j);
        end
    
     OutofPocket3(i,j,1)=OutofPocket3(i,j,1) + OVoopCP3(i,j) + RXoopCP3(i,j);   
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% Generate OOP draws for Plan 2 - PPO500 - Year 1

for i = 1:nIs
    for j = 1:Sim
        
        if income1(i,1) == 1 
            if Hosp1(i,j) < 5500 && Hosp1(i,j) >=500 && MentalHealth1(i,j) ==0
               OutofPocket1(i,j,2) = 500 + (Hosp1(i,j)-500)*0.2;
            end
            if Hosp1(i,j) >= 5500 && MentalHealth1(i,j) ==0
               OutofPocket1(i,j,2) = 1500;
            end
            if MentalHealth1(i,j) > 0 && Hosp1(i,j) < 5500 && (Hosp1(i,j) + MentalHealth1(i,j)) >=500
               OutofPocket1(i,j,2) = 500 + (Hosp1(i,j)-500*(1-MentalHealth1(i,j)/(MentalHealth1(i,j) +Hosp1(i,j) +0.001)))*0.2 + (MentalHealth1(i,j)-500*MentalHealth1(i,j)/(MentalHealth1(i,j) +Hosp1(i,j) +0.001))*0.5;  
            end
            if MentalHealth1(i,j) > 0 && Hosp1(i,j) >= 5500 
               OutofPocket1(i,j,2) = 1500 + (MentalHealth1(i,j)-500*MentalHealth1(i,j)/(MentalHealth1(i,j) +Hosp1(i,j) +0.001))*0.5;  
            end  
        end
        
        if (income1(i,1) == 2 | income1(i,1)==3) 
            if Hosp1(i,j) < 13000 && Hosp1(i,j) >=500 && MentalHealth1(i,j) ==0
               OutofPocket1(i,j,2) = 500 + (Hosp1(i,j)-500)*0.2;
            end
            if Hosp1(i,j) >= 13000 && MentalHealth1(i,j) ==0
               OutofPocket1(i,j,2) = 3000;
            end
            if MentalHealth1(i,j) > 0 && Hosp1(i,j) < 13000 && (Hosp1(i,j) + MentalHealth1(i,j)) >=500
               OutofPocket1(i,j,2) = 500 + (Hosp1(i,j)-500*(1-MentalHealth1(i,j)/(MentalHealth1(i,j) +Hosp1(i,j) +0.001)))*0.2 + (MentalHealth1(i,j)-500*MentalHealth1(i,j)/(MentalHealth1(i,j) +Hosp1(i,j) +0.001))*0.5;  
            end
            if MentalHealth1(i,j) > 0 && Hosp1(i,j) >= 13000 
               OutofPocket1(i,j,2) = 3000 + (MentalHealth1(i,j)-500*MentalHealth1(i,j)/(MentalHealth1(i,j) +Hosp1(i,j) +0.001))*0.5;  
            end  
        end        

        if (income1(i,1) == 4 | income1(i,1)==5) 
            if Hosp1(i,j) < 18000 && Hosp1(i,j) >=500 && MentalHealth1(i,j) ==0
               OutofPocket1(i,j,2) = 500 + (Hosp1(i,j)-500)*0.2;
            end
            if Hosp1(i,j) >= 18000 && MentalHealth1(i,j) ==0
               OutofPocket1(i,j,2) = 4000;
            end
            if MentalHealth1(i,j) > 0 && Hosp1(i,j) < 18000 && (Hosp1(i,j) + MentalHealth1(i,j)) >=500
               OutofPocket1(i,j,2) = 500 + (Hosp1(i,j)-500*(1-MentalHealth1(i,j)/(MentalHealth1(i,j) +Hosp1(i,j) +0.001)))*0.2 + (MentalHealth1(i,j)-500*MentalHealth1(i,j)/(MentalHealth1(i,j) +Hosp1(i,j) +0.001))*0.5;  
            end
            if MentalHealth1(i,j) > 0 && Hosp1(i,j) >= 18000 
               OutofPocket1(i,j,2) = 4000 + (MentalHealth1(i,j)-500*MentalHealth1(i,j)/(MentalHealth1(i,j) +Hosp1(i,j) +0.001))*0.5;  
            end  
        end        
        
        if Hosp1(i,j) < 500 && MentalHealth1(i,j)==0
               OutofPocket1(i,j,2) = Hosp1(i,j);
        end 
        
        if (Hosp1(i,j) + MentalHealth1(i,j) <=500)
               OutofPocket1(i,j,2) = Hosp1(i,j) + MentalHealth1(i,j);
        end
    
     OutofPocket1(i,j,2)=OutofPocket1(i,j,2) + OVoopCP1(i,j) + RXoopCP1(i,j);      
    end
end

%%%%%%%%%%%%%%% Generate OOP draws for Plan 2 - PPO500 - Year 2

for i = 1:nIs
    for j = 1:Sim
        
        if income2(i,1) == 1 
            if Hosp2(i,j) < 5500 && Hosp2(i,j) >=500 && MentalHealth2(i,j) ==0
               OutofPocket2(i,j,2) = 500 + (Hosp2(i,j)-500)*0.2;
            end
            if Hosp2(i,j) >= 5500 && MentalHealth2(i,j) ==0
               OutofPocket2(i,j,2) = 1500;
            end
            if MentalHealth2(i,j) > 0 && Hosp2(i,j) < 5500 && (Hosp2(i,j) + MentalHealth2(i,j)) >=500
               OutofPocket2(i,j,2) = 500 + (Hosp2(i,j)-500*(1-MentalHealth2(i,j)/(MentalHealth2(i,j) +Hosp2(i,j) +0.001)))*0.2 + (MentalHealth2(i,j)-500*MentalHealth2(i,j)/(MentalHealth2(i,j) +Hosp2(i,j) +0.001))*0.5;  
            end
            if MentalHealth2(i,j) > 0 && Hosp2(i,j) >= 5500 
               OutofPocket2(i,j,2) = 1500 + (MentalHealth2(i,j)-500*MentalHealth2(i,j)/(MentalHealth2(i,j) +Hosp2(i,j) +0.001))*0.5;  
            end  
        end
        
        if (income2(i,1) == 2 | income2(i,1)==3) 
            if Hosp2(i,j) < 13000 && Hosp2(i,j) >=500 && MentalHealth2(i,j) ==0
               OutofPocket2(i,j,2) = 500 + (Hosp2(i,j)-500)*0.2;
            end
            if Hosp2(i,j) >= 13000 && MentalHealth2(i,j) ==0
               OutofPocket2(i,j,2) = 3000;
            end
            if MentalHealth2(i,j) > 0 && Hosp2(i,j) < 13000 && (Hosp2(i,j) + MentalHealth2(i,j)) >=500
               OutofPocket2(i,j,2) = 500 + (Hosp2(i,j)-500*(1-MentalHealth2(i,j)/(MentalHealth2(i,j) +Hosp2(i,j) +0.001)))*0.2 + (MentalHealth2(i,j)-500*MentalHealth2(i,j)/(MentalHealth2(i,j) +Hosp2(i,j) +0.001))*0.5;  
            end
            if MentalHealth2(i,j) > 0 && Hosp2(i,j) >= 13000 
               OutofPocket2(i,j,2) = 3000 + (MentalHealth2(i,j)-500*MentalHealth2(i,j)/(MentalHealth2(i,j) +Hosp2(i,j) +0.001))*0.5;  
            end  
        end        

        if (income2(i,1) == 4 | income2(i,1)==5) 
            if Hosp2(i,j) < 18000 && Hosp2(i,j) >=500 && MentalHealth2(i,j) ==0
               OutofPocket2(i,j,2) = 500 + (Hosp2(i,j)-500)*0.2;
            end
            if Hosp2(i,j) >= 18000 && MentalHealth2(i,j) ==0
               OutofPocket2(i,j,2) = 4000;
            end
            if MentalHealth2(i,j) > 0 && Hosp2(i,j) < 18000 && (Hosp2(i,j) + MentalHealth2(i,j)) >=500
               OutofPocket2(i,j,2) = 500 + (Hosp2(i,j)-500*(1-MentalHealth2(i,j)/(MentalHealth2(i,j) +Hosp2(i,j) +0.001)))*0.2 + (MentalHealth2(i,j)-500*MentalHealth2(i,j)/(MentalHealth2(i,j) +Hosp2(i,j) +0.001))*0.5;  
            end
            if MentalHealth2(i,j) > 0 && Hosp2(i,j) >= 18000 
               OutofPocket2(i,j,2) = 4000 + (MentalHealth2(i,j)-500*MentalHealth2(i,j)/(MentalHealth2(i,j) +Hosp2(i,j) +0.001))*0.5;  
            end  
        end        
        
        if Hosp2(i,j) < 500 && MentalHealth2(i,j)==0
               OutofPocket2(i,j,2) = Hosp2(i,j);
        end 
        
        if (Hosp2(i,j) + MentalHealth2(i,j) <=500)
               OutofPocket2(i,j,2) = Hosp2(i,j) + MentalHealth2(i,j);
        end
    
  OutofPocket2(i,j,2)=OutofPocket2(i,j,2) + OVoopCP2(i,j) + RXoopCP2(i,j);     
    end
end

%%%%%%%%%%%%%%% Generate OOP draws for Plan 2 - PPO500 - Year 3


for i = 1:nIs
    for j = 1:Sim
        
        if income2(i,1) == 1 
            if Hosp3(i,j) < 5500 && Hosp3(i,j) >=500 && MentalHealth3(i,j) ==0
               OutofPocket3(i,j,2) = 500 + (Hosp3(i,j)-500)*0.2;
            end
            if Hosp3(i,j) >= 5500 && MentalHealth3(i,j) ==0
               OutofPocket3(i,j,2) = 1500;
            end
            if MentalHealth3(i,j) > 0 && Hosp3(i,j) < 5500 && (Hosp3(i,j) + MentalHealth3(i,j)) >=500
               OutofPocket3(i,j,2) = 500 + (Hosp3(i,j)-500*(1-MentalHealth3(i,j)/(MentalHealth3(i,j) +Hosp3(i,j) +0.001)))*0.2 + (MentalHealth3(i,j)-500*MentalHealth3(i,j)/(MentalHealth3(i,j) +Hosp3(i,j) +0.001))*0.5;  
            end
            if MentalHealth3(i,j) > 0 && Hosp3(i,j) >= 5500 
               OutofPocket3(i,j,2) = 1500 + (MentalHealth3(i,j)-500*MentalHealth3(i,j)/(MentalHealth3(i,j) +Hosp3(i,j) +0.001))*0.5;  
            end  
        end
        
        if (income2(i,1) == 2 | income2(i,1)==3) 
            if Hosp3(i,j) < 13000 && Hosp3(i,j) >=500 && MentalHealth3(i,j) ==0
               OutofPocket3(i,j,2) = 500 + (Hosp3(i,j)-500)*0.2;
            end
            if Hosp3(i,j) >= 13000 && MentalHealth3(i,j) ==0
               OutofPocket3(i,j,2) = 3000;
            end
            if MentalHealth3(i,j) > 0 && Hosp3(i,j) < 13000 && (Hosp3(i,j) + MentalHealth3(i,j)) >=500
               OutofPocket3(i,j,2) = 500 + (Hosp3(i,j)-500*(1-MentalHealth3(i,j)/(MentalHealth3(i,j) +Hosp3(i,j) +0.001)))*0.2 + (MentalHealth3(i,j)-500*MentalHealth3(i,j)/(MentalHealth3(i,j) +Hosp3(i,j) +0.001))*0.5;  
            end
            if MentalHealth3(i,j) > 0 && Hosp3(i,j) >= 13000 
               OutofPocket3(i,j,2) = 3000 + (MentalHealth3(i,j)-500*MentalHealth3(i,j)/(MentalHealth3(i,j) +Hosp3(i,j) +0.001))*0.5;  
            end  
        end        

        if (income2(i,1) == 4 | income2(i,1)==5) 
            if Hosp3(i,j) < 18000 && Hosp3(i,j) >=500 && MentalHealth3(i,j) ==0
               OutofPocket3(i,j,2) = 500 + (Hosp3(i,j)-500)*0.2;
            end
            if Hosp3(i,j) >= 18000 && MentalHealth3(i,j) ==0
               OutofPocket3(i,j,2) = 4000;
            end
            if MentalHealth3(i,j) > 0 && Hosp3(i,j) < 18000 && (Hosp3(i,j) + MentalHealth3(i,j)) >=500
               OutofPocket3(i,j,2) = 500 + (Hosp3(i,j)-500*(1-MentalHealth3(i,j)/(MentalHealth3(i,j) +Hosp3(i,j) +0.001)))*0.2 + (MentalHealth3(i,j)-500*MentalHealth3(i,j)/(MentalHealth3(i,j) +Hosp3(i,j) +0.001))*0.5;  
            end
            if MentalHealth3(i,j) > 0 && Hosp3(i,j) >= 18000 
               OutofPocket3(i,j,2) = 4000 + (MentalHealth3(i,j)-500*MentalHealth3(i,j)/(MentalHealth3(i,j) +Hosp3(i,j) +0.001))*0.5;  
            end  
        end        
        
        if Hosp3(i,j) < 500 && MentalHealth3(i,j)==0
               OutofPocket3(i,j,2) = Hosp3(i,j);
        end 
        
        if (Hosp3(i,j) + MentalHealth3(i,j) <=500)
               OutofPocket3(i,j,2) = Hosp3(i,j) + MentalHealth3(i,j);
        end
    
   OutofPocket3(i,j,2)=OutofPocket3(i,j,2) + OVoopCP3(i,j) + RXoopCP3(i,j);    
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% Generate OOP draws for Plan 3 - PPO1200 

%%%%%%%%%%%%%%% NOTE: Now Pharmacy and Office Visit Apply to Deductible and Coinsurance
%%%%%%%%%%%%%%% Just like Hospital / All Other so Copayments don't matter for this plan.
%%%%%%%%%%%%%%% Mental health applies to deductible, has its own coinsurance, and 
%%%%%%%%%%%%%%% coinsurance does not apply to OOP max as in other two plans. 

%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Construct Measures for All Non-Mental Health Expenditures 
%%%%%%%%%%%%%%% per person, draw, and year, since all apply to plan OOP
%%%%%%%%%%%%%%% Mapping in the same way
%%%%%%%%%%%%%%%

HDHPC1 = zeros(nIs,Sim);
HDHPC1 = Hosp1 + OfficeVis1 + Pharmacy1; 

HDHPC2 = zeros(nIs,Sim);
HDHPC2 = Hosp2 + OfficeVis2 + Pharmacy2;

HDHPC3 = zeros(nIs,Sim);
HDHPC3 = Hosp3 + OfficeVis3 + Pharmacy3;

%%%%%%%%%% Adjustment to total because preventive expenditures are 
%%%%%%%%%% free in the HDHP plan, PPO1200, but not in other plans
%%%%%%%%%% Numbers in paper to subtract are average preventive expenditures
%%%%%%%%%% Per person in plan

for i = 1:nIs
    for j =1:Sim
         HDHPC1(i,j) = max(0,HDHPC1(i,j)- 345); 
         HDHPC2(i,j) = max(0,HDHPC2(i,j)- 345); 
         HDHPC3(i,j) = max(0,HDHPC3(i,j)- 345);  
    end
end

%%%%%%%%%%%%%%% Generate OOP draws for Plan 3 - PPO1200 - YEAR 1

for i = 1:nIs
    for j = 1:Sim
        
        if income1(i,1) == 1 
            if HDHPC1(i,j) < 5200 && HDHPC1(i,j) >=1200 && MentalHealth1(i,j) ==0
               OutofPocket1(i,j,3) = 1200 + (HDHPC1(i,j)-1200)*0.2;
            end
            if HDHPC1(i,j) >= 5200 && MentalHealth1(i,j) ==0
               OutofPocket1(i,j,3) = 2000;
            end
            if MentalHealth1(i,j) > 0 && HDHPC1(i,j) < 5200 && (HDHPC1(i,j) + MentalHealth1(i,j)) >=1200
               OutofPocket1(i,j,3) = 1200 + (HDHPC1(i,j)-1200*(1-MentalHealth1(i,j)/(MentalHealth1(i,j) +HDHPC1(i,j) +0.001)))*0.2 + (MentalHealth1(i,j)-1200*MentalHealth1(i,j)/(MentalHealth1(i,j) +HDHPC1(i,j) +0.001))*0.5;  
            end
            if MentalHealth1(i,j) > 0 && HDHPC1(i,j) >= 5200 
               OutofPocket1(i,j,3) = 2000 + (MentalHealth1(i,j)-1200*MentalHealth1(i,j)/(MentalHealth1(i,j) +HDHPC1(i,j) +0.001))*0.5;  
            end  
        end
        
        if (income1(i,1) == 2 | income1(i,1)==3) 
            if HDHPC1(i,j) < 15200 && HDHPC1(i,j) >=1200 && MentalHealth1(i,j) ==0
               OutofPocket1(i,j,3) = 1200 + (HDHPC1(i,j)-1200)*0.2;
            end
            if HDHPC1(i,j) >= 15200 && MentalHealth1(i,j) ==0
               OutofPocket1(i,j,3) = 4000;
            end
            if MentalHealth1(i,j) > 0 && HDHPC1(i,j) < 15200 && (HDHPC1(i,j) + MentalHealth1(i,j)) >=1200
               OutofPocket1(i,j,3) = 1200 + (HDHPC1(i,j)-1200*(1-MentalHealth1(i,j)/(MentalHealth1(i,j) +HDHPC1(i,j) +0.001)))*0.2 + (MentalHealth1(i,j)-1200*MentalHealth1(i,j)/(MentalHealth1(i,j) +HDHPC1(i,j) +0.001))*0.5;  
            end
            if MentalHealth1(i,j) > 0 && HDHPC1(i,j) >= 15200 
               OutofPocket1(i,j,3) = 4000 + (MentalHealth1(i,j)-1200*MentalHealth1(i,j)/(MentalHealth1(i,j) +HDHPC1(i,j) +0.001))*0.5;  
            end  
        end        

        if (income1(i,1) == 4 | income1(i,1)==5) 
            if HDHPC1(i,j) < 20200 && HDHPC1(i,j) >=1200 && MentalHealth1(i,j) ==0
               OutofPocket1(i,j,3) = 1200 + (HDHPC1(i,j)-1200)*0.2;
            end
            if HDHPC1(i,j) >= 20200 && MentalHealth1(i,j) ==0
               OutofPocket1(i,j,3) = 5000;
            end
            if MentalHealth1(i,j) > 0 && HDHPC1(i,j) < 20200 && (HDHPC1(i,j) + MentalHealth1(i,j)) >=1200
               OutofPocket1(i,j,3) = 1200 + (HDHPC1(i,j)-1200*(1-MentalHealth1(i,j)/(MentalHealth1(i,j) +HDHPC1(i,j) +0.001)))*0.2 + (MentalHealth1(i,j)-1200*MentalHealth1(i,j)/(MentalHealth1(i,j) +HDHPC1(i,j) +0.001))*0.5;  
            end
            if MentalHealth1(i,j) > 0 && HDHPC1(i,j) >= 20200 
               OutofPocket1(i,j,3) = 5000 + (MentalHealth1(i,j)-1200*MentalHealth1(i,j)/(MentalHealth1(i,j) +HDHPC1(i,j) +0.001))*0.5;  
            end  
        end        
        
        if HDHPC1(i,j) < 1200 && MentalHealth1(i,j)==0
               OutofPocket1(i,j,3) = HDHPC1(i,j);
        end 
        
        if (HDHPC1(i,j) + MentalHealth1(i,j) <=1200)
               OutofPocket1(i,j,3) = HDHPC1(i,j) + MentalHealth1(i,j);
        end
    
    end
end

%%%%%%%%%%%%%%% Generate OOP draws for Plan 3 - PPO1200 - YEAR 2

for i = 1:nIs
    for j = 1:Sim
        
        if income2(i,1) == 1 
            if HDHPC2(i,j) < 5200 && HDHPC2(i,j) >=1200 && MentalHealth2(i,j) ==0
               OutofPocket2(i,j,3) = 1200 + (HDHPC2(i,j)-1200)*0.2;
            end
            if HDHPC2(i,j) >= 5200 && MentalHealth2(i,j) ==0
               OutofPocket2(i,j,3) = 2000;
            end
            if MentalHealth2(i,j) > 0 && HDHPC2(i,j) < 5200 && (HDHPC2(i,j) + MentalHealth2(i,j)) >=1200
               OutofPocket2(i,j,3) = 1200 + (HDHPC2(i,j)-1200*(1-MentalHealth2(i,j)/(MentalHealth2(i,j) +HDHPC2(i,j) +0.001)))*0.2 + (MentalHealth2(i,j)-1200*MentalHealth2(i,j)/(MentalHealth2(i,j) +HDHPC2(i,j) +0.001))*0.5;  
            end
            if MentalHealth2(i,j) > 0 && HDHPC2(i,j) >= 5200 
               OutofPocket2(i,j,3) = 2000 + (MentalHealth2(i,j)-1200*MentalHealth2(i,j)/(MentalHealth2(i,j) +HDHPC2(i,j) +0.001))*0.5;  
            end  
        end
        
        if (income2(i,1) == 2 | income2(i,1)==3) 
            if HDHPC2(i,j) < 15200 && HDHPC2(i,j) >=1200 && MentalHealth2(i,j) ==0
               OutofPocket2(i,j,3) = 1200 + (HDHPC2(i,j)-1200)*0.2;
            end
            if HDHPC2(i,j) >= 15200 && MentalHealth2(i,j) ==0
               OutofPocket2(i,j,3) = 4000;
            end
            if MentalHealth2(i,j) > 0 && HDHPC2(i,j) < 15200 && (HDHPC2(i,j) + MentalHealth2(i,j)) >=1200
               OutofPocket2(i,j,3) = 1200 + (HDHPC2(i,j)-1200*(1-MentalHealth2(i,j)/(MentalHealth2(i,j) +HDHPC2(i,j) +0.001)))*0.2 + (MentalHealth2(i,j)-1200*MentalHealth2(i,j)/(MentalHealth2(i,j) +HDHPC2(i,j) +0.001))*0.5;  
            end
            if MentalHealth2(i,j) > 0 && HDHPC2(i,j) >= 15200 
               OutofPocket2(i,j,3) = 4000 + (MentalHealth2(i,j)-1200*MentalHealth2(i,j)/(MentalHealth2(i,j) +HDHPC2(i,j) +0.001))*0.5;  
            end  
        end        

        if (income2(i,1) == 4 | income2(i,1)==5) 
            if HDHPC2(i,j) < 20200 && HDHPC2(i,j) >=1200 && MentalHealth2(i,j) ==0
               OutofPocket2(i,j,3) = 1200 + (HDHPC2(i,j)-1200)*0.2;
            end
            if HDHPC2(i,j) >= 20200 && MentalHealth2(i,j) ==0
               OutofPocket2(i,j,3) = 5000;
            end
            if MentalHealth2(i,j) > 0 && HDHPC2(i,j) < 20200 && (HDHPC2(i,j) + MentalHealth2(i,j)) >=1200
               OutofPocket2(i,j,3) = 1200 + (HDHPC2(i,j)-1200*(1-MentalHealth2(i,j)/(MentalHealth2(i,j) +HDHPC2(i,j) +0.001)))*0.2 + (MentalHealth2(i,j)-1200*MentalHealth2(i,j)/(MentalHealth2(i,j) +HDHPC2(i,j) +0.001))*0.5;  
            end
            if MentalHealth2(i,j) > 0 && HDHPC2(i,j) >= 20200 
               OutofPocket2(i,j,3) = 5000 + (MentalHealth2(i,j)-1200*MentalHealth2(i,j)/(MentalHealth2(i,j) +HDHPC2(i,j) +0.001))*0.5;  
            end  
        end        
        
        if HDHPC2(i,j) < 1200 && MentalHealth2(i,j)==0
               OutofPocket2(i,j,3) = HDHPC2(i,j);
        end 
        
        if (HDHPC2(i,j) + MentalHealth2(i,j) <=1200)
               OutofPocket2(i,j,3) = HDHPC2(i,j) + MentalHealth2(i,j);
        end
    
    end
end

%%%%%%%%%%%%%%% Generate OOP draws for Plan 3 - PPO1200 - YEAR 3

%%%%%%%%%%%%%%% NOTE: Deductible change in this plan and year from 1200 to 1400

for i = 1:nIs
    for j = 1:Sim
        
        if income2(i,1) == 1 
            if HDHPC3(i,j) < 5200 && HDHPC3(i,j) >=1400 && MentalHealth3(i,j) ==0
               OutofPocket3(i,j,3) = 1400 + (HDHPC3(i,j)-1400)*0.2;
            end
            if HDHPC3(i,j) >= 5200 && MentalHealth3(i,j) ==0
               OutofPocket3(i,j,3) = 2000;
            end
            if MentalHealth3(i,j) > 0 && HDHPC3(i,j) < 5200 && (HDHPC3(i,j) + MentalHealth3(i,j)) >=1400
               OutofPocket3(i,j,3) = 1400 + (HDHPC3(i,j)-1400*(1-MentalHealth3(i,j)/(MentalHealth3(i,j) +HDHPC3(i,j) +0.001)))*0.2 + (MentalHealth3(i,j)-1400*MentalHealth3(i,j)/(MentalHealth3(i,j) +HDHPC3(i,j) +0.001))*0.5;  
            end
            if MentalHealth3(i,j) > 0 && HDHPC3(i,j) >= 5200 
               OutofPocket3(i,j,3) = 2000 + (MentalHealth3(i,j)-1400*MentalHealth3(i,j)/(MentalHealth3(i,j) +HDHPC3(i,j) +0.001))*0.5;  
            end  
        end
        
        if (income2(i,1) == 2 | income2(i,1)==3) 
            if HDHPC3(i,j) < 15200 && HDHPC3(i,j) >=1400 && MentalHealth3(i,j) ==0
               OutofPocket3(i,j,3) = 1400 + (HDHPC3(i,j)-1400)*0.2;
            end
            if HDHPC3(i,j) >= 15200 && MentalHealth3(i,j) ==0
               OutofPocket3(i,j,3) = 4000;
            end
            if MentalHealth3(i,j) > 0 && HDHPC3(i,j) < 15200 && (HDHPC3(i,j) + MentalHealth3(i,j)) >=1400
               OutofPocket3(i,j,3) = 1400 + (HDHPC3(i,j)-1400*(1-MentalHealth3(i,j)/(MentalHealth3(i,j) +HDHPC3(i,j) +0.001)))*0.2 + (MentalHealth3(i,j)-1400*MentalHealth3(i,j)/(MentalHealth3(i,j) +HDHPC3(i,j) +0.001))*0.5;  
            end
            if MentalHealth3(i,j) > 0 && HDHPC3(i,j) >= 15200 
               OutofPocket3(i,j,3) = 4000 + (MentalHealth3(i,j)-1400*MentalHealth3(i,j)/(MentalHealth3(i,j) +HDHPC3(i,j) +0.001))*0.5;  
            end  
        end        

        if (income2(i,1) == 4 | income2(i,1)==5) 
            if HDHPC3(i,j) < 20200 && HDHPC3(i,j) >=1400 && MentalHealth3(i,j) ==0
               OutofPocket3(i,j,3) = 1400 + (HDHPC3(i,j)-1400)*0.2;
            end
            if HDHPC3(i,j) >= 20200 && MentalHealth3(i,j) ==0
               OutofPocket3(i,j,3) = 5000;
            end
            if MentalHealth3(i,j) > 0 && HDHPC3(i,j) < 20200 && (HDHPC3(i,j) + MentalHealth3(i,j)) >=1400
               OutofPocket3(i,j,3) = 1400 + (HDHPC3(i,j)-1400*(1-MentalHealth3(i,j)/(MentalHealth3(i,j) +HDHPC3(i,j) +0.001)))*0.2 + (MentalHealth3(i,j)-1400*MentalHealth3(i,j)/(MentalHealth3(i,j) +HDHPC3(i,j) +0.001))*0.5;  
            end
            if MentalHealth3(i,j) > 0 && HDHPC3(i,j) >= 20200 
               OutofPocket3(i,j,3) = 5000 + (MentalHealth3(i,j)-1400*MentalHealth3(i,j)/(MentalHealth3(i,j) +HDHPC3(i,j) +0.001))*0.5;  
            end  
        end        
        
        if HDHPC3(i,j) < 1400 && MentalHealth3(i,j)==0
               OutofPocket3(i,j,3) = HDHPC3(i,j);
        end 
        
        if (HDHPC3(i,j) + MentalHealth3(i,j) <=1400)
               OutofPocket3(i,j,3) = HDHPC3(i,j) + MentalHealth3(i,j);
        end     
    
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Now, for each plan, individual, year and expenditure draw we have 
%%%%%%%%%%%%%%%%%%%% OOP expenditure under each plan. the enxt part of the code aggregates
%%%%%%%%%%%%%%%%%%%% individuals into family units, which is the relevant unit for decision 
%%%%%%%%%%%%%%%%%%%% making. The plans have family deductible caps and family OOP max caps which 
%%%%%%%%%%%%%%%%%%%% bind if multiple family members reach their individual deductibles and 
%%%%%%%%%%%%%%%%%%%% individual OOP maxes. These are the family restrictions that bind and require 
%%%%%%%%%%%%%%%%%%%% us to do more than just adding up individual OOP expenditure draws at the family level
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% For the plans here, as discussed in Table 2, PPO1200 incorporates Family restrictions if 
%%%%%%%%%%%%%%%%%%%% 2 people in family hit deductible / OOP max. In other plans family caps start to bind 
%%%%%%%%%%%%%%%%%%%% when 3 people hit these individual limits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K=7000;

FamilyOOP1 = zeros(K,Sim,3);
FamilyOOP2 = zeros(K,Sim,3);
FamilyOOP3 = zeros(K,Sim,3);

FamilyHosp1 = zeros(K,Sim,3);
FamilyHosp2 = zeros(K,Sim,3);
FamilyHosp3 = zeros(K,Sim,3);

FamilyDed1 = zeros(K,Sim,3);
FamilyDed2 = zeros(K,Sim,3);
FamilyDed3 = zeros(K,Sim,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Following loop (approximately 400 lines of code) goes through population,
%%%% identifies people as being in families, and incrementally adds people in given 
%%%% family to incremental expenditures. Part of incremental process is to determine if 
%%%% family has reached deductible or OOP max FAMILY cap by having many individuals 
%%%% hit their individual caps for given plans. 
%%%%
%%%% End result is FAMILY-PLAN-YEAR specific distribution of OOP expenditures, which
%%%% is then main input into plan choice models. 
%%%%
%%%% Note: One reason this loop is lengthy is because it has to be run for different 
%%%% income tier groups because OOP maxes differ by income tier groups for each plan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IND = zeros(K,1);

for i = 1:nIs
    for j = 1:Sim
        
        
        if j ==1
            IND(famid(i,1),1)= IND(famid(i,1),1) + 1;
        end
    
        FamilyDed1(famid(i,1),j,1) = FamilyDed1(famid(i,1),j,1) + min(250,Hosp1(i,j))*(Hosp1(i,j)/(Hosp1(i,j) + MentalHealth1(i,j)+0.0001)) + min(250,MentalHealth1(i,j))*(MentalHealth1(i,j)/(Hosp1(i,j) + MentalHealth1(i,j)+ 0.0001)) ;
        FamilyDed1(famid(i,1),j,2) = FamilyDed1(famid(i,1),j,2) + min(500,Hosp1(i,j))*(Hosp1(i,j)/(Hosp1(i,j) + MentalHealth1(i,j)+0.0001)) + min(500,MentalHealth1(i,j))*(MentalHealth1(i,j)/(Hosp1(i,j) + MentalHealth1(i,j)+ 0.0001)) ;
        FamilyDed1(famid(i,1),j,3) = FamilyDed1(famid(i,1),j,3) + min(1200,HDHPC1(i,j))*(HDHPC1(i,j)/(HDHPC1(i,j) + MentalHealth1(i,j)+0.0001)) + min(1200,MentalHealth1(i,j))*(MentalHealth1(i,j)/(HDHPC1(i,j) + MentalHealth1(i,j)+ 0.0001)) ;
        
        FamilyDed2(famid(i,1),j,1) = FamilyDed2(famid(i,1),j,1) + min(250,Hosp2(i,j))*(Hosp2(i,j)/(Hosp2(i,j) + MentalHealth2(i,j)+0.0001)) + min(250,MentalHealth2(i,j))*(MentalHealth2(i,j)/(Hosp2(i,j) + MentalHealth2(i,j)+ 0.0001)) ;
        FamilyDed2(famid(i,1),j,2) = FamilyDed2(famid(i,1),j,2) + min(500,Hosp2(i,j))*(Hosp2(i,j)/(Hosp2(i,j) + MentalHealth2(i,j)+0.0001)) + min(500,MentalHealth2(i,j))*(MentalHealth2(i,j)/(Hosp2(i,j) + MentalHealth2(i,j)+ 0.0001)) ;
        FamilyDed2(famid(i,1),j,3) = FamilyDed2(famid(i,1),j,3) + min(1200,HDHPC2(i,j))*(HDHPC2(i,j)/(HDHPC2(i,j) + MentalHealth2(i,j)+0.0001)) + min(1200,MentalHealth2(i,j))*(MentalHealth2(i,j)/(HDHPC2(i,j) + MentalHealth2(i,j)+ 0.0001)) ;
        
        FamilyDed3(famid(i,1),j,1) = FamilyDed3(famid(i,1),j,1) + min(250,Hosp3(i,j))*(Hosp3(i,j)/(Hosp3(i,j) + MentalHealth3(i,j)+0.0001)) + min(250,MentalHealth3(i,j))*(MentalHealth3(i,j)/(Hosp3(i,j) + MentalHealth3(i,j)+ 0.0001)) ;
        FamilyDed3(famid(i,1),j,2) = FamilyDed3(famid(i,1),j,2) + min(500,Hosp3(i,j))*(Hosp3(i,j)/(Hosp3(i,j) + MentalHealth3(i,j)+0.0001)) + min(500,MentalHealth3(i,j))*(MentalHealth3(i,j)/(Hosp3(i,j) + MentalHealth3(i,j)+ 0.0001)) ;
        FamilyDed3(famid(i,1),j,3) = FamilyDed3(famid(i,1),j,3) + min(1400,HDHPC3(i,j))*(HDHPC3(i,j)/(HDHPC3(i,j) + MentalHealth3(i,j)+0.0001)) + min(1400,MentalHealth3(i,j))*(MentalHealth3(i,j)/(HDHPC3(i,j) + MentalHealth3(i,j)+ 0.0001)) ;
               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        FamilyHosp1(famid(i,1),j,1) = FamilyHosp1(famid(i,1),j,1) + OutofPocket1(i,j,1) - OVoopCP1(i,j) - RXoopCP1(i,j) - (MentalHealth1(i,j)- min(250,MentalHealth1(i,j))*(MentalHealth1(i,j)/(Hosp1(i,j) + MentalHealth1(i,j) + 0.0001)))*.5;
        FamilyHosp1(famid(i,1),j,2) = FamilyHosp1(famid(i,1),j,2) + OutofPocket1(i,j,2) - OVoopCP1(i,j) - RXoopCP1(i,j) - (MentalHealth1(i,j)- min(500,MentalHealth1(i,j))*(MentalHealth1(i,j)/(Hosp1(i,j) + MentalHealth1(i,j) + 0.0001)))*.5;
        FamilyHosp1(famid(i,1),j,3) = FamilyHosp1(famid(i,1),j,3) + OutofPocket1(i,j,3) - (MentalHealth1(i,j)- min(1200,MentalHealth1(i,j))*(MentalHealth1(i,j)/(Hosp1(i,j) + MentalHealth1(i,j) + 0.0001)))*.5;
        
        FamilyHosp2(famid(i,1),j,1) = FamilyHosp2(famid(i,1),j,1) + OutofPocket2(i,j,1) - OVoopCP2(i,j) - RXoopCP2(i,j) - (MentalHealth2(i,j)- min(250,MentalHealth2(i,j))*(MentalHealth2(i,j)/(Hosp2(i,j) + MentalHealth2(i,j) + 0.0001)))*.5;
        FamilyHosp2(famid(i,1),j,2) = FamilyHosp2(famid(i,1),j,2) + OutofPocket2(i,j,2) - OVoopCP2(i,j) - RXoopCP2(i,j) - (MentalHealth2(i,j)- min(500,MentalHealth2(i,j))*(MentalHealth2(i,j)/(Hosp2(i,j) + MentalHealth2(i,j) + 0.0001)))*.5;
        FamilyHosp2(famid(i,1),j,3) = FamilyHosp2(famid(i,1),j,3) + OutofPocket2(i,j,3) - (MentalHealth2(i,j)- min(1200,MentalHealth2(i,j))*(MentalHealth2(i,j)/(Hosp2(i,j) + MentalHealth2(i,j) + 0.0001)))*.5;
        
        FamilyHosp3(famid(i,1),j,1) = FamilyHosp3(famid(i,1),j,1) + OutofPocket3(i,j,1) - OVoopCP3(i,j) - RXoopCP3(i,j) - (MentalHealth3(i,j)- min(250,MentalHealth3(i,j))*(MentalHealth3(i,j)/(Hosp3(i,j) + MentalHealth3(i,j) + 0.0001)))*.5;
        FamilyHosp3(famid(i,1),j,2) = FamilyHosp3(famid(i,1),j,2) + OutofPocket3(i,j,2) - OVoopCP3(i,j) - RXoopCP3(i,j) - (MentalHealth3(i,j)- min(500,MentalHealth3(i,j))*(MentalHealth3(i,j)/(Hosp3(i,j) + MentalHealth3(i,j) + 0.0001)))*.5;
        FamilyHosp3(famid(i,1),j,3) = FamilyHosp3(famid(i,1),j,3) + OutofPocket3(i,j,3) - (MentalHealth3(i,j)- min(1200,MentalHealth3(i,j))*(MentalHealth3(i,j)/(Hosp3(i,j) + MentalHealth3(i,j) + 0.0001)))*.5;
        
        if IND(famid(i,1),1) <=2   
        
        FamilyOOP1(famid(i,1),j,1) = FamilyOOP1(famid(i,1),j,1) + OutofPocket1(i,j,1);
        FamilyOOP1(famid(i,1),j,2) = FamilyOOP1(famid(i,1),j,2) + OutofPocket1(i,j,2);
        FamilyOOP1(famid(i,1),j,3) = FamilyOOP1(famid(i,1),j,3) + OutofPocket1(i,j,3);
        
        FamilyOOP2(famid(i,1),j,1) = FamilyOOP2(famid(i,1),j,1) + OutofPocket2(i,j,1);
        FamilyOOP2(famid(i,1),j,2) = FamilyOOP2(famid(i,1),j,2) + OutofPocket2(i,j,2);
        FamilyOOP2(famid(i,1),j,3) = FamilyOOP2(famid(i,1),j,3) + OutofPocket2(i,j,3);
        
        FamilyOOP3(famid(i,1),j,1) = FamilyOOP3(famid(i,1),j,1) + OutofPocket3(i,j,1);
        FamilyOOP3(famid(i,1),j,2) = FamilyOOP3(famid(i,1),j,2) + OutofPocket3(i,j,2);
        FamilyOOP3(famid(i,1),j,3) = FamilyOOP3(famid(i,1),j,3) + OutofPocket3(i,j,3);
    
        end
        
        if IND(famid(i,1),1) >=3   
        
            if (income1(i,1) ==1)
                 
                if FamilyHosp1(famid(i,1),j,1) > 3000
                    FamilyOOP1(famid(i,1),j,1) = FamilyOOP1(famid(i,1),j,1) + OutofPocket1(i,j,1) - (FamilyHosp1(famid(i,1),j,1)-3000);
                    FamilyHosp1(famid(i,1),j,1) = 3000;
                 end
                 
                 if FamilyHosp1(famid(i,1),j,1) <= 3000 && FamilyDed1(famid(i,1),j,1) < 750
                     FamilyOOP1(famid(i,1),j,1) = FamilyOOP1(famid(i,1),j,1) + OutofPocket1(i,j,1);
                 end
       
                 if FamilyHosp1(famid(i,1),j,1) < 3000 && FamilyDed1(famid(i,1),j,1) >= 750
                     FamilyOOP1(famid(i,1),j,1) = FamilyOOP1(famid(i,1),j,1) + OutofPocket1(i,j,1) - 0.9*(FamilyDed1(famid(i,1),j,1)-750);
                     FamilyDed1(famid(i,1),j,1)= 750;
                 end
            end
            
            if (income2(i,1) ==1)
                if FamilyHosp2(famid(i,1),j,1) > 3000
                    FamilyOOP2(famid(i,1),j,1) = FamilyOOP2(famid(i,1),j,1) + OutofPocket2(i,j,1) - (FamilyHosp2(famid(i,1),j,1)-3000);
                    FamilyHosp2(famid(i,1),j,1) = 3000;
                 end
                 
                 if FamilyHosp2(famid(i,1),j,1) <= 3000 && FamilyDed2(famid(i,1),j,1) < 750
                     FamilyOOP2(famid(i,1),j,1) = FamilyOOP2(famid(i,1),j,1) + OutofPocket2(i,j,1);
                 end
       
                 if FamilyHosp2(famid(i,1),j,1) < 3000 && FamilyDed2(famid(i,1),j,1) >= 750
                     FamilyOOP2(famid(i,1),j,1) = FamilyOOP2(famid(i,1),j,1) + OutofPocket2(i,j,1) - 0.9*(FamilyDed2(famid(i,1),j,1)-750);
                     FamilyDed2(famid(i,1),j,1)= 750;
                 end
              
                 
                if FamilyHosp3(famid(i,1),j,1) > 3000
                    FamilyOOP3(famid(i,1),j,1) = FamilyOOP3(famid(i,1),j,1) + OutofPocket3(i,j,1) - (FamilyHosp3(famid(i,1),j,1)-3000);
                    FamilyHosp3(famid(i,1),j,1) = 3000;
                 end
                 
                 if FamilyHosp3(famid(i,1),j,1) <= 3000 && FamilyDed3(famid(i,1),j,1) < 750
                     FamilyOOP3(famid(i,1),j,1) = FamilyOOP3(famid(i,1),j,1) + OutofPocket3(i,j,1);
                 end
       
                 if FamilyHosp3(famid(i,1),j,1) < 3000 && FamilyDed3(famid(i,1),j,1) >= 750
                     FamilyOOP3(famid(i,1),j,1) = FamilyOOP3(famid(i,1),j,1) + OutofPocket3(i,j,1) - 0.9*(FamilyDed3(famid(i,1),j,1)-750);
                     FamilyDed3(famid(i,1),j,1)= 750;
                 end
            end
            if (income1(i,1) ==1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
                 
                 if FamilyHosp1(famid(i,1),j,2) > 4500
                    FamilyOOP1(famid(i,1),j,2) = FamilyOOP1(famid(i,1),j,2) + OutofPocket1(i,j,2) - (FamilyHosp1(famid(i,1),j,2)-4500);
                    FamilyHosp1(famid(i,1),j,2)=4500;
                 end
                 
                 if FamilyHosp1(famid(i,1),j,2) <= 4500 && FamilyDed1(famid(i,1),j,2) < 1500
                     FamilyOOP1(famid(i,1),j,2) = FamilyOOP1(famid(i,1),j,2) + OutofPocket1(i,j,2);
                 end
       
                 if FamilyHosp1(famid(i,1),j,2) < 4500 && FamilyDed1(famid(i,1),j,2) >= 1500
                     FamilyOOP1(famid(i,1),j,2) = FamilyOOP1(famid(i,1),j,2) + OutofPocket1(i,j,2) - 0.8*(FamilyDed1(famid(i,1),j,2)-1500);
                     FamilyDed1(famid(i,1),j,2)= 1500;
                 end
            end
            if (income2(i,1) ==1)
                 if FamilyHosp2(famid(i,1),j,2) > 4500
                    FamilyOOP2(famid(i,1),j,2) = FamilyOOP2(famid(i,1),j,2) + OutofPocket2(i,j,2) - (FamilyHosp2(famid(i,1),j,2)-4500);
                    FamilyHosp2(famid(i,1),j,2)=4500;
                 end
            
                 if FamilyHosp2(famid(i,1),j,2) <= 4500 && FamilyDed2(famid(i,1),j,2) < 1500
                     FamilyOOP2(famid(i,1),j,2) = FamilyOOP2(famid(i,1),j,2) + OutofPocket2(i,j,2);
                 end
       
                 if FamilyHosp2(famid(i,1),j,2) < 4500 && FamilyDed2(famid(i,1),j,2) >= 1500
                     FamilyOOP2(famid(i,1),j,2) = FamilyOOP2(famid(i,1),j,2) + OutofPocket2(i,j,2) - 0.8*(FamilyDed2(famid(i,1),j,2)-1500);
                     FamilyDed2(famid(i,1),j,2)= 1500;
                 end

                 
                 if FamilyHosp3(famid(i,1),j,2) > 4500
                    FamilyOOP3(famid(i,1),j,2) = FamilyOOP3(famid(i,1),j,2) + OutofPocket3(i,j,2) - (FamilyHosp3(famid(i,1),j,2)-4500);
                    FamilyHosp3(famid(i,1),j,2)=4500;
                 end
                 
                 if FamilyHosp3(famid(i,1),j,2) <= 4500 && FamilyDed3(famid(i,1),j,2) < 1500
                     FamilyOOP3(famid(i,1),j,2) = FamilyOOP3(famid(i,1),j,2) + OutofPocket3(i,j,2);
                 end
       
                 if FamilyHosp3(famid(i,1),j,2) < 4500 && FamilyDed3(famid(i,1),j,2) >= 1500
                     FamilyOOP3(famid(i,1),j,2) = FamilyOOP3(famid(i,1),j,2) + OutofPocket3(i,j,2) - 0.8*(FamilyDed3(famid(i,1),j,2)-1500);
                     FamilyDed3(famid(i,1),j,2)= 1500;
                 end
            end
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (income1(i,1) ==1)        
                     
                               
                 if FamilyHosp1(famid(i,1),j,3) >= 6000         
                    FamilyOOP1(famid(i,1),j,3) = FamilyOOP1(famid(i,1),j,3) + OutofPocket1(i,j,3)- (FamilyHosp1(famid(i,1),j,3)-6000);
                    FamilyHosp1(famid(i,1),j,3)= 6000;
                 end
                 
                 if FamilyHosp1(famid(i,1),j,3) < 6000 && FamilyDed1(famid(i,1),j,3) < 2400
                     FamilyOOP1(famid(i,1),j,3) = FamilyOOP1(famid(i,1),j,3) + OutofPocket1(i,j,3);
                 end
       
                 if FamilyHosp1(famid(i,1),j,3) < 6000 && FamilyDed1(famid(i,1),j,3) >= 2400
                     FamilyOOP1(famid(i,1),j,3) = FamilyOOP1(famid(i,1),j,3) + OutofPocket1(i,j,3) - 0.8*(FamilyDed1(famid(i,1),j,3)-2400);
                     FamilyDed1(famid(i,1),j,3)=2400;
                 end
            end
            if (income2(i,1) ==1)
                 if FamilyHosp2(famid(i,1),j,3) >= 6000         
                    FamilyOOP2(famid(i,1),j,3) = FamilyOOP2(famid(i,1),j,3) + OutofPocket2(i,j,3)- (FamilyHosp2(famid(i,1),j,3)-6000);
                    FamilyHosp2(famid(i,1),j,3)= 6000;
                 end
                 
                 if FamilyHosp2(famid(i,1),j,3) < 6000 && FamilyDed2(famid(i,1),j,3) < 2400
                     FamilyOOP2(famid(i,1),j,3) = FamilyOOP2(famid(i,1),j,3) + OutofPocket2(i,j,3);
                 end
       
                 if FamilyHosp2(famid(i,1),j,3) < 6000 && FamilyDed2(famid(i,1),j,3) >= 2400
                     FamilyOOP2(famid(i,1),j,3) = FamilyOOP2(famid(i,1),j,3) + OutofPocket2(i,j,3) - 0.8*(FamilyDed2(famid(i,1),j,3)-2400);
                     FamilyDed2(famid(i,1),j,3)=2400;
                 end
                 
                 if FamilyHosp3(famid(i,1),j,3) >= 6000         
                    FamilyOOP3(famid(i,1),j,3) = FamilyOOP3(famid(i,1),j,3) + OutofPocket3(i,j,3)- (FamilyHosp3(famid(i,1),j,3)-6000);
                    FamilyHosp3(famid(i,1),j,3)= 6000;
                 end
                 
                 if FamilyHosp3(famid(i,1),j,3) < 6000 && FamilyDed3(famid(i,1),j,3) < 2800
                     FamilyOOP3(famid(i,1),j,3) = FamilyOOP3(famid(i,1),j,3) + OutofPocket3(i,j,3);
                 end
       
                 if FamilyHosp3(famid(i,1),j,3) < 6000 && FamilyDed3(famid(i,1),j,3) >= 2800
                     FamilyOOP3(famid(i,1),j,3) = FamilyOOP3(famid(i,1),j,3) + OutofPocket3(i,j,3) - 0.8*(FamilyDed3(famid(i,1),j,3)-2800);
                     FamilyDed3(famid(i,1),j,3)=2400;
                 end
            
            end   
                   
            if (income1(i,1) ==2 | income1(i,1) ==3)
                 
                if FamilyHosp1(famid(i,1),j,1) > 5000
                    FamilyOOP1(famid(i,1),j,1) = FamilyOOP1(famid(i,1),j,1) + OutofPocket1(i,j,1) - (FamilyHosp1(famid(i,1),j,1)-5000);
                    FamilyHosp1(famid(i,1),j,1) = 5000;
                 end
                 
                 if FamilyHosp1(famid(i,1),j,1) <= 5000 && FamilyDed1(famid(i,1),j,1) < 750
                     FamilyOOP1(famid(i,1),j,1) = FamilyOOP1(famid(i,1),j,1) + OutofPocket1(i,j,1);
                 end
       
                 if FamilyHosp1(famid(i,1),j,1) < 5000 && FamilyDed1(famid(i,1),j,1) >= 750
                     FamilyOOP1(famid(i,1),j,1) = FamilyOOP1(famid(i,1),j,1) + OutofPocket1(i,j,1) - 0.9*(FamilyDed1(famid(i,1),j,1)-750);
                     FamilyDed1(famid(i,1),j,1)= 750;
                 end
            end
            
            if (income2(i,1) ==2 | income2(i,1) ==3)
                if FamilyHosp2(famid(i,1),j,1) > 5000
                    FamilyOOP2(famid(i,1),j,1) = FamilyOOP2(famid(i,1),j,1) + OutofPocket2(i,j,1) - (FamilyHosp2(famid(i,1),j,1)-5000);
                    FamilyHosp2(famid(i,1),j,1) = 5000;
                 end
                 
                 if FamilyHosp2(famid(i,1),j,1) <= 5000 && FamilyDed2(famid(i,1),j,1) < 750
                     FamilyOOP2(famid(i,1),j,1) = FamilyOOP2(famid(i,1),j,1) + OutofPocket2(i,j,1);
                 end
       
                 if FamilyHosp2(famid(i,1),j,1) < 5000 && FamilyDed2(famid(i,1),j,1) >= 750
                     FamilyOOP2(famid(i,1),j,1) = FamilyOOP2(famid(i,1),j,1) + OutofPocket2(i,j,1) - 0.9*(FamilyDed2(famid(i,1),j,1)-750);
                     FamilyDed2(famid(i,1),j,1)= 750;
                 end
              
                 
                if FamilyHosp3(famid(i,1),j,1) > 5000
                    FamilyOOP3(famid(i,1),j,1) = FamilyOOP3(famid(i,1),j,1) + OutofPocket3(i,j,1) - (FamilyHosp3(famid(i,1),j,1)-5000);
                    FamilyHosp3(famid(i,1),j,1) = 5000;
                 end
                 
                 if FamilyHosp3(famid(i,1),j,1) <= 5000 && FamilyDed3(famid(i,1),j,1) < 750
                     FamilyOOP3(famid(i,1),j,1) = FamilyOOP3(famid(i,1),j,1) + OutofPocket3(i,j,1);
                 end
       
                 if FamilyHosp3(famid(i,1),j,1) < 5000 && FamilyDed3(famid(i,1),j,1) >= 750
                     FamilyOOP3(famid(i,1),j,1) = FamilyOOP3(famid(i,1),j,1) + OutofPocket3(i,j,1) - 0.9*(FamilyDed3(famid(i,1),j,1)-750);
                     FamilyDed3(famid(i,1),j,1)= 750;
                 end
            end
            if (income1(i,1) ==2 | income1(i,1) ==3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
                 
                 if FamilyHosp1(famid(i,1),j,2) > 7000
                    FamilyOOP1(famid(i,1),j,2) = FamilyOOP1(famid(i,1),j,2) + OutofPocket1(i,j,2) - (FamilyHosp1(famid(i,1),j,2)-7000);
                    FamilyHosp1(famid(i,1),j,2)=7000;
                 end
                 
                 if FamilyHosp1(famid(i,1),j,2) <= 7000 && FamilyDed1(famid(i,1),j,2) < 1500
                     FamilyOOP1(famid(i,1),j,2) = FamilyOOP1(famid(i,1),j,2) + OutofPocket1(i,j,2);
                 end
       
                 if FamilyHosp1(famid(i,1),j,2) < 7000 && FamilyDed1(famid(i,1),j,2) >= 1500
                     FamilyOOP1(famid(i,1),j,2) = FamilyOOP1(famid(i,1),j,2) + OutofPocket1(i,j,2) - 0.8*(FamilyDed1(famid(i,1),j,2)-1500);
                     FamilyDed1(famid(i,1),j,2)= 1500;
                 end
            end
            if (income2(i,1) ==2 | income2(i,1) ==3)
                 if FamilyHosp2(famid(i,1),j,2) > 7000
                    FamilyOOP2(famid(i,1),j,2) = FamilyOOP2(famid(i,1),j,2) + OutofPocket2(i,j,2) - (FamilyHosp2(famid(i,1),j,2)-7000);
                    FamilyHosp2(famid(i,1),j,2)=7000;
                 end
            
                 if FamilyHosp2(famid(i,1),j,2) <= 7000 && FamilyDed2(famid(i,1),j,2) < 1500
                     FamilyOOP2(famid(i,1),j,2) = FamilyOOP2(famid(i,1),j,2) + OutofPocket2(i,j,2);
                 end
       
                 if FamilyHosp2(famid(i,1),j,2) < 7000 && FamilyDed2(famid(i,1),j,2) >= 1500
                     FamilyOOP2(famid(i,1),j,2) = FamilyOOP2(famid(i,1),j,2) + OutofPocket2(i,j,2) - 0.8*(FamilyDed2(famid(i,1),j,2)-1500);
                     FamilyDed2(famid(i,1),j,2)= 1500;
                 end

                 
                 if FamilyHosp3(famid(i,1),j,2) > 7000
                    FamilyOOP3(famid(i,1),j,2) = FamilyOOP3(famid(i,1),j,2) + OutofPocket3(i,j,2) - (FamilyHosp3(famid(i,1),j,2)-7000);
                    FamilyHosp3(famid(i,1),j,2)=7000;
                 end
                 
                 if FamilyHosp3(famid(i,1),j,2) <= 7000 && FamilyDed3(famid(i,1),j,2) < 1500
                     FamilyOOP3(famid(i,1),j,2) = FamilyOOP3(famid(i,1),j,2) + OutofPocket3(i,j,2);
                 end
       
                 if FamilyHosp3(famid(i,1),j,2) < 7000 && FamilyDed3(famid(i,1),j,2) >= 1500
                     FamilyOOP3(famid(i,1),j,2) = FamilyOOP3(famid(i,1),j,2) + OutofPocket3(i,j,2) - 0.8*(FamilyDed3(famid(i,1),j,2)-1500);
                     FamilyDed3(famid(i,1),j,2)= 1500;
                 end
            end
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (income1(i,1) ==2 | income1(i,1) ==3)        
                     
                               
                 if FamilyHosp1(famid(i,1),j,3) >= 8000         
                    FamilyOOP1(famid(i,1),j,3) = FamilyOOP1(famid(i,1),j,3) + OutofPocket1(i,j,3)- (FamilyHosp1(famid(i,1),j,3)-8000);
                    FamilyHosp1(famid(i,1),j,3)= 8000;
                 end
                 
                 if FamilyHosp1(famid(i,1),j,3) < 8000 && FamilyDed1(famid(i,1),j,3) < 2400
                     FamilyOOP1(famid(i,1),j,3) = FamilyOOP1(famid(i,1),j,3) + OutofPocket1(i,j,3);
                 end
       
                 if FamilyHosp1(famid(i,1),j,3) < 8000 && FamilyDed1(famid(i,1),j,3) >= 2400
                     FamilyOOP1(famid(i,1),j,3) = FamilyOOP1(famid(i,1),j,3) + OutofPocket1(i,j,3) - 0.8*(FamilyDed1(famid(i,1),j,3)-2400);
                     FamilyDed1(famid(i,1),j,3)=2400;
                 end
            end
            if (income2(i,1) ==2 | income2(i,1) ==3)
                 if FamilyHosp2(famid(i,1),j,3) >= 8000         
                    FamilyOOP2(famid(i,1),j,3) = FamilyOOP2(famid(i,1),j,3) + OutofPocket2(i,j,3)- (FamilyHosp2(famid(i,1),j,3)-8000);
                    FamilyHosp2(famid(i,1),j,3)= 8000;
                 end
                 
                 if FamilyHosp2(famid(i,1),j,3) < 8000 && FamilyDed2(famid(i,1),j,3) < 2400
                     FamilyOOP2(famid(i,1),j,3) = FamilyOOP2(famid(i,1),j,3) + OutofPocket2(i,j,3);
                 end
       
                 if FamilyHosp2(famid(i,1),j,3) < 8000 && FamilyDed2(famid(i,1),j,3) >= 2400
                     FamilyOOP2(famid(i,1),j,3) = FamilyOOP2(famid(i,1),j,3) + OutofPocket2(i,j,3) - 0.8*(FamilyDed2(famid(i,1),j,3)-2400);
                     FamilyDed2(famid(i,1),j,3)=2400;
                 end
                 
                 if FamilyHosp3(famid(i,1),j,3) >= 8000         
                    FamilyOOP3(famid(i,1),j,3) = FamilyOOP3(famid(i,1),j,3) + OutofPocket3(i,j,3)- (FamilyHosp3(famid(i,1),j,3)-8000);
                    FamilyHosp3(famid(i,1),j,3)= 8000;
                 end
                 
                 if FamilyHosp3(famid(i,1),j,3) < 8000 && FamilyDed3(famid(i,1),j,3) < 2800
                     FamilyOOP3(famid(i,1),j,3) = FamilyOOP3(famid(i,1),j,3) + OutofPocket3(i,j,3);
                 end
       
                 if FamilyHosp3(famid(i,1),j,3) < 8000 && FamilyDed3(famid(i,1),j,3) >= 2800
                     FamilyOOP3(famid(i,1),j,3) = FamilyOOP3(famid(i,1),j,3) + OutofPocket3(i,j,3) - 0.8*(FamilyDed3(famid(i,1),j,3)-2800);
                     FamilyDed3(famid(i,1),j,3)=2400;
                 end
            
            end           
                 
        
            if (income1(i,1) ==4 | income1(i,1) ==5)
                 
                if FamilyHosp1(famid(i,1),j,1) > 8000
                    FamilyOOP1(famid(i,1),j,1) = FamilyOOP1(famid(i,1),j,1) + OutofPocket1(i,j,1) - (FamilyHosp1(famid(i,1),j,1)-8000);
                    FamilyHosp1(famid(i,1),j,1) = 8000;
                 end
                 
                 if FamilyHosp1(famid(i,1),j,1) <= 8000 && FamilyDed1(famid(i,1),j,1) < 750
                     FamilyOOP1(famid(i,1),j,1) = FamilyOOP1(famid(i,1),j,1) + OutofPocket1(i,j,1);
                 end
       
                 if FamilyHosp1(famid(i,1),j,1) < 8000 && FamilyDed1(famid(i,1),j,1) >= 750
                     FamilyOOP1(famid(i,1),j,1) = FamilyOOP1(famid(i,1),j,1) + OutofPocket1(i,j,1) - 0.8*(FamilyDed1(famid(i,1),j,1)-750);
                     FamilyDed1(famid(i,1),j,1)= 750;
                 end
            end
            if (income2(i,1) ==4 | income2(i,1) ==5)
                if FamilyHosp2(famid(i,1),j,1) > 8000
                    FamilyOOP2(famid(i,1),j,1) = FamilyOOP2(famid(i,1),j,1) + OutofPocket2(i,j,1) - (FamilyHosp2(famid(i,1),j,1)-8000);
                    FamilyHosp2(famid(i,1),j,1) = 8000;
                 end
                 
                 if FamilyHosp2(famid(i,1),j,1) <= 8000 && FamilyDed2(famid(i,1),j,1) < 750
                     FamilyOOP2(famid(i,1),j,1) = FamilyOOP2(famid(i,1),j,1) + OutofPocket2(i,j,1);
                 end
       
                 if FamilyHosp2(famid(i,1),j,1) < 8000 && FamilyDed2(famid(i,1),j,1) >= 750
                     FamilyOOP2(famid(i,1),j,1) = FamilyOOP2(famid(i,1),j,1) + OutofPocket2(i,j,1) - 0.8*(FamilyDed2(famid(i,1),j,1)-750);
                     FamilyDed2(famid(i,1),j,1)= 750;
                 end
              
                 
                if FamilyHosp3(famid(i,1),j,1) > 8000
                    FamilyOOP3(famid(i,1),j,1) = FamilyOOP3(famid(i,1),j,1) + OutofPocket3(i,j,1) - (FamilyHosp3(famid(i,1),j,1)-8000);
                    FamilyHosp3(famid(i,1),j,1) = 8000;
                 end
                 
                 if FamilyHosp3(famid(i,1),j,1) <= 8000 && FamilyDed3(famid(i,1),j,1) < 750
                     FamilyOOP3(famid(i,1),j,1) = FamilyOOP3(famid(i,1),j,1) + OutofPocket3(i,j,1);
                 end
       
                 if FamilyHosp3(famid(i,1),j,1) < 8000 && FamilyDed3(famid(i,1),j,1) >= 750
                     FamilyOOP3(famid(i,1),j,1) = FamilyOOP3(famid(i,1),j,1) + OutofPocket3(i,j,1) - 0.8*(FamilyDed3(famid(i,1),j,1)-750);
                     FamilyDed3(famid(i,1),j,1)= 750;
                 end
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
            if (income1(i,1) ==4 | income1(i,1) ==5)  
                 if FamilyHosp1(famid(i,1),j,2) > 9000
                    FamilyOOP1(famid(i,1),j,2) = FamilyOOP1(famid(i,1),j,2) + OutofPocket1(i,j,2) - (FamilyHosp1(famid(i,1),j,2)-9000);
                    FamilyHosp1(famid(i,1),j,2)=9000;
                 end
                 
                 if FamilyHosp1(famid(i,1),j,2) <= 9000 && FamilyDed1(famid(i,1),j,2) < 1500
                     FamilyOOP1(famid(i,1),j,2) = FamilyOOP1(famid(i,1),j,2) + OutofPocket1(i,j,2);
                 end
       
                 if FamilyHosp1(famid(i,1),j,2) < 9000 && FamilyDed1(famid(i,1),j,2) >= 1500
                     FamilyOOP1(famid(i,1),j,2) = FamilyOOP1(famid(i,1),j,2) + OutofPocket1(i,j,2) - 0.8*(FamilyDed1(famid(i,1),j,2)-1500);
                     FamilyDed1(famid(i,1),j,2)= 1500;
                 end
            end
            if (income2(i,1) ==4 | income2(i,1) ==5)
                 if FamilyHosp2(famid(i,1),j,2) > 9000
                    FamilyOOP2(famid(i,1),j,2) = FamilyOOP2(famid(i,1),j,2) + OutofPocket2(i,j,2) - (FamilyHosp2(famid(i,1),j,2)-9000);
                    FamilyHosp2(famid(i,1),j,2)=9000;
                 end
            
                 if FamilyHosp2(famid(i,1),j,2) <= 9000 && FamilyDed2(famid(i,1),j,2) < 1500
                     FamilyOOP2(famid(i,1),j,2) = FamilyOOP2(famid(i,1),j,2) + OutofPocket2(i,j,2);
                 end
       
                 if FamilyHosp2(famid(i,1),j,2) < 9000 && FamilyDed2(famid(i,1),j,2) >= 1500
                     FamilyOOP2(famid(i,1),j,2) = FamilyOOP2(famid(i,1),j,2) + OutofPocket2(i,j,2) - 0.8*(FamilyDed2(famid(i,1),j,2)-1500);
                     FamilyDed2(famid(i,1),j,2)= 1500;
                 end

                 
                 if FamilyHosp3(famid(i,1),j,2) > 9000
                    FamilyOOP3(famid(i,1),j,2) = FamilyOOP3(famid(i,1),j,2) + OutofPocket3(i,j,2) - (FamilyHosp3(famid(i,1),j,2)-9000);
                    FamilyHosp3(famid(i,1),j,2)=9000;
                 end
                 
                 if FamilyHosp3(famid(i,1),j,2) <= 9000 && FamilyDed3(famid(i,1),j,2) < 1500
                     FamilyOOP3(famid(i,1),j,2) = FamilyOOP3(famid(i,1),j,2) + OutofPocket3(i,j,2);
                 end
       
                 if FamilyHosp3(famid(i,1),j,2) < 9000 && FamilyDed3(famid(i,1),j,2) >= 1500
                     FamilyOOP3(famid(i,1),j,2) = FamilyOOP3(famid(i,1),j,2) + OutofPocket3(i,j,2) - 0.8*(FamilyDed3(famid(i,1),j,2)-1500);
                     FamilyDed3(famid(i,1),j,2)= 1500;
                 end
            end

                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       
               if (income1(i,1) ==4 | income1(i,1) ==5)      
                               
                 if FamilyHosp1(famid(i,1),j,3) >= 10000         
                    FamilyOOP1(famid(i,1),j,3) = FamilyOOP1(famid(i,1),j,3) + OutofPocket1(i,j,3)- (FamilyHosp1(famid(i,1),j,3)-10000);
                    FamilyHosp1(famid(i,1),j,3)= 10000;
                 end
                 
                 if FamilyHosp1(famid(i,1),j,3) < 10000 && FamilyDed1(famid(i,1),j,3) < 2400
                     FamilyOOP1(famid(i,1),j,3) = FamilyOOP1(famid(i,1),j,3) + OutofPocket1(i,j,3);
                 end
       
                 if FamilyHosp1(famid(i,1),j,3) < 10000 && FamilyDed1(famid(i,1),j,3) >= 2400
                     FamilyOOP1(famid(i,1),j,3) = FamilyOOP1(famid(i,1),j,3) + OutofPocket1(i,j,3) - 0.8*(FamilyDed1(famid(i,1),j,3)-2400);
                     FamilyDed1(famid(i,1),j,3)=2400;
                 end
               end
               if (income2(i,1) ==4 | income2(i,1) ==5)
                 
                 if FamilyHosp2(famid(i,1),j,3) >= 10000         
                    FamilyOOP2(famid(i,1),j,3) = FamilyOOP2(famid(i,1),j,3) + OutofPocket2(i,j,3)- (FamilyHosp2(famid(i,1),j,3)-10000);
                    FamilyHosp2(famid(i,1),j,3)= 10000;
                 end
                 
                 if FamilyHosp2(famid(i,1),j,3) < 10000 && FamilyDed2(famid(i,1),j,3) < 2400
                     FamilyOOP2(famid(i,1),j,3) = FamilyOOP2(famid(i,1),j,3) + OutofPocket2(i,j,3);
                 end
       
                 if FamilyHosp2(famid(i,1),j,3) < 10000 && FamilyDed2(famid(i,1),j,3) >= 2400
                     FamilyOOP2(famid(i,1),j,3) = FamilyOOP2(famid(i,1),j,3) + OutofPocket2(i,j,3) - 0.8*(FamilyDed2(famid(i,1),j,3)-2400);
                     FamilyDed2(famid(i,1),j,3)=2400;
                 end
                 
                 if FamilyHosp3(famid(i,1),j,3) >= 10000         
                    FamilyOOP3(famid(i,1),j,3) = FamilyOOP3(famid(i,1),j,3) + OutofPocket3(i,j,3)- (FamilyHosp3(famid(i,1),j,3)-10000);
                    FamilyHosp3(famid(i,1),j,3)= 10000;
                 end
                 
                 if FamilyHosp3(famid(i,1),j,3) < 10000 && FamilyDed3(famid(i,1),j,3) < 2800
                     FamilyOOP3(famid(i,1),j,3) = FamilyOOP3(famid(i,1),j,3) + OutofPocket3(i,j,3);
                 end
       
                 if FamilyHosp3(famid(i,1),j,3) < 10000 && FamilyDed3(famid(i,1),j,3) >= 2800
                     FamilyOOP3(famid(i,1),j,3) = FamilyOOP3(famid(i,1),j,3) + OutofPocket3(i,j,3) - 0.8*(FamilyDed3(famid(i,1),j,3)-2800);
                     FamilyDed3(famid(i,1),j,3)=2400;
                 end
            end            
        
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% The remainder of this file is pure data processing. Because variables are originally in individual 
%%%%%%%%%%%%%% form, and we want to aggregate everything into a family decision unit, the goal is to have a set of variables
%%%%%%%%%%%%%% describing each family to be used in the choice model. Thus, this part transforms 
%%%%%%%%%%%%%% all variables into family level variables, and removes unnecessary individual variables. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PPO250OOP1 = zeros(1,Sim);
PPO250OOP2 = zeros(1,Sim);
PPO250OOP3 = zeros(1,Sim);

PPO500OOP1 = zeros(1,Sim);
PPO500OOP2 = zeros(1,Sim);
PPO500OOP3 = zeros(1,Sim);

PPO1200OOP1 = zeros(1,Sim);
PPO1200OOP2 = zeros(1,Sim);
PPO1200OOP3 = zeros(1,Sim);

Inc1 = zeros(1,1);
Inc2 = zeros(1,1);

Tier0 = zeros(1,1);
Tier1 = zeros(1,1);
Tier2 = zeros(1,1);

FamidX = zeros(1,1);
Famsize = zeros(1,1);

J =1;
for k = 1:K
    if IND(k,1)>0
       PPO250OOP1(J,:) = reshape(FamilyOOP1(k,:,1),1,Sim);
       PPO250OOP2(J,:) = reshape(FamilyOOP2(k,:,1),1,Sim);
       PPO250OOP3(J,:) = reshape(FamilyOOP3(k,:,1),1,Sim);
       
       PPO500OOP1(J,:) = reshape(FamilyOOP1(k,:,2),1,Sim);
       PPO500OOP2(J,:) = reshape(FamilyOOP2(k,:,2),1,Sim);
       PPO500OOP3(J,:) = reshape(FamilyOOP3(k,:,2),1,Sim);
    
       PPO1200OOP1(J,:) = reshape(FamilyOOP1(k,:,3),1,Sim);
       PPO1200OOP2(J,:) = reshape(FamilyOOP2(k,:,3),1,Sim);
       PPO1200OOP3(J,:) = reshape(FamilyOOP3(k,:,3),1,Sim);      
       
       FamidX(J,1) = k;
       Famsize(J,1) = IND(k,1);
       J = J + 1;
    end
end

Inc1 = zeros(J-1,1);
Inc2 = zeros(J-1,1);

Tier0 = zeros(J-1,1);
Tier1 = zeros(J-1,1);
Tier2 = zeros(J-1,1);
 
FSAY1 = zeros(J-1,1);
FSAY2 = zeros(J-1,1);
FSAY3 = zeros(J-1,1);

managerX = zeros(J-1,1);
quantsoph = zeros(J-1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Inc1(1,1) = income1(1,1);
Inc2(1,1) = income2(1,1);
       
Tier0(1,1) = HTier0(1,1);
Tier1(1,1) = HTier1(1,1);
Tier2(1,1) = HTier2(1,1);
       
FSAY1(1,1) = FSA1(1,1);
FSAY2(1,1) = FSA2(1,1);
FSAY3(1,1) = FSA3(1,1);

managerX(1,1) = manager(1,1);
quantsoph(1,1) = QuantSoph(1,1);

J=1;

for i = 2:nIs
     
     if famid(i,1) ~= famid(i-1,1)
         J = J+1;
     end    
       
     Inc1(J,1) = income1(i,1);
     Inc2(J,1) = income2(i,1);
       
     Tier0(J,1) = HTier0(i,1);
     Tier1(J,1) = HTier1(i,1);
     Tier2(J,1) = HTier2(i,1);
	 
     FSAY1(J,1) = FSA1(i,1);
     FSAY2(J,1) = FSA2(i,1);
     FSAY3(J,1) = FSA3(i,1);

     managerX(J,1) = manager(i,1);
     quantsophF(J,1) = QuantSoph(i,1);
end

%%%%%%%%%%%%%% Sample67 etc. are family specific variables so they go here
%%%%%%%%%%%%%% above as they are going to be exactly the same for every
%%%%%%%%%%%%%% family member so they can be filled in in this more simple
%%%%%%%%%%%%%% way. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Make Age and Gender
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Matrices by Family
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ages = zeros(J,max(IND(k,1)));
Genders = zeros(J,max(IND(k,1)));

ChronCtY1 = zeros(J,max(IND(k,1)));
ChronCtY2 = zeros(J,max(IND(k,1)));

Ages(1,1) = age(1,1);
Genders(1,1) = gender(1,1);

FAMY1PPO250 = zeros(J,Sim, max(IND(k,1)));
FAMY2PPO250 = zeros(J,Sim, max(IND(k,1)));
FAMY3PPO250 = zeros(J,Sim, max(IND(k,1)));

FAMY1PPO500 = zeros(J,Sim, max(IND(k,1)));
FAMY2PPO500 = zeros(J,Sim, max(IND(k,1)));
FAMY3PPO500 = zeros(J,Sim, max(IND(k,1)));

FAMY1PPO1200 = zeros(J,Sim, max(IND(k,1)));
FAMY2PPO1200 = zeros(J,Sim, max(IND(k,1)));
FAMY3PPO1200 = zeros(J,Sim, max(IND(k,1)));

ChronCtY1(1,1) = ChronCt1(1,1);
ChronCtY2(1,1) = ChronCt2(1,1);

FAMY1PPO250(1,:,1) = reshape(OutofPocket1(1,:,1),1,Sim);
FAMY2PPO250(1,:,1) = reshape(OutofPocket2(1,:,1),1,Sim);
FAMY3PPO250(1,:,1) = reshape(OutofPocket3(1,:,1),1,Sim);

FAMY1PPO500(1,:,1) = reshape(OutofPocket1(1,:,2),1,Sim);
FAMY2PPO500(1,:,1) = reshape(OutofPocket2(1,:,2),1,Sim);
FAMY3PPO500(1,:,1) = reshape(OutofPocket3(1,:,2),1,Sim);

FAMY1PPO1200(1,:,1) = reshape(OutofPocket1(1,:,3),1,Sim);
FAMY2PPO1200(1,:,1) = reshape(OutofPocket2(1,:,3),1,Sim);
FAMY3PPO1200(1,:,1) = reshape(OutofPocket3(1,:,3),1,Sim);

J = 1;

for i = 2:nIs
     
     if famid(i,1) ~= famid(i-1,1)
         J = J+1;
         Q = 1;
     end    
       
    Ages(J,Q) = age(i,1);
    Genders(J,Q) = gender(i,1);
    
    ChronCtY1(J,Q) = ChronCt1(i,1);
    ChronCtY2(J,Q) = ChronCt2(i,1);
    
    FAMY1PPO1200(J,:,Q) = reshape(OutofPocket1(i,:,3),1,Sim);
    FAMY2PPO1200(J,:,Q) = reshape(OutofPocket2(i,:,3),1,Sim);
    FAMY3PPO1200(J,:,Q) = reshape(OutofPocket3(i,:,3),1,Sim);
    
    FAMY1PPO500(J,:,Q) = reshape(OutofPocket1(i,:,2),1,Sim);
    FAMY2PPO500(J,:,Q) = reshape(OutofPocket2(i,:,2),1,Sim);
    FAMY3PPO500(J,:,Q) = reshape(OutofPocket3(i,:,2),1,Sim);
    
    FAMY1PPO250(J,:,Q) = reshape(OutofPocket1(i,:,1),1,Sim);
    FAMY2PPO250(J,:,Q) = reshape(OutofPocket2(i,:,1),1,Sim);
    FAMY3PPO250(J,:,Q) = reshape(OutofPocket3(i,:,1),1,Sim);
   
    Q = Q+1;
end

Total1 = zeros(J,Sim);
Total2 = zeros(J,Sim);
Total3 = zeros(J,Sim);

Total1(1,:) = TotalY1(1,:);
Total2(1,:) = TotalY2(1,:);
Total3(1,:) = TotalY3(1,:);

J = 1; 

for i = 2:nIs
   
    if famid(i,1) ~= famid(i-1,1)
        J = J+1;
    end
	
    Total1(J,:) = Total1(J,:) + TotalY1(i,:);    
    Total2(J,:) = Total2(J,:) + TotalY2(i,:);
    Total3(J,:) = Total3(J,:) + TotalY3(i,:);
end

H = size(Total1,1);

PlanPaid1 = zeros(H,Sim,3);
PlanPaid2 = zeros(H,Sim,3);
PlanPaid3 = zeros(H,Sim,3);

for i = 1:H
    for p = 1:3
       if p ==1
          PlanPaid1(i,:,1) = Total1(i,:) - PPO250OOP1(i,:); 
          PlanPaid2(i,:,1) = Total2(i,:) - PPO250OOP2(i,:);
          PlanPaid3(i,:,1) = Total3(i,:) - PPO250OOP3(i,:);
       end
       if p ==2
          PlanPaid1(i,:,2) = Total1(i,:) - PPO500OOP1(i,:); 
          PlanPaid2(i,:,2) = Total2(i,:) - PPO500OOP2(i,:);
          PlanPaid3(i,:,2) = Total3(i,:) - PPO500OOP3(i,:);
       end
       if p ==3
          PlanPaid1(i,:,3) = Total1(i,:) - PPO1200OOP1(i,:); 
          PlanPaid2(i,:,3) = Total2(i,:) - PPO1200OOP2(i,:);
          PlanPaid3(i,:,3) = Total3(i,:) - PPO1200OOP3(i,:);
       end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Clear all residual variables that are unnecessary    %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% to choice model estimation. 					     %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
clear j k nIs i
clear income1 income2 famid X HDHPC1 HDHPC2 HDHPC3 Sim SIGMA1 SIGMA2 SIGMA3 RXoopCP1 RXoopCP2 RXoopCP3 OVoopCP1 OVoopCP2 OVoopCP3
clear RXP1 RXP2 RXP3 R1 R2 R3 Pharmacy1 Pharmacy2 Pharmacy3 OutofPocket1 OutofPocket2 OutofPocket3 OfficeVis1 OfficeVis2 OfficeVis3 OVP1 OVP2 OVP3
clear MentalHealth1 MentalHealth2 MentalHealth3 MU J K IND IncL1 HplanY1 HplanY2 HplanY3 Hosp1 Hosp2 Hosp3 HTier0 HTier1 HTier2
clear FSA1 FSA2 FSA3 manager QuantSoph age gender ChronCt1 ChronCt2 Q TotalY1 TotalY2 TotalY3
clear HP1 HP2 HP3 FamilyOOP1 FamilyOOP2 FamilyOOP3 FamilyHosp1 FamilyHosp2 FamilyHosp3 FamilyDed1 FamilyDed2 FamilyDed3
clear PrZOV1 OVb1 OVc1 PrZRX1 RXAgeCoeffb1 RXGenderCoeffb1 RXb1 RXAgeCoeffc1 RXGenderCoeffc1 RXc1 PrZHP1 PrHHP1 Hb1 Hc1 PrZMH1 
clear MHAgeCoeffb1 MHGenderCoeffb1 MHb1 MHAgeCoeffc1 MHGenderCoeffc1 MHc1 CorrRXH1 CorrOVH1 CorrOVRX1 PrZOV2 OVb2 OVc2 PrZRX2 RXAgeCoeffb2 
clear RXGenderCoeffb2 RXb2 RXAgeCoeffc2 RXGenderCoeffc2 RXc2 PrZHP2 PrHHP2 Hb2 Hc2 PrZMH2 MHAgeCoeffb2 MHGenderCoeffb2 MHb2 MHAgeCoeffc2 
clear MHGenderCoeffc2 MHc2 CorrRXH2 CorrOVH2 CorrOVRX2 PrZOV3 OVb3 OVc3 PrZRX3 RXAgeCoeffb3 RXGenderCoeffb3 RXb3 RXAgeCoeffc3 RXGenderCoeffc3 
clear RXc3 PrZHP3 PrHHP3 Hb3 Hc3 PrZMH3 MHAgeCoeffb3 MHGenderCoeffb3 MHb3 MHAgeCoeffc3 MHGenderCoeffc3 MHc3 CorrRXH3 CorrOVH3 CorrOVRX3
clear p ans famsize N quantsoph Individual H 
clear FAMY1PPO250 FAMY1PPO500 FAMY1PPO1200 FAMY2PPO250 FAMY2PPO500 FAMY2PPO1200 FAMY3PPO250 FAMY3PPO500 FAMY3PPO1200

pack      

%%%%%%%%%%%%%%%% Save Simulated Data with Cost Model Simulations and Other Processing for Choice Model Estimation %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Use as Input Into Simulated Choice Code that simulates choices for different health plans with   %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% simulated data. Those simulated choices will be the choices used in choice model estimation 	  %%%%%%%%%%%%%%%%%%%%%%%%

save 'ASIN-ChoiceModelData.mat'
       
       
       
       