% This program computes the random coefficeints discrete choice model described 
% the "A Research Assistant's Guide to Discrete Choice Models of Demand," NBER technical
% paper #221, and "Measuring Market Power in the Ready-to-Eat Cereal Industry," NBER WP
% #6387. 

% Written by Aviv Nevo, May 1998.
global invA ns x1 x2 s_jt IV vfull dfull theta1 theti thetj cdid cdindex

% load data. see description in readme-ps2 datasets.txt
load ps2
load iv

ns = 20;       % number of simulated "indviduals" per market %
nmkt = 100;     % number of markets = (# of cities)*(# of quarters)  %
nbrn = 5;     % number of brands per market. if the numebr differs by market this requires some "accounting" vector %
ninstr = 4;   % number of price IVs

% To run original files, uncomment below 5 lines
% load ps2_old
% load iv_old
% nmkt = 94;
% nbrn = 24;
% ninstr = 20;

IV = [iv(:,2:ninstr+1) x1(:,2:nbrn+1)];
clear iv

% this vector relates each observation to the market it is in %
cdid = kron([1:nmkt]',ones(nbrn,1));    
% this vector provides for each index the of the last observation %
% in the data used here all brands appear in all markets. if this %
% is not the case the two vectors, cdid and cdindex, have to be   % 
% created in a different fashion but the rest of the program works fine.%
cdindex = [nbrn:nbrn:nbrn*nmkt]';       


% starting values. zero elements in the following matrix correspond to %
% coeff that will not be max over,i.e are fixed at zero. % 
%            sigma    income   income^2    age       child
theta2w=    [0.3302   5.4819         0    0.2037         0;
             2.4526  15.8935    -1.2000        0    2.6342;
             0.0163  -0.2506         0    0.0511         0;
             0.2441   1.2650         0   -0.8091         0];

% create a vector of the non-zero elements in the above matrix, and the %
% corresponding row and column indices. this facilitates passing values % 
% to the functions below. %
[theti, thetj, theta2]=find(theta2w);
          
horz=['    mean       sigma      income   income^2    age    child'];
vert=['constant  ';
      'price     ';
      'sugar     ';
      'mushy     '];

% create weight matrix
invA = inv([IV'*IV]);

% Logit results and save the mean utility as initial values for the search below

% compute the outside good market share by market
temp = cumsum(s_jt);
sum1 = temp(cdindex,:);
sum1(2:size(sum1,1),:) = diff(sum1);
outshr = 1.0 - sum1(cdid,:);

y = log(s_jt) - log(outshr);
mid = x1'*IV*invA*IV';
t = inv(mid*x1)*mid*y;
mvalold = x1*t;
oldt2 = zeros(size(theta2));
mvalold = exp(mvalold);

save mvalold mvalold oldt2
clear mid y outshr t oldt2 mvalold temp sum1

vfull = v(cdid,:);
dfull = demogr(cdid,:);

load foptions-acw.mat
options = foptions;
options(2) = 0.1;
options(3) = 0.01;

tic
% the following line computes the estimates using a Quasi-Newton method % 
% with an *analytic* gradient %
options = optimset('GradObj','on'); % Gradient is supplied
[theta2,FVAL,EXITFLAG,OUTPUT,GRAD,HESSIAN] = fminunc({@gmmobj,@gradobj},theta2,options);
options

% the following line computes the estimates using a simplex search method
%theta2 = fmins('gmmobj',theta2)
comp_t = toc/60;

% computing the s.e.
vcov = var_cov(theta2);
se = sqrt(diag(vcov));

theta2w = full(sparse(theti,thetj,theta2));
t = size(se,1) - size(theta2,1);
se2w = full(sparse(theti,thetj,se(t+1:size(se,1))));

% the MD estimates
omega = inv(vcov(2:25,2:25));
xmd = [x2(1:24,1) x2(1:24,3:4)];
ymd = theta1(2:25);

beta = inv(xmd'*omega*xmd)*xmd'*omega*ymd;
resmd = ymd - xmd*beta;
semd = sqrt(diag(inv(xmd'*omega*xmd)));
mcoef = [beta(1); theta1(1); beta(2:3)];
semcoef = [semd(1); se(1); semd];

Rsq = 1-((resmd-mean(resmd))'*(resmd-mean(resmd)))/((ymd-mean(ymd))'*(ymd-mean(ymd)));
Rsq_G = 1-(resmd'*omega*resmd)/((ymd-mean(ymd))'*omega*(ymd-mean(ymd)));
Chisq = size(id,1)*resmd'*omega*resmd;

diary results
disp(horz)
disp('  ')
for i=1:size(theta2w,1)
     disp(vert(i,:))
     disp([mcoef(i) theta2w(i,:)])
     disp([full(semcoef(i)) se2w(i,:)])
end

disp(['GMM objective:  ' num2str(FVAL)])
disp(['MD R-squared:  ' num2str(Rsq)])
disp(['MD weighted R-squared:  ' num2str(Rsq_G)])
disp(['# of objective function evaluations:  ' num2str(OUTPUT.iterations)])
disp(['run time (minutes):  ' num2str(comp_t)])
diary off
