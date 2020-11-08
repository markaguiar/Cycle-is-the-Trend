%solve_uhlig
%solves model using Harald Uhlig's toolkit



addpath c:\toolkit4_1; %This is the directory where Uhlig's toolkit is located.

m1=gamma*(1-sigma)-1;
m2=(1-gamma)*(1-sigma);



% solving for Steady State

Qbar=beta*(Gbar^m1);
YKbar=((1/Qbar)-(1-delta))/(1-alpha);
CYbar=1+(1-Gbar-delta)*(1/YKbar)-(1-Gbar*Qbar)*BYbar;
Nbar=(alpha*gamma)/(CYbar-gamma*CYbar+alpha*gamma);
Kbar=(((Gbar^alpha)*(Nbar^alpha))/YKbar)^(1/alpha);
Ybar=(Gbar^(alpha))*(Kbar^(1-alpha))*(Nbar^alpha);
Cbar=CYbar*Ybar;
Xbar=(Gbar-1+delta)*Kbar;
Bbar=BYbar*Ybar;
Lbar=1-Nbar;
nxbar=(Ybar-Cbar-Xbar)/Ybar;

%Declaring the matrices.

%Endogenous state variables "x(t)":  k(t+1), b(t+1)
%Endogenous other variables "y(t)": c(t),l(t),n(t),inv(t),gdp(t),q(t),nx(t)
%Exogenous state variables "z(t)": g(t), zed(t)
% Switch to that notation.  Find matrices for format
% 0 = AA x(t) + BB x(t-1) + CC y(t) + DD z(t)
% 0 = E_t [ FF x(t+1) + GG x(t) + HH x(t-1) + JJ y(t+1) + KK y(t) + LL z(t+1) + MM z(t)]
% z(t+1) = NN z(t) + epsilon(t+1) with E_t [ epsilon(t+1) ] = 0,

%Size of States
EnS = 2; %number of endogenous states
EnV = 7; %number of endogenous variables
ExS = 2; %number of exogenous states

%Number of Equations
ND = EnV;%number of deterministic equations
NF = EnS;%number of forward looking equations
NZ = ExS ;%number of exogenous shock equations

%Normalize Matrices
AA=zeros(ND,EnS);
BB=AA;
CC=zeros(ND,EnV);
DD = zeros(ND,ExS);
FF = zeros(NF,EnS);
GG = FF;
HH = FF;
JJ = zeros(NF,EnV);
KK = JJ;
LL = zeros(NF,ExS);
MM=LL;
NN=zeros(NZ,ExS);
Sigma=zeros(ExS,ExS);

%row number of variables
row_k=1;
row_b=2;
row_c=1;
row_l=2;
row_n=3;
row_inv=4;
row_gdp=5;
row_q=6;
row_nx=7;
row_g=1;
row_zed=2;

%Deterministic Equations
% 0 = AA x(t) + BB x(t-1) + CC y(t) + DD z(t)

%Labor-Leisure
if gamma<1;
	CC(1,row_gdp)=1;
	CC(1,row_n)=-1;
	CC(1,row_c)=-1;
   CC(1,row_l)=1;
elseif gamma==1;
   CC(1,row_l)=1;
end;


%Budget Constraint
CC(2,row_gdp)=Ybar;
CC(2,row_c)=-Cbar;
CC(2,row_inv)=-Xbar;
CC(2,row_q)=Gbar*Qbar*Bbar;
BB(2,row_b)=-Bbar;
AA(2,row_b)=Gbar*Qbar*Bbar;
DD(2,row_g)=Gbar*Qbar*Bbar;

%Investment
CC(3,row_inv)=-Xbar;
AA(3,row_k)=Kbar*Gbar;
BB(3,row_k)=-Kbar*(1-delta);
DD(3,row_g)=Kbar*Gbar;

%Technology
CC(4,row_gdp)=-1;
BB(4,row_k)=1-alpha;
CC(4,row_n)=alpha;
DD(4,row_g)=alpha;
DD(4,row_zed)=1;

%Labor=1-Leisure (time budget constraint)
CC(5,row_l)=Lbar;
CC(5,row_n)=Nbar;

%Interest Rate
CC(6,row_q)=-1;
AA(6,row_b)=-psi*Bbar*Qbar;

%Net Exports
CC(7,row_gdp)=1-nxbar;
CC(7,row_inv)=-Xbar/Ybar;
CC(7,row_c)=-Cbar/Ybar;
CC(7,row_nx)=-1;



% EXPECTATIONAL EQUATIONS:
% 0 = E_t [ FF x(t+1) + GG x(t) + HH x(t-1) + JJ y(t+1) + KK y(t) + LL z(t+1) + MM z(t)]

%K' FOC:
JJ(1,row_c)=m1;
JJ(1,row_l)=m2;
JJ(1,row_gdp)=beta*Gbar^m1*(1-alpha)*Ybar/Kbar;
LL(1,row_g)=beta*Gbar^(m1+2)*phi;
FF(1,row_k)=beta*Gbar^(m1+2)*phi;
GG(1,row_k)=-(beta*Gbar^m1*((1-alpha)*Ybar/Kbar+phi*Gbar^2)+phi*Gbar);
KK(1,row_c)=-m1;
KK(1,row_l)=-m2;
MM(1,row_g)=(m1-phi*Gbar);
HH(1,row_k)=phi*Gbar;

%B' FOC:
JJ(2,row_c)=m1;
JJ(2,row_l)=m2;
MM(2,row_g)=m1;
KK(2,row_c)=-m1;
KK(2,row_l)=-m2;
KK(2,row_q)=-1;

% AUTOREGRESSIVE MATRIX FOR z(t)
% z(t+1) = NN z(t) + epsilon(t+1) with E_t [ epsilon(t+1) ] = 0,

NN(row_g,row_g)=rhog;
NN(row_zed,row_zed)=rhoz;

Sigma(row_g,row_g) =sigmag^2;
Sigma(row_zed,row_zed)=sigmaz^2;




[l_equ,m_states] = size(AA);
[l_equ,n_endog ] = size(CC);
[l_equ,k_exog  ] = size(DD);

  
% Starting the calculations:

[l_equ,m_states] = size(AA);
[l_equ,n_endog ] = size(CC);
[l_equ,k_exog  ] = size(DD);



warnings='';
options;
solve;

% Form expanded expressions
% See program_notes.pdf for explanations
% x(t) will now be (k(t),b(t),Gamma(t));
% y(t) will now be c(t),l(t),n(t),inv(t),gdp(t),q(t),nx(t),C(t),Inv(t),GDP(t);
% where C(t)=c(t)+Gamma(t-1), etc.
% z(t) as before

%Size of States
EnS = 3; %number of endogenous states
EnV = 10; %number of endogenous variables

m_states=EnS;
n_endog=EnV;
k_exog=ExS;
l_equ=EnV;


row_C=8;
row_Inv=9;
row_GDP=10;

VARNAMES = ...
['capital              ';...
 'bonds                ';...
 'Gamma                ';...
 'consumption          ';...
 'leisure              ';...
 'labor                ';...
 'investment           ';...
 'output               ';...
 'interest (q)         ';...
 'trade balance/gdp(nx)';...
 'Consumption          ';...
 'Investment           ';...
 'Output               ';...
 'g                    ';...
 'zed                  '];



PP_0=PP;
QQ_0=QQ;
RR_0=RR;
SS_0=SS;
WW_0=WW;

PP = [PP_0,zeros(size(PP_0,1),1); zeros(1,size(PP_0,2)),1];
QQ = [QQ_0;[1,0]];

RR=[RR_0,zeros(size(RR,1),1)];
RR(row_C,:)=[RR_0(row_c,:),1];
RR(row_Inv,:)=[RR_0(row_inv,:),1];
RR(row_GDP,:)=[RR_0(row_gdp,:),1];

SS(row_C,:)=SS_0(row_c,:);
SS(row_Inv,:)=SS_0(row_inv,:);
SS(row_GDP,:)=SS_0(row_gdp,:);


WW = [ eye(m_states)         , zeros(m_states,k_exog);...
       RR*pinv(PP)           , (SS-RR*pinv(PP)*QQ) ;...
       zeros(k_exog,m_states), eye(k_exog)            ];
 
 
HP_SELECT  = 1:(m_states+n_endog+k_exog); % Selecting the variables for the HP Filter calculations.
N_GRIDPOINTS=400;