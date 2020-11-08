% environment.m
% this program sets the default parameters for the model

clear all;

global alpha beta sigma delta BYbar gamma psi Gbar0 phi0 rhoz0 rhog0 use_uhlig impulse_response;

% Setting parameters:
alpha=0.68;     % labor share
beta=1/1.02;    % discount factor  for quarterly   
sigma=2;        % risk aversion
delta=0.1/2;      % depreciation 
BYbar=0.1;      % steady state bond holdings  

gamma=.36;        % consumption coef in Cobb-Douglas utility

psi=0.001;       % risk premium parameter on debt

use_uhlig=0;     %if ==1, use altered version of Uhlig's moments.m program

impulse_response=0; %if ==1, plots impulse responses when using key_moments.m

%default parameters
Gbar0 = 1.006; 
rhoz0 = 0.95;
rhog0 = 0.01;
phi0  = 4;



%column number for moments when using gmm.m
sd_y=1;
sd_dy=2;
sd_i=3;
sd_c=4;
sd_nx=5;
rho_y=6;
rho_dy=7;
rho_nx=8;
rho_c=9;
rho_i=10;
mu_g=11;


