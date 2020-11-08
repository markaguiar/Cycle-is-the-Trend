% This program calculates moments of the model given parameters
% It calls "solve_uhlig.m" to solve the linearized model and then
% calculates moments as described in program_notes.pdf
% moments are collected in the vector "moments_out" and 
% consist of (sigma(y), sigma(dy), sigma(inv), sigma(c), sigma(nx)...
% corr(y,y'), corr(dy,dy'), corr(nx,y), corr(c,y), corr(inv,y), corr(n,y))
% where everything besides dy is HP-filtered.

solve_uhlig; %this program needs to be in the current directory or you need to add the correct path


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Set Notation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%x_t=(khat_t+1,bhat_t+1)'
%zed=(g z)'
%x_t=PP0 x_t-1 + QQ0 zed_t
%s=state variables = (khat,bhat,g,z)'
% s_t+1 = M s_t + epsilon;
% epsilon = (0,0,epsilon_g, epsilon_z)';
% Let Sigma_tilde = E( epsilon epsilon')

M=[PP_0,QQ_0;zeros(2,2),NN];

%F =(chat,n,Ihat,nx,yhat)'
row_c_F=1;
row_n_F=2;
row_I_F=3;
row_nx_F=4;
row_y_F=5;

%F = Hs_t = H1 x(t-1) + H2 zed(t)
% where H1 and H2 are derived from RR and SS 
% row location of consumption, labor, investment, output, and the trade balance (net exports)
% in RR and SS (and in y)
ic=1;
in=3;
idelk=4;
iy=5;
itby=7;



H1 = [RR(ic,1:2);RR(in,1:2);RR(idelk,1:2);RR(itby,1:2);RR(iy,1:2)];

H2 = [SS(ic,:);SS(in,:);SS(idelk,:);SS(itby,:);SS(iy,:)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Get covariance of s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gamm0 = E(ss')
% vec(gamm0) = (I-kron(M,M))^{-1} vec(Sigma_tilde);

Sigma_tilde = [zeros(2,4);zeros(2,2),Sigma];
if psi==0;
   warning off;
end;
gamm0 = (eye(16)-kron(M,M))^(-1)*reshape(Sigma_tilde,16,1);
if psi==0;
   warning on;
end;
gamm0 = reshape(gamm0, 4,4);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Get covariance of first differences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%sbar = (x_t-x_t-1 + g_t-1, w_t, w_t-1)';
%sbar_t = Mbar sbar_t-1 + errorbar
%sbar_t = B0 s_t + B1 s_t-1

%theta1 extracts first row of a vector
theta1=[ones(2,1),zeros(2,1)];

Mbar = [PP_0, QQ_0+theta1, -(QQ_0+PP_0*theta1); zeros(2),NN,zeros(2);zeros(2),eye(2),zeros(2)];

B0=[eye(4);zeros(2,4)];
B1 =[-eye(2),theta1;zeros(2,4);zeros(2),eye(2)];

gammbar0 = B0*gamm0*B0' + B1*gamm0*B1' + B0*M*gamm0*B1'+B1*gamm0'*M'*B0';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Get covariance of first differences of F
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thetaf=[[1;0;1;0;1],zeros(5,1)];

%dFbar = Hbar sbar 

Hbar = [H1, H2, (thetaf-H2-H1*theta1)];

Vardf=[Hbar*gammbar0*Hbar',Hbar*Mbar*gammbar0*Hbar'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%HP Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=500;
NN_n=2*N;

jp=(1:1:NN_n)';
jm=(NN_n:-1:1)';
hpap = -(0.894.^jp).*(0.0561*cos(jp*0.112)+0.0558*sin(jp*0.112));
hpam = -(0.894.^jm).*(0.0561*cos(jm*0.112)+0.0558*sin(jm*0.112));
b=[hpam;1-0.0561;hpap];
btilde=cumsum(b);
btilde=btilde(NN_n-N+1:2*N+1+NN_n-N); %makes btilde run from j=-N to j=N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Get covariance of filtered F
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncorr=1; % ncorr is the number of autocorrelations we will calculate for F
M_n=100;

nf=size(Hbar,1); % nf is the size of F
ns=size(gammbar0);% ns is the size of sbar 


gammf=kron(zeros(1,ncorr+1),zeros(ns));
facv=zeros(nf,(ncorr+1)*nf);


for i=0:ncorr;
   for k=0:M_n;
      if k==0;
         gammj=gammbar0;
         gammf(:,i*ns+1:(i+1)*ns)=gammf(:,i*ns+1:(i+1)*ns)+...
            gammj*(btilde(1:2*N+1-i)'*btilde(1+i:2*N+1));
      else;
         gammj=Mbar*gammj;
         if k<=i;
            gammf(:,i*ns+1:(i+1)*ns)=gammf(:,i*ns+1:(i+1)*ns)+...
               gammj*(btilde(1:2*N+1+k-i)'*btilde(1+i-k:2*N+1))+...
               gammj'*(btilde(1:2*N+1-k-i)'*btilde(1+i+k:2*N+1));
         else;
            gammf(:,i*ns+1:(i+1)*ns)=gammf(:,i*ns+1:(i+1)*ns)+...
               gammj*(btilde(1:2*N+1-k+i)'*btilde(1-i+k:2*N+1))+...
               gammj'*(btilde(1:2*N+1-k-i)'*btilde(1+i+k:2*N+1));
         end; %end if k<=i
      end; %end if k==0
   end; %end k
   
   facv(:,i*nf+1:(i+1)*nf)=Hbar*gammf(:,i*ns+1:(i+1)*ns)*Hbar';
   
end;% end i


stdev=sqrt(diag(facv(1:nf,1:nf)));
tcorr=kron(stdev,stdev');
facr=facv./kron(ones(1,ncorr+1),tcorr);
facr(find(abs(facv)<1e-16))=0;%%%%%%%%%%%%%

stdevdF=sqrt(diag(Vardf(1:nf,1:nf)));
tcorrdF=kron(stdevdF,stdevdF');
facrdF=Vardf./kron(ones(1,ncorr+1),tcorrdF);




moments_out=[stdev(row_y_F);stdevdF(row_y_F);stdev(row_I_F);stdev(row_c_F);stdev(row_nx_F);...
      facr(row_y_F,2*row_y_F);facrdF(row_y_F,2*row_y_F);...
      facr(row_y_F,row_nx_F);facr(row_y_F,row_c_F);facr(row_y_F,row_I_F);facr(row_y_F,row_n_F)];

   


