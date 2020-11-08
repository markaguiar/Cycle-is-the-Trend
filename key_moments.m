% key_moments.m
% extract key moments given parameters  xparam=[sigmaz,sigmag,rhoz,rhog,Gbar,phi];
% if rhoz, rhog, Gbar, phi are not included in xparam, default values are used (defined
% in environment.m).   Need to maintain the above order, so, for example,
% if xparam has a length  three, this sets rhog, Gbar, and phi to the default values, etc.
% if the global variable "impulse_response" is set to one, this will also plot IR's.
% Moments are M=[
%  Stdev of filtered output*100
%  Stdev of unfiltered growth of output*100
%  Stdev of filtered investment/sd output
%  Stdev of filtered consumption/sd output 
%  Stdev of filtered net exports/gdp/sd output
%  Corr of lagged filtered output with output
%  Corr of lagged growth of output with growth of output
%  Corr of filtered net exports/gdp with output
%  Corr of filtered consumtion with output
%  Corr of filtered investment with output
%  Corr of employment with output]
% sigmaz=standard deviation in percentage terms of z shock
% sigmag=standard deviation in percentage terms of g shock
% rhoz=persistence to shock in z  
% rhog=persistence to shock in g 
% Gbar=long run mean of growth
% phi=capital adjustment cost parameter
% If rhoz and rhog are omitted, default values rhoz0 and rhog0 are used.
% If Gbar and phi are omitted, default values Gbar0 and phi0 are used.

function M=key_moments(xparam);

global alpha beta sigma delta BYbar gamma psi Gbar0 phi0 rhog0 rhoz0 use_uhlig impulse_response;

sigmaz=xparam(1);  	   
sigmag=xparam(2);   		

rhoz=rhoz0;
rhog=rhog0;

if length(xparam)==3;
   rhoz=rhoz0;
   rhog=xparam(3);
end;


if length(xparam)>3;
rhoz=xparam(3);        
rhog=xparam(4);       
end;

if length(xparam)<=4;
   Gbar=Gbar0;
   phi=phi0;
elseif length(xparam)==5;
   Gbar=Gbar0;
   phi=xparam(5);
elseif length(xparam)==6;
   Gbar=xparam(5);
   phi=xparam(6);
end;

if use_uhlig==0;
   calculate_moments;
   mm=moments_out;
   M=[mm(1:2,:)*100;mm(3:5,:)./mm(1,1);mm(6:end,:)];
elseif use_uhlig==1;
	calculate_moments_uhlig;

	M=[mm(13,1)*100;mm(13,4)*100;mm(12,1)/mm(13,1);mm(11,1)/mm(13,1);mm(10,1)/mm(13,1);mm(13,3);mm(13,5);mm(10:12,2);mm(6,2)];
end;

if impulse_response==1;
   IMP_SELECT = [4:10];
   impresp;
end;


