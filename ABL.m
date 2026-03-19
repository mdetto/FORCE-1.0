function [t,h,th,Tair,Tcan,vpd,wt,wq,RH,GPP,ea,ra] = ABL(forcing,ANN,varargin)

% Inputs:
%      Tair, RH, ustar: air temperature, RH and friction velocity at reference level (zm)
%      SW0, LW0:  total incoming shortwave and longwave radiation
%      ze: solar zenith angle (rad)
%      fD: fraction of diffuce radiation
%      psiS: soil water potential


%% default parameters
param.gT = 4.8e-3;				% temperature lapse rate (C m-1)
param.gq = -3.0e-6;		        % humidity lapse rate (Kg Kg-1 m-1)
param.ws =  0;	    			% subsidence (m/s)
param.beta = 0.2;				% enntrainment ratio
param.p = 101325;				% atmospheric pressure (Pa)
param.Hc  = 30;                 % canopy heigth (m)
param.h0 = 100;		    		% initial ABL height
param.th0 = 24;		    		% initial surface temperature
param.RH0 = 70;		    		% initial RH
param.t0 = 7.0;		    		% initial time
param.t1 = 17;		    		% final time


%% user-defined parameters
if nargin>2
if isstruct(varargin{1})
    param = varargin{1};
else
    for i=1:(nargin-2)/2
        param.(varargin{2*(i-1)+1}) = varargin{2*i};
    end
end
end

% ABL Simulate ABL growth
% Detailed explanation goes here

%% parameters
net    =  ANN.net;
me.wt  = ANN.me.wt;
me.wq  = ANN.me.wq;
me.tc  = ANN.me.tc;
me.GPP = ANN.me.GPP;
sd.wt  = ANN.sd.wt;
sd.wq  = ANN.sd.wq;
sd.GPP = ANN.sd.GPP;
sd.tc  = ANN.sd.tc;

ws   = param.ws;
gT   = param.gT;
gq   = param.gq; 
beta = param.beta; 
p    = param.p;
% z0H = param.z0H;
Hc = param.Hc;
d0 = 0.65*Hc;        % zero-plane displacement (m)
%% initial conditions
opts1 = optimset('TolX',1e-2,'display','off');
% opts2 = optimoptions('fsolve','Display','off');

h0   = param.h0;
th0  = param.th0;
RH0 = param.RH0;
esat = 610.94*exp((17.625*th0)./(th0+243.04));
ea0 = RH0/100*esat;
q0 = 0.622*ea0./(p-0.378*ea0);
% xhat = fsolve(@(x) ...
%     InitialCondition(x),[th0 q0 7],options);
% 
% t0 = xhat(3);
t0 = param.t0;
t1 = param.t1;

% % check initial condition
% SW = ppval(forcing.SW0,t0*3600);
% LW = ppval(forcing.LW0,t0*3600);
% ze = ppval(forcing.ze,t0*3600)/180*pi;
% U = ppval(forcing.ubar,t0*3600);
% ustar = ppval(forcing.ustar,t0*3600);
% fD = ppval(forcing.fD,t0*3600);
% 
% ra = log((0.1*h0)/z0H)/(0.4*ustar);
% 
% xhat = fsolve(@(x,x1,x2,x3,x4,x5,x6,x7,x8) ...
%     SurfaceCoupling(x,th0,q0,SW,LW,U,ze,fD,ra),[th0 q0],options);
% 
% 
% Ta = xhat(1);
% qa = xhat(2);
% 
% ea = p*qa/(0.622+0.378*qa);
% esat = 610.94*exp((17.625*Ta)./(Ta+243.04));
% RH = ea/esat*100;
% 
% ymod = net([Ta,RH,SW,LW,U,ze,fD].');
% 
% wt = ymod(1)*sd.wt + me.wt;
% wq = ymod(2)*sd.wq + me.wq;
% Tcan = ymod(3,:).'*sd.tc + me.tc;

%% solver


t = linspace(t0,t1,50)*3600; %sec
dt = t(2)-t(1);

 y0 = [h0+50 th0 q0];
[t,y] = ode45(@ABL_growth,t,y0);

h   = y(:,1);
th  = y(:,2);
q   = y(:,3);

SW = forcing.SW0(t);
LW = forcing.LW0(t);
ze = forcing.ze(t)/180*pi;
U = forcing.ubar(t);
ustar = forcing.ustar(t);
fD = forcing.fD(t);


Tair = zeros(length(t),1);
qair = zeros(length(t),1);
zSL = Hc+0.1*h;
% zSL = min(100,Hc+0.1*h);

for i=1:length(t)

xhat = fsolve(@(x,x1,x2,x3,x4,x5,x6,x7,x8,x9) ...
    SurfaceCoupling(x,zSL(i),th(i),q(i),SW(i),LW(i),U(i),ze(i),fD(i),ustar(i)),[th(i) q(i)],opts1);

Tair(i) = xhat(1);
qair(i) = xhat(2);

end

ea = p*qair./(0.622+0.378*qair);
esat = 610.94*exp((17.625*Tair)./(Tair+243.04));
RH = ea./esat*100;

ymod = net([Tair,RH,SW,LW,U,ze,fD].');

wt   = ymod(1,:).'*sd.wt  + me.wt;
wq   = ymod(2,:).'*sd.wq  + me.wq;
GPP  = ymod(3,:).'*sd.GPP + me.GPP;
Tcan = ymod(4,:).'*sd.tc  + me.tc;

TK = Tair+273.15;
wtv  = wt+0.61*TK.*wq;

% check resistance
Tv = TK.*(1+0.61*qair);
wtv  = wt+0.61*Tv.*wq;
L = -ustar.^3.*Tv./(0.4*9.81.*wtv);

ra = zeros(length(t),1);
for i=1:length(t)
if L(i)<0
    x1 = 1+(1-16*(zSL(i) -d0)/L(i)).^0.5;
    x2 = 1+(1-16*(Hc-d0)/L(i)).^0.5;
    ra(i) = (log((zSL(i)-d0)/(Hc-d0)) - 2*log(x1/x2))/(0.4*ustar(i));
else
    ra(i) = log((zSL(i)-d0)/(Hc-d0) + 5*(zSL(i)-Hc)/L(i))/(0.4*ustar(i));
end
end
% 
% subplot(121)
% plot(Tair-th,wt.*ra,'.');
% subplot(122)
% plot(qair-q,wq.*ra,'.');

vpd = esat-ea;
t = t/3600; % back to hours

%% ABL growth
function dy = ABL_growth(t,y)
t/3600
h   = y(1);
th  = y(2);
q   = y(3);

dth = th0+gT*(h-h0)-th;
dq =  q0 +gq*(h-h0)-q;

SW = forcing.SW0(t);
LW = forcing.LW0(t);
ze = forcing.ze(t)/180*pi;
U = forcing.ubar(t);
ustar = forcing.ustar(t);
fD = forcing.fD(t);

zSL = Hc+0.1*h;
% zSL = min(100,Hc+0.1*h);
xhat = fsolve(@(x,x1,x2,x3,x4,x5,x6,x7,x8,x9) ...
    SurfaceCoupling(x,zSL,th,q,SW,LW,U,ze,fD,ustar),[th q],opts1);

Ta = xhat(1);
qa = xhat(2);

ea = p*qa/(0.622+0.378*qa);
esat = 610.94*exp((17.625*Ta)./(Ta+243.04));
RH = ea/esat*100;

ymod = net([Ta,RH,SW,LW,U,ze,fD].');

wt = ymod(1)*sd.wt + me.wt;
wq = ymod(2)*sd.wq + me.wq;

TK = Ta+273.15;
wtv  = wt+0.61*TK*wq;
thv = (th+273.15)*(1+0.61*q);          % virtual temperature (K)
dthv = dth+0.61*(q*dth+thv*dq+dth*dq); % Eq. (5.7) in Arellano et al., 2015 - Atmospheric Boundary Layer

if wtv>0  %positive bouyancy
    we = (beta * wtv + 5*ustar^3*thv/(9.81*h))./dthv;
       % we = (beta * wtv)./dthv;
else
    we = (5*ustar^3*thv/(9.81*h))./dthv;
    % we = 0;
end

dy(1,:) = we + ws;
dy(2,:) = (wt+we*dth)./h;
dy(3,:) = (wq+we*dq )./h;

end


%% ABL - surface couplings
function F = SurfaceCoupling(x,z,th,q,SW,LW,U,ze,fD,ustar)


Ta(1) = x(1);
qa(1) = x(2);

ea = p*qa/(0.622+0.378*qa);
esat = 610.94*exp((17.625*Ta)./(Ta+243.04));
RH = ea/esat*100;

ymod = net([Ta,RH,SW,LW,U,ze,fD].');

wt = ymod(1)*sd.wt + me.wt;
wq = ymod(2)*sd.wq + me.wq;

Tv = (Ta+273.15)*(1+0.61*qa);
wtv  = wt+0.61*Tv.*wq;
L = -ustar.^3.*Tv./(0.4*9.81.*wtv);

if L<0
    x1 = 1+(1-16*(z -d0)/L).^0.5;
    x2 = 1+(1-16*(Hc-d0)/L).^0.5;
    ra = (log((z-d0)/(Hc-d0)) - 2*log(x1/x2))/(0.4*ustar);
else
    ra = log((z-d0)/(Hc-d0) + 5*(z-Hc)/L)/(0.4*ustar);
end
F(1) = abs((Ta-th) - ra*wt)/th;
F(2) = abs((qa-q ) - ra*wq)/q;

% ra = log((zsl-d0)/z0H)/(0.4*ustar);
% ra = log(z/z0H)/(k*ustar);

end

function F = InitialCondition(x)

Ta(1) = x(1);
qa(1) = x(2);
t = x(3)*3600;

ea = p*qa/(0.622+0.378*qa);
esat = 610.94*exp((17.625*Ta)./(Ta+243.04));
RH = ea/esat*100;

SW = forcing.SW0(t);
LW = forcing.LW0(t);
ze = forcing.ze(t)/180*pi;
U = forcing.ubar(t);
ustar = forcing.ustar(t);
fD = forcing.fD(t);

% ra = log((0.1*h0)/z0H)/(0.4*ustar);


ymod = net([Ta,RH,SW,LW,U,ze,fD].');

wt = ymod(1)*sd.wt + me.wt;
wq = ymod(2)*sd.wq + me.wq;
Tcan = ymod(4,:).'*sd.tc + me.tc;

TK = Ta+273.15;
wtv  = wt+0.61*TK*wq;
L = -ustar^3*TK/(0.4*9.81*wtv);
if L<0
    x1 = 1+(1-16*(z-d0)/L).^0.5;
    x2 = 1+(1-16*(Hc-d0)/L).^0.5;
    ra = (log((z-d0)/(Hc-d0)) - 2*log(x1/x2))/(0.4*ustar);
else
    ra = log((z-d0)/(Hc-d0))/(0.4*ustar);
end

F(1) = abs((Ta-th0) - ra*wt)/th0;
F(2) = abs((qa-q0 ) - ra*wq)/q0;
F(3) = abs(Tcan-Ta)/Ta;

% TK = Ta+273.15;
% F(3)  = abs(wt-0.61*TK*wq);

end

end