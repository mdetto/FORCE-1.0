function output = FORCE(Tair,RH,SW0,LW0,ustar,ze,fD,psiS,varargin)

% Inputs:
%      Tair, RH, ustar: air temperature, RH and friction velocity at reference level (zm)
%      SW0, LW0:  total incoming shortwave and longwave radiation
%      ze: solar zenith angle (rad)
%      fD: fraction of diffuce radiation
%      psiS: soil water potential 

%% default parameters
param.LAI = 6;					% Leaf Area Index
param.emissivity = 0.97;		% leaf emissivity
param.Topt =  42;				% optimal temp for Vcmax (C)
param.gmin = 0.01;				% min stomata conductance (mol m-2 s-1)
param.gmax = 0.5;				% max stomata conductance (mol m-2 s-1)
param.Vcmax0 = 25;				% Vcmax top of the canopy (TOC)
param.cw0 = 15;			        % parameter of carbon cost of water (TOC)
param.zeta = 0.0;				% increase in cw with canopy depth
param.heta = 0.45;				% decrease Vcmax with canopy depth (defualt = 34/75)
param.keta = 0.0;               % increase Kmax with canopy depth (<1)
param.d = 0.15;					% charateristic leaf size (m)
param.rho_VIS = 0.11;			% leaf reflectance (VIS)
param.rho_NIR = 0.45;			% leaf reflectance (NIR)
param.ME0 = 60;                 % mean leaf angle TOC (Deg)
param.zc = 30;					% canopy height (m)
param.Kmax = 4;					% plant conductance per unit of leaf area (mmol m-2 s-1 Mpa-1)
param.p50  = -1.5;				% P50
param.a  = 1;					% shape parameter of loss of conductivity
param.k = 0.5;					% wind extintion coeff
param.OptType = false;          % only use light-limited
param.graph = false;			% plot vertical distrivutions (true/false)

%% user-defined parameters
for i=1:(nargin-7)/2
    param.(varargin{2*(i-1)+1}) = varargin{2*i};
end
%% assign parameters
emissivity = param.emissivity;  %leaf emissivity
LAI = param.LAI;
Topt =  param.Topt; %optimal temp for photosysntesis

% zeta = param.zeta;
% heta = param.heta;
keta = param.keta;
gmin = param.gmin;
gmax = param.gmax;

Vcmax0 = param.Vcmax0; %Vcmax top of the canopy
cw0 = param.cw0; % carbon cost of water
d = param.d; %charateristic leaf size
rho_VIS = param.rho_VIS;
rho_NIR = param.rho_NIR;
Kmax0 = param.Kmax*1e-3;
p50 = param.p50;
a = param.a;
type = false;
%% paramterization

% wind and aerodynamic properties
% zm = 42;

z0m = 1.5;
% z0H = 0.2*z0m;
d0 = 0.7*param.zc;
U0 = ustar/0.4.*log((param.zc-d0)/z0m);
k = param.k; %wind extinction coeff within the canpy


% physical air properties
p = 101325;
cp = 29.3;   % (J mol-1 K-1)
lambda = 18.015*(2500.8-2.36*Tair+0.00016*Tair.^2-0.00006*Tair.^3); % (J mol-1)

sigma = emissivity*5.6703e-8;
esat = 610.94*exp((17.625*Tair)./(Tair+243.04));
ec = RH/100*esat;
D = (esat-ec)/p; 


rsw_s = 0.15;   %soil reflectance (shortwave)
rlw_s = 0.05;   %soil reflectance (longwave)
Tsoil = 28;
Rg = (1-rlw_s)*sigma*(Tsoil+273.15).^4; %soil emission

omega = 1-emissivity;
delta = 1-emissivity;
alfa = 1 - omega;                      % Absorption coefficient for incoming diffuse radiation

[k0,J0,x0] = Gfunction(param.ME0,ze,LAI);

% %speherical leaf angle distribution
% x0 = linspace(0,LAI);
% J0 = 1/3*ones(100,1);
% k0 = 1/2/cos(ze)*ones(100,1);

% gamma = 1/2*(omega+J*delta);           % Backward scattering coefficient for incoming diffuse radiation

%% short-wave radiation
if SW0>0
[A,x,Rs,Idn_sw,Iup_sw] = RT2S_SW(LAI,SW0,fD,[rho_VIS rho_NIR],k0,J0,x0);
% [A,x,Rs,Idn_sw,Iup_sw] = RT2S_analytical(LAI,SW0,ze,fD,[rho_VIS rho_NIR],'GRoss',1/2,'J',1/3);
A0 = A(:,1)+A(:,2);
PAR = 4*A(:,1);
WS0 = polyfit(x,log(A0.'),4);
WS1 = polyfit(x,log(PAR.'),4);
else
    WS0 = [0 -inf];
    WS1 = [0 -inf];
    x = linspace(0,LAI,100);
end

%% long-wave
% ,'FJacobian',@fjac
options = bvpset('RelTol',1e-3,'AbsTol',1e-3,...
    'BCJacobian',@fjacbc,'Vectorized','on','NMax', 1e+5, 'Stats', 'On');
% tic

solinit = bvpinit([0 LAI],[LW0,Rg]);


temp = bvp5c(@TwoStream,@bcfun,solinit,options);


% CPUTime=toc5
% display(['Solver converged, CPU Time: ' num2str(CPUTime,'%2.2f')])

sol = deval(temp,x);
Idn(:,1) = sol(1,:);
Iup(:,1) = sol(2,:);

%% absorbed radiation
output.x = x;
output.gs = zeros(size(x));
output.T = zeros(size(x));
output.A = zeros(size(x));
output.h = zeros(size(x));
output.evap = zeros(size(x));
output.an = zeros(size(x));
output.Rd = zeros(size(x));
cl = zeros(size(x));
aRad = exp(polyval(WS0,x));
aPAR = exp(polyval(WS1,x));
U = U0*exp(-k*x);
w = sqrt(U/d);
An=zeros(2,1);Evap=zeros(2,1);Tleaf=zeros(2,1);Hsen=zeros(2,1);Resp=zeros(2,1);
for i=1:length(x)

    input.gHb_forced = 0.135*w(i);
    input.gvb_forced = 0.147*w(i);
    input.gcb_forced = 0.110*w(i);


    input.Vcmax = Vcmax0*(1 - 0.07*x(i));
    input.Rd    = 0.015*Vcmax0*(1 - 0.11*x(i));
    input.Jmax = 1.96*input.Vcmax;
    input.cw  = cw0;
    input.Rabs  = alfa.*(Iup(i) + Idn(i)) + aRad(i);
    input.PAR   = aPAR(i);
  
%% linearization between gmin and gmax
[An(1),Evap(1),Tleaf(1),Hsen(1),Resp(1)] = CarbonGain(gmin,input,type);
[An(2),Evap(2),Tleaf(2),Hsen(2),Resp(2)] = CarbonGain(gmax,input,type);

Amax = (Evap(1)-Evap(2))/(Evap(1)/An(1)-Evap(2)/An(2));
E0 =   Evap(2)/An(2)*Amax-Evap(2);
T1 = (Tleaf(1)-Tleaf(2))/(Evap(1)-Evap(2));
T0 = Tleaf(1)-T1*Evap(1);
 	
% if i==1
% 	gx = linspace(gmin,gmax);
% 	Anx = zeros(100,1);Evapx = zeros(100,1);Tleafx = zeros(100,1);
% 	for j=1:100
% 	[Anx(j),Evapx(j),Tleafx(j)] = CarbonGain(gx(j),input,type);
% 	end
% 	
% end

Kmax = Kmax0*(1 - keta*x(i));
dC_dE = MarCost(Evap);
Emax = Kmax/a*log(exp(a*(psiS-p50))+1);


if Emax>Evap(2)

	if Amax*E0/(E0+Evap(1))^2 > dC_dE(1) && Amax*E0/(E0+Evap(2))^2 < dC_dE(2)
		
		Estar = fzero(@(x) x+E0-sqrt(Amax*E0./MarCost(x)),[Evap(1) Evap(2)]);
% 		output.T(i) = fzero3(Estar);
		output.T(i) = T0+T1*Estar;
		output.an(i) = Amax*Estar./(E0+Estar);
		output.evap(i) = Estar;
		DT = output.T(i)-Tair;
		gHb	= 1.4*input.gHb_forced + 0.05*(abs(DT)/d).^0.25;
		gvb = 1.4*input.gvb_forced + 1.09*(abs(DT)/d).^0.25;
		DL = (610.94*exp((17.625*output.T(i))./(output.T(i)+243.04))-ec)/p;
		output.h(i) =  2*cp.*gHb*DT;
		output.gs(i) = gvb/(gvb*DL/Estar-1);
		output.Rd(i)  =  input.Rd*exp((output.T(i)+273.15-298)./(8.414*(output.T(i)+273.15)*298)*46390);
		
		
	elseif Amax*E0/(E0+Evap(1))^2 < dC_dE(1)
		
		output.T(i) = Tleaf(1);
		output.an(i) = Amax*Evap(1)./(E0+Evap(1));
		output.evap(i) = Evap(1);
		output.h(i) = Hsen(1);
		output.gs(i) = gmin;
		output.Rd(i) = Resp(1);
		
	elseif Amax*E0/(E0+Evap(2))^2 > dC_dE(2)
		
		output.T(i) = Tleaf(2);
		output.an(i) = Amax*Evap(2)./(E0+Evap(2));
		output.evap(i) = Evap(2);
		output.h(i) = Hsen(2);
		output.gs(i) = gmax;
		output.Rd(i) = Resp(2);
	end
	
else %% Emax<Evpa(2)
	if Amax*E0/(E0+Evap(1))^2 > dC_dE(1) && Amax*E0/(E0+Emax)^2 < input.cw*a/Kmax/4
		
		Estar = fzero(@(x) x+E0-sqrt(Amax*E0./MarCost(x)),[Evap(1) Emax]);
% 		output.T(i) = fzero3(Estar);
		output.T(i) = T0+T1*Estar;
		output.an(i) = Amax*Estar./(E0+Estar);
		output.evap(i) = Estar;
		DT = output.T(i)-Tair;
		gHb	= 1.4*input.gHb_forced + 0.05*(abs(DT)/d).^0.25;
		gvb = 1.4*input.gvb_forced + 1.09*(abs(DT)/d).^0.25;
		DL = (610.94*exp((17.625*output.T(i))./(output.T(i)+243.04))-ec)/p;
		output.h(i) =  2*cp.*gHb*DT;
		output.gs(i) = gvb/(gvb*DL/Estar-1);
		output.Rd(i)  =  input.Rd*exp((output.T(i)+273.15-298)./(8.414*(output.T(i)+273.15)*298)*46390);
		
		
	elseif Amax*E0/(E0+Evap(1))^2 < dC_dE(1)
		
		output.T(i) = Tleaf(1);
		output.an(i) = Amax*Evap(1)./(E0+Evap(1));
		output.evap(i) = Evap(1);
		output.h(i) = Hsen(1);
		output.gs(i) = gmin;
		output.Rd(i) = Resp(1);
		
	elseif Amax*E0/(E0+Emax)^2 > input.cw*a/Kmax/4
		
		Estar = sqrt(Amax*E0*4*Kmax/input.cw/a)-E0;
		output.an(i) = Amax*Estar./(E0+Emax);
		output.evap(i) = Estar;
%  		output.T(i) = fzero3(Emax);
		output.T(i) = T0+T1*Estar;
		DT = output.T(i)-Tair;
		gHb	= 1.4*input.gHb_forced + 0.05*(abs(DT)/d).^0.25;
		gvb = 1.4*input.gvb_forced + 1.09*(abs(DT)/d).^0.25;
		DL = (610.94*exp((17.625*output.T(i))./(output.T(i)+243.04))-ec)/p;
		output.h(i) =  2*cp.*gHb*DT;
		output.gs(i) = gvb/(gvb*DL/Estar-1);
		output.Rd(i)  =  input.Rd*exp((output.T(i)+273.15-298)./(8.414*(output.T(i)+273.15)*298)*46390);
		
	end
end


end


%% ecosystem fluxes and surface layer bulk-transfer reletionships
output.H  = trapz(x,output.h);
output.LE = trapz(x,lambda.*output.evap);
output.GPP = trapz(x,output.an+output.Rd);
output.cl = min(x(cl==2)); %height of carbon/light limited
output.RsIn = SW0;
output.RlOut = Iup(1);
output.RlIn = LW0;
output.Tskin = (Iup(1)/sigma).^(1/4)-273.15;
if SW0>0
    output.RsOut = Iup_sw(1,1)+Iup_sw(1,2);
    output.G = (1-rlw_s)*Idn(100) - Rg + (1-rsw_s)*(Rs(100,1)+Rs(100,2)+Idn_sw(100,1)+Idn_sw(100,2));
else
    output.RsOut = 0;
    output.G = (1-rlw_s)*Idn(100) - Rg;
end

output.ubar = U0;
%% check energy balance closure
% [output.H+output.LE ...
%     output.RsIn-output.RsOut+output.RlIn-output.RlOut-output.G]
%
% TK = Tair+273.15;
% rhoa_mol = p./(8.414*TK);
% wt = output.H./(cp.*rhoa_mol);
% wr = output.LE./lambda./rhoa_mol;
% wtv  = wt+0.61*TK*wr;
% L = -ustar^3*TK/(0.4*9.81*wtv);
% FH = integral(@(z) 1./sqrt(1-16*z/L)./z,zc-d0,zm-d0);
% output.Ta = Tair - wt/(0.4*ustar).*FH;
% ea = ec - wr*p/(0.4*ustar).*FH;
% output.qa = 0.622*ea./(p-0.378*ea);
% output.wq = wr*rhoa_mol*1e+3;
% output.wt = wt;
% output.wtv = wtv;
%% plottings

if param.graph
    
v=subplot(221);
plot(output.gs,x);
% lin = yline(output.cl);lin.LineStyle='--';
set(v,'ydir','reverse')
xlabel('stomata conductance (mol m^-^2 s^-^1)')
ylabel('canopy depth (LAI)')
hold on

h=subplot(222);
plot(output.T,x);
lin1 = xline(Tair);lin1.LineStyle='--';
% lin2 = xline(Topt);lin2.LineStyle='--';lin2.Color='r';
xlabel('leaf temeperature (^oC)')
set(h,'ydir','reverse')
hold on


z=subplot(223);
plot(output.an,x);
% lin = yline(output.cl);lin.LineStyle='--';
set(z,'ydir','reverse')
xlabel('Net Phothosyntesis (\mumol m^-^2 s^-^1)')
ylabel('canopy depth (LAI)')
hold on

g=subplot(224);
plot(output.evap,x);
xlabel('transpiration (mmol m^{-2} s^{-1})')
set(g,'ydir','reverse')
hold on

end


%% two-stream equations and leaf energy budget
function dY = TwoStream(x,Y)

% vertical gradients  
J = spline(x0,J0,x);
gamma = 1/2*(omega+J*delta);            % Backward scattering coefficient for incoming diffuse radiation
U = U0*exp(-k*x/LAI);                   % wind
aRad = exp(polyval(WS0,x));             % absorbed short-wave radiation
aPAR = exp(polyval(WS1,x));             % absorbed PAR   
Vcmax = Vcmax0*(1 - 0.075*x);           % Vcmax
Jmax = 1.96*Vcmax;           % Jmax
Rd = 0.015*Vcmax0*(1 - 0.11*x);         % leaf dark respiration


T = zeros(size(x));
An=zeros(2,1);Evap=zeros(2,1);Tleaf=zeros(2,1);
for jj=1:length(x)
    w = sqrt(U(jj)/d);
    input.gHb_forced  = 0.135*w;
    input.gvb_forced  = 0.147*w;
    input.gcb_forced  = 0.110*w;
    input.Rabs =  alfa.*(Y(1,jj)+Y(2,jj)) + aRad(jj);
    input.PAR = aPAR(jj);
    input.Vcmax = Vcmax(jj);
    input.Jmax = Jmax(jj);
    input.Rd = Rd(jj);
    input.cw = cw0;

	Kmax = Kmax0*(1 - keta*x(jj));
%% linearization
[An(1),Evap(1),Tleaf(1)] = CarbonGain(gmin,input,type);
[An(2),Evap(2),Tleaf(2)] = CarbonGain(gmax,input,type);


Amax = (Evap(1)-Evap(2))/(Evap(1)/An(1)-Evap(2)/An(2));
E0 =   Evap(2)/An(2)*Amax-Evap(2);
T1 = (Tleaf(1)-Tleaf(2))/(Evap(1)-Evap(2));
T0 = Tleaf(1)-T1*Evap(1);

dC_dE = MarCost(Evap);

Emax = Kmax/a*log(exp(a*(psiS-p50))+1);

if Emax>Evap(2)
if Amax*E0/(E0+Evap(1))^2 > dC_dE(1) && Amax*E0/(E0+Evap(2))^2 < dC_dE(2)
	
    Estar = fzero(@(x) x+E0-sqrt(Amax*E0./MarCost(x)),[Evap(1) Evap(2)]);
%     T(jj) = fzero3(Estar);
	T(jj) = T0+T1*Estar;
elseif Amax*E0/(E0+Evap(1))^2 < dC_dE(1)
	
	T(jj) = Tleaf(1);
elseif Amax*E0/(E0+Evap(2))^2 > dC_dE(2)
	T(jj) = Tleaf(2);
end

else
	
	if Amax*E0/(E0+Evap(1))^2 > dC_dE(1) && Amax*E0/(E0+Emax)^2 < input.cw*a/Kmax/4
	
    Estar = fzero(@(x) x+E0-sqrt(Amax*E0./MarCost(x)),[Evap(1) Emax]);
%     T(jj) = fzero3(Estar);
	T(jj) = T0+T1*Estar;
elseif Amax*E0/(E0+Evap(1))^2 < dC_dE(1)
	
	T(jj) = Tleaf(1);
elseif Amax*E0/(E0+Emax)^2 >input.cw*a/Kmax/4
	Estar = sqrt(Amax*E0*4*Kmax/input.cw/a)-E0;
	T(jj) = T0+T1*Estar;
	end
end
% [~,~,T(jj)] = CarbonGain(0.1,input,type);

end
LOE = sigma*(T+273.15).^4;
dY(1,:) =  -(alfa + gamma).*Y(1,:) + gamma.*Y(2,:) + LOE;
dY(2,:) =   (alfa + gamma).*Y(2,:) - gamma.*Y(1,:) - LOE;


end

%% Carbon Gain (as function of gs)
function [An,E,TL,h,Rd,cl] = CarbonGain(gs,input,type)


[Tmin,ymin] = fminbnd(@(x) LeafEnegyBudget(x,gs),Tair-10,Tair);

if ymin<0
    TL = fzero(@(x) LeafEnegyBudget(x,gs),[Tair-10 Tmin]);
elseif ymin>0
    try
    	TL = fzero(@(x) LeafEnegyBudget(x,gs),[Tmin Tair+10]);
    catch
    	TL = fzero(@(x) LeafEnegyBudget(x,gs),[Tmin Tair+15]);
    end
end

DT = TL-Tair;
gHb	= 1.4*input.gHb_forced + 0.05*(abs(DT)/d).^0.25;
gvb = 1.4*input.gvb_forced + 1.09*(abs(DT)/d).^0.25;
gcb = 1.4*input.gcb_forced + 0.75*(abs(DT)/d).^0.25;

h = 2*cp.*gHb*DT;

gv  = gvb*gs/(gvb+gs);
gcs = gs/1.6;
gc  = gcs.*gcb./(gcs+gcb);

[Ac,Aj,Rd] = FarqhuarModel(gc,input.Vcmax,input.Jmax,input.Rd,Topt,TL,input.PAR);
if type
	An = Aj; % always light limited
else

	[An,cl] = min([Ac Aj]);
end

DL = (610.94*exp((17.625*TL)./(TL+243.04))-ec)/p;
E = gv.*DL;

end

%% dC/dE
function y = MarCost(x)
	
S = (exp(a*(psiS-p50))+1).*exp(-a*x/Kmax);
y = input.cw*a/Kmax*S./(S+1).^2;
end

%%
function y = LeafEnegyBudget(x,gs)
DT = x-Tair;
gHb	= 1.4*input.gHb_forced + 0.05*(abs(DT)/d).^0.25;
gvb = 1.4*input.gvb_forced + 1.09*(abs(DT)/d).^0.25;
gv = gvb*gs/(gvb+gs);

DL = (610.94*exp((17.625*x)./(x+243.04))-ec)/p;
y = input.Rabs - 2*sigma*(x+273.15).^4 - 2*cp.*gHb*DT - lambda.*gv*DL;

% if y==-inf || y==inf || imag(y)~=0 || isnan(y)
% 	999
% end
end



%% boundary conditions
function res = bcfun(ya,yb)

    res = [ya(1)-LW0;yb(2) - rlw_s*yb(1)-Rg];
end   

%% Jacobians
function [dbcdya,dbcdyb] = fjacbc(~,~)
    
dbcdya=[1 0 
        0 0];

dbcdyb=[0   0
        -rlw_s 1];

end


end