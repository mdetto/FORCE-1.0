function [Ac,Aj,Rd] = FarqhuarModel(gs,Vcmax25,Rd25,Topt,TL,PAR)
% Net Photosyntesis with stomata optimization coupled with hydraulic transport
%
% Parameters
% Vcmax25:  maximum carboxillation velocity at reference temperature [mumol/m2/s]
% gs: stomata conductance to CO2 [mol m-2 s-1]
% Topt: optimal temperature for Vcamax [Celsius[
%
% environmetal inputs 
% TL:  leaf temperature (C)
% Ca:  ambient CO2 concentration (mubar)
% O:   ambient O2 concentration (mbar)
% PAR: absorbed photosyntetic active radiation (mumol m-2 s-1) 
%
% Outputs:
% Ac, Aj, Rd: Carbon limites, light limited net photosyntesis, dark respiration [mumol m-2 s-1] 
%
% values Medlyn et al. (2002). Temperature response of parameters of a 
% biochemically based model of photosynthesis. II. A review of experimental data. 
% Plant, Cell & Environment, 25, 1167–1179.
% Bernacchi, C. J., Singsaas, E. L., Portis Jr, A. R., Pimentel, C., & Long, S. P. (2001).
% Improved temperature response functions for models of Rubisco-limited photosynthesis. 
% Plant, Cell and Environment, 24(2), 253–259.

%% environmental paramters
Ca  = 400;            % mubar
O  =  210;            % mbar
TK =  TL+273.15;      % Kelvin
% PAR = 0;            % (mumol m-2 s-1) 
%% parametrization
% Vcmax25 = 40;
R = 8.414;            % J mol-1 K-1  
Kc25 = 404.9;         % mubar
E.Kc = 79430;         % J mol-1   
Ko25 = 278.4;         % mbar
E.Ko = 36380;         % J mol-1 
Gstar25 = 42.75;      % mumol/mol
E.Gstar = 37830;      % mubar
% Rd25 = 0.015*Vcmax25; % mumol m-2 s-1
E.Rd = 46390;         % J mol-1 

% Vcmax (from Medelyn et al., 2002)
Ha.Vcmax = 66560;     % J mol-1  
Hd.Vcmax = 200000;    % J mol-1  
% S.Vcmax = 637;        % J mol-1  
S.Vcmax = Hd.Vcmax./(Topt + 273.15)  + R*log(Ha.Vcmax/(Hd.Vcmax-Ha.Vcmax)); % J mol-1 K-1
% Topt = Hd.Vcmax./(S.Vcmax-R*log(Ha.Vcmax/(Hd.Vcmax-Ha.Vcmax)))-273.15;

% Jmax (from Medelyn et al., 2002)
% Jmax25 = exp(1.01+0.89*log(Vcmax25)); % mumol/m2/s (Walker at al., 2014, table 4)
Jmax25 = 1.96*Vcmax25; % mumol/m2/s (Lamour at al., 2023, table 2)

Ha.Jmax = 43965;      % J mol-1  
Hd.Jmax = 200000;     % J mol-1  
% S.Jmax = 638;         % J mol-1 K-1
ToptJ = Topt-5; % usually Topt for Jmax < Topt for Vcmax (Crous et al., 2021, Tansley Review)
S.Jmax = Hd.Jmax./(ToptJ + 273.15)  + R*log(Ha.Jmax/(Hd.Jmax-Ha.Jmax)); % J mol-1 K-1
theta = .7;           % unitless
alpha = 0.36;         % unitless

% Mesophyll conductance
gm25 = 0.05; %(mumol m-2 s-1 Pa-1)
Ha.gm = 45e+3; %J mol-1
%% temperature functions
Kc  =  Kc25*exp((TK-298)./(R*TK*298)*E.Kc);
Ko  =  Ko25*exp((TK-298)./(R*TK*298)*E.Ko);
Gstar  =  Gstar25*exp((TK-298)./(R*TK*298)*E.Gstar);
Rd  =  Rd25*exp((TK-298)./(R*TK*298)*E.Rd);

Vcmax  =  Vcmax25*exp((TK-298)./(R*TK*298)*Ha.Vcmax).*(1+exp((298*S.Vcmax-Hd.Vcmax)/(298*R)))...
       ./(1+exp((S.Vcmax*TK-Hd.Vcmax)./(R*TK)));
   

Jmax   =   Jmax25*exp((TK-298)./(R*TK*298)*Ha.Jmax).*(1+exp((298*S.Jmax-Hd.Jmax)/(298*R)))...
          ./(1+exp((S.Jmax*TK-Hd.Jmax)./(R*TK)));
Km = Kc.*(1+O./Ko); 

%% no temperature dependence
% Kc  =  Kc25;
% Ko  =  Ko25;
% Gstar  =  Gstar25;
% Rd  =  Rd25;
% Vcmax  =  Vcmax25;
% Jmax   =   Jmax25;
% Km = Kc.*(1+O./Ko); 
%%
gm = gm25*exp(Ha.gm./(R*TK));
% gm = gm25;
if gm == inf
    gc = gs;
else
    gc = gs.*gm./(gs+gm);
end
%% carbon limited
a  =  (Vcmax-Rd)/2;
b  =  (Ca+Km)/2;
c  =  Rd./2.*(Ca+Km) + Vcmax./2.*(2*Gstar-Ca+Km);

Ac = a + b.*gc - sqrt(b.^2.*gc.^2+c.*gc+a.^2);


if PAR>0
%% light limited

J  =  (alpha.*PAR + Jmax - sqrt((alpha.*PAR + Jmax).^2 - 4*alpha*Jmax.*theta.*PAR))./(2*theta);
a  =  J/8-Rd/2;
b  =  Ca/2+Gstar;
c  =  Rd./2.*(Ca  + 2*Gstar) + J./2.*(Gstar - Ca/4);

Aj = a + b.*gc - sqrt(b.^2.*gc.^2+c.*gc+a.^2);
else
    Aj = -Rd;
end

%% co-limitation Vico et al., (2013)
% elseif strcmp(dat.model,'colimit')
%     
% J  =  (alpha.*dat.I + Jmax - sqrt((alpha.*dat.I+Jmax).^2 - 4*alpha*Jmax.*theta.*dat.I))./(2*theta);
% a  =  J/8-Rd/2;
% b  =  Ca+Km.*J./Vcmax/4;
% c  =  Rd./2.*(Ca+Km.*J./Vcmax/4) + J./8.*(2*Gstar-Ca + Km.*J./Vcmax/4);
% 
% 
% An  =  (2*a.*b-c)./(2*b) - sqrt(a.^2-(c./(2*b)).^2).*sqrt(LD./(2*b-LD));
% if  real(An)<0
%     gmax = 0;
% else
%     gmax  =  sqrt((2*a.*b-c).*(2*a.*b+c))./(2*b.^2).*(b-LD)./sqrt(2*b.*LD-LD.^2) - c./(2*b.^2);
% end
% 
% g0 = min(gmax,glim);
% if gmax>0
%     
%     A3 = exp(-gamma*D/Kmax);
%     A4 = (lambda0*D).^(gamma*psi0);
%     gs = fzero(@(x) ...
%         ((b - (b^2*x + c/2)/sqrt(b^2*x^2+c*x+a^2))).^(gamma*psi0) - ...
%         (A1*A3.^(1.6*x)-A2).*A4,...
%         [0 g0],opts);
% end
% 
% An = a + b.*gs - sqrt(b.^2.*gs.^2+c.*gs+a.^2);
% gs0 = 1.6*gs;
% psiL = p50 + 1/gamma*log((exp(gamma*(psiS-p50)) + 1)*exp(-gamma*gs0*D/Kmax) - 1);
% end


% check analitycal solutions
% gs  =  sqrt((a.*b-c).*(a.*b+c))./b.^2.*(b-2*LD)./sqrt(b.*LD-LD.^2) - 2*c./b.^2;
% cc(:,1)  =  Ca - Ac(:,1)./gs;
% Ac(:,2)  =  Vcmax./(Km+cc).*(cc-Gstar)-Rd;

% gj  =  sqrt((a.*b-c).*(a.*b+c))./b.^2.*(b-2*LD)./sqrt(b.*LD-LD.^2) - 2*c./b.^2;
% cj  =  Ca - Aj(:,1)./gj;
% Aj(:,2)  =  J./(4*cj+8*Gstar).*(cj-Gstar)-Rd;



%% Medelyn et al, 2002

% dat  =  readtable('C:\Users\mdetto\Dropbox (Smithsonian)\paper\Rubisco\Medelyn et al, 2002.xlsx');
% R = 8.414*1e-3;
% Ha.Jmax    =  median(dat.Ha);
% Hd.Jmax    =  median(dat.Hd);
% Topt.Jmax  =  median(dat.Topt);
% S.Jmax     =  median(R*log(dat.Ha./(dat.Hd-dat.Ha)) + dat.Hd./(273+dat.Topt));
% Ha.Vcmax   =  median(dat.Ha_1);
% Hd.Vcmax   =  median(dat.Hd_1);
% S.Vcmax    =  median(R*log(dat.Ha_1./(dat.Hd_1-dat.Ha_1)) + dat.Hd_1./(273+dat.Topt_1));

