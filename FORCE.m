function [output,param] = FORCE(Tair,RH,SW0,LW0,U0,ze,fD,psiS,varargin)

% Inputs:
%      Tair, RH, ustar: air temperature, RH and friction velocity at reference level (zm)
%      SW0, LW0:  total incoming shortwave and longwave radiation
%      ze: solar zenith angle (rad)
%      fD: fraction of diffuce radiation
%      psiS: soil water potential

%% default parameters
param.LAI = 6;	    			% Leaf Area Index
param.emissivity = 0.97;		% leaf emissivity
param.Topt =  42;				% optimal temp for Vcmax (C)
param.gmin = 0.01;				% min stomata conductance (mol m-2 s-1)
param.gmax = 0.5;				% max stomata conductance (mol m-2 s-1)
param.Vcmax0 = 30;				% Vcmax top of the canopy (TOC)
param.cw0 = 2;			        % parameter of carbon cost of water (TOC)
param.heta = 0.075;				% decrease Vcmax with canopy depth
param.zeta = 0.11;				% decrease Rd with canopy depth
param.keta = 0.0;               % increase Kmax with canopy depth (<1)
param.d = 0.15;					% charateristic leaf size (m)
param.rhoL_VIS = 0.07;			% leaf reflectance (VIS)
param.rhoL_NIR = 0.30;			% leaf reflectance (NIR)
param.tauL_VIS = 0.05;			% leaf trasmittance (VIS)
param.tauL_NIR = 0.25;			% leaf trasmittance (NIR)
param.rhoW_VIS = 0.30;			% wood reflectance (VIS)
param.rhoW_NIR = 0.60;			% wood reflectance (NIR)
param.MLA0 = 60;                % mean leaf angle TOC (Deg)
param.zc = 30;					% canopy height (m)
param.Kmax = 5; 			    % plant conductance per unit of leaf area (mmol m-2 s-1 Mpa-1)
param.p50  = -1.5;				% P50
param.tlp  = -2.8;				% turgor loss point
param.fw  = 0.10;				% fraction of woody elements (Olivas et al., 2013)
param.a  = 0.5;					% shape parameter of loss of conductivity
param.kw = 0.5;					% wind extintion coeff
param.graph = false;			% plot vertical distrivutions (true/false)
param.m = 3;		        	% number sunlit classes
param.CI = 0.9;		        	% clumpig index (structure parameter)

%% user-defined parameters
if isstruct(varargin{1})
    param = varargin{1};
else
    for i=1:(nargin-7)/2
        param.(varargin{2*(i-1)+1}) = varargin{2*i};
    end
end
%% assign parameters
emissivity = param.emissivity;  %leaf emissivity
LAI = param.LAI;
Topt =  param.Topt; %optimal temp for photosysntesis
fw = param.fw;
zeta = param.zeta;
heta = param.heta;
keta = param.keta;
gmin = param.gmin;
gmax = param.gmax;


Vcmax0 = param.Vcmax0; %Vcmax top of the canopy

d = param.d; %charateristic leaf size
rho.L(1) = param.rhoL_VIS;
rho.L(2) = param.rhoL_NIR;
rho.W(1) = param.rhoW_VIS;
rho.W(2) = param.rhoW_NIR;
tau.L(1) = param.tauL_VIS;
tau.L(2) = param.tauL_NIR;

%% photosyntetic and hydraulics perameters
Kmax0 = param.Kmax*1e-3;
p50 = param.p50;
tlp = param.tlp;
cw0 = param.cw0; % carbon cost of water
a = param.a;
%% paramterization

% wind and aerodynamic properties
kw = param.kw; %wind extinction coeff within the canpy


% physical air properties
pair = 101325;
cp = 29.3;   % (J mol-1 K-1)
lambda = 18.015*(2500.8-2.36*Tair+0.00016*Tair.^2-0.00006*Tair.^3); % (J mol-1)

sigma = emissivity*5.6703e-8;
esat = 610.94*exp((17.625*Tair)./(Tair+243.04));
ec = RH/100*esat;

rsw_s = 0.15;   % soil reflectance (shortwave)
rlw_s = 0.05;   % soil reflectance (longwave)
Tsoil = 25;     % surface soil temperature
Rg = (1-rlw_s)*5.6703e-8*(Tsoil+273.15).^4; %soil emission
m = param.m;CI=param.CI;
[k0,J0,kstar0,x0,zL,zW] = TwoStreamCoeff(param.MLA0,ze,LAI,fw,CI,m);

%% VIS & NIR
alfaL = 1 - rho.L - tau.L; 
alfaW = 1 - rho.W; 
rho = (1-fw)*rho.L+fw*rho.W;
tau = (1-fw)*tau.L;

[~,x,Rs,Idn_sw,Iup_sw,gap] = RT2S_SW(LAI,SW0,fD,rho,tau,k0,J0,kstar0,x0);
   
%% check absorbance closure
% R0 = 0.5*SW0*(1-fD)/cos(ze);
% Id = (Idn_sw + Iup_sw).*kstar0.';
% a0 = zeros(100,2*m+2);
% a0(:,1)   =   alfaL(1)*(Id(:,1)) + alfaL(2)*(Id(:,2));
% a0(:,m+2) =   alfaW(1)*(Id(:,1)) + alfaW(2)*(Id(:,2));
% for c=1:m
%     a0(:,c+1)     =   a0(:,1)   + (alfaL(1)+alfaL(2))*zL(:,c)*R0;
%     a0(:,m+1+c+1) =   a0(:,m+2) + (alfaW(1)+alfaW(2))*zW(:,c)*R0;
% end
% 
% 
% P = zeros(100,2*m+2);
% P(:,1) = (1-gap.*CI)*(1-fw);
% P(:,m+2) = (1-gap.*CI)*fw*CI;
% for c=2:m+1
% P(:,c)     = gap.*CI.*(1-fw)/m;
% P(:,m+1+c) = gap.*CI.*fw/m;
% end
% 
% clf
% plot(sum(a0.*P,2),x)
% set(gca,'ydir','reverse')
% hold on
% plot(sum(A,2),x)

%% absorbed radiation
R0 = 0.5*SW0*(1-fD)/cos(ze);
Id = (Idn_sw + Iup_sw).*kstar0.';

a0.L = zeros(5,m+1);
a0.W = zeros(5,m+1);
a1.L = zeros(5,m+1);
a0.L(:,1) =   polyfit(x,log(alfaL(1)*Id(:,1) + alfaL(2)*Id(:,2)),4);
a0.W(:,1) =   polyfit(x,log(alfaW(1)*Id(:,1) + alfaW(2)*Id(:,2)),4);
a1.L(:,1) = polyfit(x,log(4*alfaL(1)*Id(:,1)),4);
for c=1:m
    a0.L(:,c+1) =   polyfit(x,log(alfaL(1)*(Id(:,1) + zL(:,c)*R0) + alfaL(2)*(Id(:,2) + zL(:,c)*R0)),4);
    a0.W(:,c+1) =   polyfit(x,log(alfaW(1)*(Id(:,1) + zW(:,c)*R0) + alfaW(2)*(Id(:,2) + zW(:,c)*R0)),4);
    a1.L(:,c+1) = polyfit(x,log(4*alfaL(1)*(Id(:,1) + zL(:,c)*R0)),4);
end


 PS0 = polyfit(x,log(gap*CI),4); %CI taken into account


%% TIR
% ,'FJacobian',@fjac
options = bvpset('RelTol',1e-3,'AbsTol',1e-3,...
    'BCJacobian',@fjacbc,'Vectorized','on','NMax', 1e+5, 'Stats', 'On');
% tic
solinit = bvpinit([0 LAI],[LW0,Rg]);

temp = bvp5c(@TwoStream,@bcfun,solinit,options);
% CPUTime=toc5
% display(['Solver converged, CPU Time: ' num2str(CPUTime,'%2.2f')])

sol = deval(temp,x);
Idn_lw(:,1) = sol(1,:);
Iup_lw(:,1) = sol(2,:);

%% compute final outputs
output.x = x;
output.gs   = zeros(length(x),m+1);
output.A    = zeros(length(x),m+1);
output.evap = zeros(length(x),m+1);
output.an   = zeros(length(x),m+1);
output.Rd   = zeros(length(x),m+1);
output.psiL   = zeros(length(x),m+1);
output.T    = zeros(length(x),2*m+2);
output.h    = zeros(length(x),2*m+2);
U  = U0*exp(-kw*x);
w  = sqrt(U/d);
An = zeros(2,1);Evap=zeros(2,1);Tleaf=zeros(2,1);
Abs_lw0 = emissivity.*(Iup_lw+Idn_lw).*kstar0.';
 for i=1:length(x)

    input.gHb_forced = 0.135*w(i);
    input.gvb_forced = 0.147*w(i);
    input.gcb_forced = 0.110*w(i);
    input.Vcmax = Vcmax0*(1 - heta*x(i));
    input.Rd    = 0.015*Vcmax0*(1 - zeta*x(i));
    input.Jmax = 1.96*input.Vcmax;
    cw  = cw0;
    Kmax = Kmax0*(1 - keta*x(i));
    % PSI_max = Kmax/a*log(exp(-a*p50)+1);
    PSI_s = Kmax/a*log(exp(a*(psiS-p50))+1);
    PSI_p = Kmax/a*log(exp(a*(tlp-p50))+1);

   for c=1:m+1
            %% leaf energy budget
            input.aRad =  Abs_lw0(i) + exp(polyval(a0.L(:,c),x(i))); 
            input.aPAR = exp(polyval(a1.L(:,c),x(i)));
% if c==4 && (i==1 || i==2)
%     ggs = linspace(gmin,gmax);
%              for ii=1:100
%                  [An(ii),Evap(ii),Tleaf(ii)] = Photo(ggs(ii),input);
%              end
% 
%             Amax = (Evap(1)-Evap(end))/(Evap(1)/An(1)-Evap(end)/An(end));
%             E0 =   Evap(end)/An(end)*Amax-Evap(end);
%             clf
% plot(Evap,An,'-')
% hold on;
% plot(Evap,Amax.*Evap./(E0+Evap),'-')
%              % subplot(2,2,i)
%              % plot(Evap,Amax*E0./(E0+Evap).^2);hold on
%              % plot(Evap,cw*PSI_p./(PSI_s-Evap-PSI_p).^2,'--')
%              pause(.1)
% end
            %% linearization
            [An(1),Evap(1),Tleaf(1)] = Photo(gmin,input);
          if Evap(1)>0
            [An(2),Evap(2),Tleaf(2)] = Photo(gmax,input);
             

            Amax = (Evap(2)-Evap(1))/(Evap(2)/An(2)-Evap(1)/An(1));
            E0 =   (An(2)-An(1))/(An(1)/Evap(1)-An(2)/Evap(2));
            T1 = (Tleaf(1)-Tleaf(2))/(Evap(1)-Evap(2));
            T0 = Tleaf(1)-T1*Evap(1);

            dC_dE0 = cw*PSI_p./(PSI_s-Evap(1)-PSI_p).^2;
            dA_dE0 = Amax.*E0/(E0+Evap(1)).^2;

            if dC_dE0>dA_dE0 || PSI_s-PSI_p<Evap(1)
                Estar = Evap(1);
            else
                D = sqrt(Amax.*E0./cw/PSI_p);
                Estar1 = (D*(PSI_s-PSI_p)-E0)./(D+1);
                if D>1
                    Estar2 = (D*(PSI_s-PSI_p)+E0)./(D-1);
                    Estar = min([Estar1 Estar2 Evap(2)]);
                else
                    Estar = min([Estar1 Evap(2)]);
                end
            end
            % xline(Estar)

        output.T(i,c) = T0+T1*Estar;
        output.an(i,c) = Amax*Estar./(E0+Estar);
        output.evap(i,c) = Estar;
        output.psiL(i,c) = p50+1./a*log(exp(a*(psiS-p50-Estar/Kmax))+exp(-a*Estar./Kmax)-1);
        DT = output.T(i,c)-Tair;
        gHb	= 1.4*input.gHb_forced + 0.05*(abs(DT)/d).^0.25;
        gvb = 1.4*input.gvb_forced + 1.09*(abs(DT)/d).^0.25;
        DL = (610.94*exp((17.625*output.T(i,c))./(output.T(i,c)+243.04))-ec)/pair;
        output.h(i,c)  = input.aRad - 2*sigma*(output.T(i,c)+273.15).^4 - lambda.*Estar;
        % output.h(i,c) =  2*cp.*gHb*DT;
        output.gs(i,c) = gvb/(gvb*DL/Estar-1);
        output.Rd(i,c)  =  input.Rd*exp((output.T(i,c)+273.15-298)./(8.414*(output.T(i,c)+273.15)*298)*46390);


          else %condenasation
        output.T(i,c) = Tleaf(1);

        gcb = 1.4*input.gcb_forced + 0.75*(abs(DT)/d).^0.25;
        gcs = gmax/1.6;
        gc  = gcs.*gcb./(gcs+gcb);
        [Ac,Aj,Rd] = FarqhuarModel(gc,input.Vcmax,input.Jmax,input.Rd,Topt,Tleaf(1),input.aPAR);
        output.an(i,c) = min([Ac Aj]);
        output.evap(i,c) = Evap(1);
        output.psiL(i,c) = psiS;
        DT = output.T(i,c)-Tair;
        gHb	= 1.4*input.gHb_forced + 0.05*(abs(DT)/d).^0.25;
        gvb = 1.4*input.gvb_forced + 1.09*(abs(DT)/d).^0.25;
        DL = (610.94*exp((17.625*output.T(i,c))./(output.T(i,c)+243.04))-ec)/pair;
        output.h(i,c) =  2*cp.*gHb*DT;
        output.gs(i,c) = gmax;
        output.Rd(i,c)  =  input.Rd*exp((output.T(i,c)+273.15-298)./(8.414*(output.T(i,c)+273.15)*298)*46390);
          end

        %% wood energy budget
        input.aRad =  Abs_lw0(i) + exp(polyval(a0.W(:,c),x(i))); 
        output.T(i,m+1+c) = fzero(@(x) WoodEnegyBudget(x),[Tair-15 Tair+15]);
        DT = output.T(i,m+1+c)-Tair;
        gHb	= 1.4*input.gHb_forced + 0.05*(abs(DT)/d).^0.25;
        output.h(i,m+1+c) =  2*cp.*gHb*DT;
    end
end


    P = zeros(length(x),2*m+2);
    P(:,1) = (1-gap*CI)*(1-fw);
    P(:,m+2) = (1-gap*CI)*fw;
    for c=2:m+1
    P(:,c) = gap*CI*(1-fw)/m;
    P(:,m+1+c) = gap*CI*fw/m;
    end

    P = P./sum(P,2); %ensure sum(P,2)=1


%% ecosystem fluxes
output.H  = trapz(x,sum(output.h.*P,2));
output.LE = trapz(x,lambda.*sum(output.evap.*P(:,1:m+1),2));
output.GPP = trapz(x,sum((output.an+output.Rd).*P(:,1:m+1),2));
output.NPP = trapz(x,sum((output.an).*P(:,1:m+1),2));
output.RsIn = SW0;
output.RlOut = Iup_lw(1);
output.RlIn = LW0;
output.RsOut = Iup_sw(1,1)+Iup_sw(1,2);
output.G = (1-rlw_s)*Idn_lw(100) - Rg + (1-rsw_s)*(Rs(100,1)+Rs(100,2)+Idn_sw(100,1)+Idn_sw(100,2));
output.Rnet = output.RsIn-output.RsOut + output.RlIn-output.RlOut;

%% averaged profiles
output.TL = sum(output.T(:,1:m+1).*P(:,1:m+1),2)./sum(P(:,1:m+1),2); %only phytoelements
output.TC = sum(output.T.*P,2)./sum(P,2); %including woody elements
output.gs_avg = sum(output.gs.*P(:,1:m+1),2)./sum(P(:,1:m+1),2);
output.psiL_avg = sum(output.psiL.*P(:,1:m+1),2)./sum(P(:,1:m+1),2);
output.evap_avg = sum(output.evap.*P(:,1:m+1),2);
output.an_avg = sum(output.an.*P(:,1:m+1),2);

% output.evapstd = sqrt(sum(output.evap.^2.*P,2)-output.evapm.^2);
% output.Tstd = sqrt(sum(output.T.^2.*P,2)-output.Tm.^2);

%% other outputs
% output.Tskin = (Iup_lw(1)/sigma).^(1/4)-273.15;
% output.ubar = U0;
% output.aPAR = exp(polyval(WS1,x));
%% check energy balance closure
[output.H+output.LE ...
 output.Rnet-output.G]
%% Surface layer couplong variables
TK = Tair+273.15;
rhoa_mol = pair./(8.414*TK);
wt = output.H./(cp.*rhoa_mol);
wr = output.LE./lambda./rhoa_mol;
wtv  = wt+0.61*TK*wr;
% L = -ustar^3*TK/(0.4*9.81*wtv);
% FH = integral(@(z) 1./sqrt(1-16*z/L)./z,zc-d0,zm-d0);
% output.Ta = Tair - wt/(0.4*ustar).*FH;
% ea = ec - wr*pair/(0.4*ustar).*FH;
% output.qa = 0.622*ea./(pair-0.378*ea); %specific humidity (Kg/Kg)
rhoa = (3.5*pair+1.3*ec)/TK; % air density (g/m3)
output.wq = output.LE./lambda*18.015/rhoa;  % water flux (Kg/kg m/s)
output.wt = wt;                             % temperature flux (K m/s
output.wtv = wtv;                           % virtual temperature flux (K m/s)
%% plottings

if param.graph

    v=subplot(221);
    % plot(output.gs,x,'-','Color',[1 1 1]*0.5);
    % hold on
    plot(output.gs_avg,x,'-k','linewidth',2)
    set(v,'ydir','reverse')
    xlabel('stomata conductance (mol m^-^2 s^-^1)')
    ylabel('canopy depth (LAI)')
    hold on

    h=subplot(222);
    % plot(output.T,x,'-','Color',[1 1 1]*0.5);
    % hold on
    plot(output.T_avg,x,'-r','linewidth',2)
    lin1 = xline(Tair);lin1.LineStyle='--';
    xlabel('leaf temeperature (^oC)')
    set(h,'ydir','reverse')
    hold on


    z=subplot(223);
    % plot(output.an,x,'-','Color',[1 1 1]*0.5);
    % hold on
    plot(output.an_avg,x,'k-','linewidth',2)
    set(z,'ydir','reverse')
    xlabel('Net Phothosyntesis (\mumol m^-^2 s^-^1)')
    ylabel('canopy depth (LAI)')
    title(['GPP: ' num2str(output.GPP,3)])
    hold on

    g=subplot(224);
    % plot(output.evap,x,'-','Color',[1 1 1]*0.5);
    % hold on
    plot(output.evap_avg,x,'b-','linewidth',2)
    xlabel('transpiration (mmol m^{-2} s^{-1})')
    set(g,'ydir','reverse')
    title(['LE: ' num2str(output.LE,3) ' H: ' num2str(output.H,3)])
    hold on

end


%% two-stream equations and leaf energy budget
function dY = TwoStream(x,Y)

    % vertical gradients
    J = spline(x0,J0,x);
    kstar = spline(x0,kstar0,x);              % inverse of the average optical depth for diffuse raduation
    alfa = emissivity.*kstar;                 % abserbed coefficient for diffuse radiation
    beta = 1/2*(1-emissivity).*(J+1).*kstar;  % Backward scattering coefficient for diffuse radiation

    U = U0*exp(-kw*x);                      % wind
    Vcmax = Vcmax0*(1 - heta*x);            % Vcmax
    Jmax = 1.96*Vcmax;                      % Jmax
    Rd = 0.015*Vcmax0*(1 - zeta*x);         % leaf dark respiration 
    Kmax = Kmax0*(1 - keta*x);
    fsun = exp(polyval(PS0,x));      % fraction of sunlit leaves (it account for CI)
    Abs_lw = alfa.*(Y(1,:)+Y(2,:));
    PP = zeros(2*m+2,length(x));
    T = zeros(2*m+2,length(x));
    An=zeros(2,1);Evap=zeros(2,1);Tleaf=zeros(2,1);
    for ii=1:length(x)

        w = sqrt(U(ii)/d);
        input.gHb_forced  = 0.135*w;
        input.gvb_forced  = 0.147*w;
        input.gcb_forced  = 0.110*w;
        input.Vcmax = Vcmax(ii);
        input.Jmax = Jmax(ii);
        input.Rd = Rd(ii);
        cw  = cw0;
        PSI_s = Kmax(ii)/a*log(exp(a*(psiS-p50))+1);
        PSI_p = Kmax(ii)/a*log(exp(a*(tlp-p50))+1);

        for cc=1:m+1
            %% leaf energy budget
            input.aRad =  Abs_lw(ii) + exp(polyval(a0.L(:,cc),x(ii))); 
            input.aPAR = exp(polyval(a1.L(:,cc),x(ii)));

            %% linearization
            [An(1),Evap(1),Tleaf(1)] = Photo(gmin,input);
         
          if Evap(1)>0
            [An(2),Evap(2),Tleaf(2)] = Photo(gmax,input);

            Amax = (Evap(2)-Evap(1))/(Evap(2)/An(2)-Evap(1)/An(1));
            E0 =   (An(2)-An(1))/(An(1)/Evap(1)-An(2)/Evap(2));
            T1 = (Tleaf(1)-Tleaf(2))/(Evap(1)-Evap(2));
            T0 = Tleaf(1)-T1*Evap(1);

                    % if E0<0
                    %     999

            % ggs = linspace(gmin,gmax);
            % An = zeros(size(ggs));
            % Evap = zeros(size(ggs));
            % Tleaf = zeros(size(ggs));
            % for jkj=1:100
            % 
            %         [An(jkj),Evap(jkj),Tleaf(jkj)] = Photo(ggs(jkj),input);
            % end
            % 
            % plot(Evap,An)
            % hold on;
            % plot(Evap,Amax.*Evap./(E0+Evap),'-')
            % pause

                    % end

            dC_dE0 = cw*PSI_p./(PSI_s-Evap(1)-PSI_p).^2;
            dA_dE0 = Amax.*E0/(E0+Evap(1)).^2;


            if dC_dE0>dA_dE0 || PSI_s-PSI_p<Evap(1)
                Estar = Evap(1);
            else
                D = sqrt(Amax.*E0./cw/PSI_p);
                Estar1 = (D*(PSI_s-PSI_p)-E0)./(D+1);
                if D>1 
                    Estar2 = (D*(PSI_s-PSI_p)+E0)./(D-1);
                    Estar = min([Estar1 Estar2 Evap(2)]);
                else
                    Estar = min(Estar1,Evap(2));
                end
            end

            T(cc,ii) = T0+T1*Estar;
            
            else  % condensation case
             T(cc,ii) = Tleaf(1);
            end

            %% wood energy budget
            input.aRad =  Abs_lw(ii) + exp(polyval(a0.W(:,cc),x(ii))); 
            T(m+1+cc,ii) = fzero(@(x) WoodEnegyBudget(x),[Tair-15 Tair+15]);
        
            %% fraction of each class
            if cc==1
            PP(1,ii) = (1-fsun(ii))*(1-fw);
            PP(m+2,ii) = (1-fsun(ii))*fw;
            else
            PP(cc,ii) = fsun(ii)*(1-fw)/m;
            PP(m+1+cc,ii) = fsun(ii)*fw/m;
            end
        
        end

    end

    PP = PP./sum(PP); %ensure sum(P)=1

    RL_emit = sum(sigma*(T+273.15).^4.*PP);
    dY(1,:) =  -(alfa + beta).*Y(1,:) + beta.*Y(2,:) + RL_emit;
    dY(2,:) =   (alfa + beta).*Y(2,:) - beta.*Y(1,:) - RL_emit;


end

%% Photosysntesis
function [An,E,TL,h,Rd] = Photo(gs,input)


y0 = LeafEnegyBudget(Tair,gs);
[Tmin,ymin] = fminbnd(@(x) LeafEnegyBudget(x,gs),Tair-10,Tair);
if y0<0 || ymin<0
        try
            TL = fzero(@(x) LeafEnegyBudget(x,gs),[Tmin-10 Tmin]);
        catch
            TL = fzero(@(x) LeafEnegyBudget(x,gs),[Tmin-14 Tmin]);
        end
elseif y0>0 && ymin>0
        try
        	TL = fzero(@(x) LeafEnegyBudget(x,gs),[Tair Tair+10]);
        catch

        	TL = fzero(@(x) LeafEnegyBudget(x,gs),[Tair Tair+30]);
        end
end



    DT = TL-Tair;
    DL = (610.94*exp((17.625*TL)./(TL+243.04))-ec)/pair;
    % gHb = 1.4*input.gHb_forced;
    % gvb = 1.4*input.gvb_forced;
    % gcb = 1.4*input.gcb_forced;

    gHb	= 1.4*input.gHb_forced + 0.05*(abs(DT)/d).^0.25;
    gvb = 1.4*input.gvb_forced + 1.09*(abs(DT)/d).^0.25;
    gcb = 1.4*input.gcb_forced + 0.75*(abs(DT)/d).^0.25;


    h = 2*cp.*gHb*DT;
   if DL>0 
    gv  = gvb*gs/(gvb+gs);
   else
    gv  = gvb;
   end

    gcs = gs/1.6;
    gc  = gcs.*gcb./(gcs+gcb);

    [Ac,Aj,Rd] = FarqhuarModel(gc,input.Vcmax,input.Jmax,input.Rd,Topt,TL,input.aPAR);
    % An = min([Ac Aj]);
    NDA = (Ac-Aj)./(Ac+Aj);
    weight = 1./(1+exp(-NDA*10));
    An = (1-weight)*Ac+weight*Aj;

    E = gv.*DL;

end

%% leaf energy budget
function y = LeafEnegyBudget(x,gs)
    DT = x-Tair;
    DL = (610.94*exp((17.625*x)./(x+243.04))-ec)/pair;
    % gHb = 1.4*input.gHb_forced;
    % gvb = 1.4*input.gvb_forced;
    gHb	= 1.4*input.gHb_forced + 0.05*(abs(DT)/d).^0.25;
    gvb = 1.4*input.gvb_forced + 1.09*(abs(DT)/d).^0.25;
        if DL>0
            gv = gvb*gs/(gvb+gs);
        else
            gv = gvb; %condensation
        end

    y = input.aRad - 2*sigma*(x+273.15).^4 - 2*cp.*gHb*DT - lambda.*gv*DL;

end

%% wood energy budget
function y = WoodEnegyBudget(x)
    DT = x-Tair;
    DL = (610.94*exp((17.625*x)./(x+243.04))-ec)/pair;
    gHb	= 1.4*input.gHb_forced + 0.05*(abs(DT)/d).^0.25;
    if DL>0

     y = input.aRad - 2*sigma*(x+273.15).^4 - 2*cp.*gHb*DT;
    else
     gvb = 1.4*input.gvb_forced + 1.09*(abs(DT)/d).^0.25;
     y = input.aRad - 2*sigma*(x+273.15).^4 - 2*cp.*gHb*DT - lambda.*gvb*DL;
    end

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
