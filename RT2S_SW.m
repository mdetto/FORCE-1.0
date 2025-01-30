function [A,x,R0,Idn,Iup,gap,albedo] = RT2S_SW(LAI,Rs,f,RHO,k0,J0,x0)
% Radiative transfer model for shortwave radiation
% Description: the model is a numerical solution of a two-stream approximation
% It is run for visible and NIR separately, assuming 50% of incoming radiation is vis
% and 50% is NIR. It assumes homogeneous canopy with spherical leaf angle distribution
% (different distributions with veritcal variation can be implemented)
% The parametrization of the scatter coefficients are taken from Yuan et al., (2017, JAMES - table 2)
%
% Inputs:
%      LAI: leaf area index
%      Rs:  total incoming radiation at the top of the canopy
%      ze:  solar zenit angle (rad)
%      f:   fraction of diffuse radiation

% Outputs:
%      A:  absorbed radiation
%      x:  canopy depth (cumulative LAI)
%      R0: direct radiation
%      Idn/Iup: diffuse radiation (up and down)
%      albedo: the sum of reflected vis and NIR divided by I0_tot
options_1 = odeset('RelTol',1e-6,'AbsTol',1e-6,'Jacobian',@jac1,'Vectorized','on');%,...'InitialSlope',-k0(1));

options_2 = bvpset('RelTol',1e-6,'AbsTol',1e-6,'FJacobian',@fjac,...
    'BCJacobian',@fjacbc,'Vectorized','on','NMax', 1e+5, 'Stats', 'On');

% % spherical leaf angle distribution
% G = 1/2;
% J = 1/3;
%%%% Visible
% rho = 0.10;  % leaf reflectance
% tau = 0.05;  % leaf transmittance
% rg  = 0.15;  % soil reflectance
%%%% NIR
% rho = 0.40;  % leaf reflectance
% tau = 0.25;  % leaf transmittance
% rg  = 0.15;  % soil reflectance

TAU = [0.05 0.25]; % leaf transmittance
RG = [0.15 0.30];  % soil reflectance
s = [0.5 0.5]; %fraction of incoming radiation in band VIS and NIR 
N = 100;
x = linspace(0,LAI,N);
Idn = zeros(N,2);
Iup = zeros(N,2);
R0 = zeros(N,2);
A = zeros(N,2);

% gap fraction
temp1 = ode15s(@Direct,[0 LAI],1,options_1);
gap = deval(temp1,x).';
for i=1:2

rho = RHO(i);  % leaf reflectance
tau = TAU(i);  % leaf transmittance
rg  = RG(i);  % soil reflectance
omega = rho+tau;
delta = rho-tau;
alfa = 1 - omega;  % Absorption coefficient for incoming diffuse radiation
I0 = s(i)*Rs;
%% numerical solver

Rbottom = gap(N)*(1-f)*I0;
solinit = bvpinit([0 LAI],[f*I0,f*I0]);
temp = bvp5c(@TwoStream,@bcfun,solinit,options_2);


sol = deval(temp,x);
Idn(:,i) = sol(1,:);
Iup(:,i) = sol(2,:);
R0(:,i) = gap*(1-f)*I0;

%% absorbed radiation
k  = spline(x0,k0,x).';  
A(:,i) =  alfa.*(Iup(:,i) + Idn(:,i) + k.*R0(:,i));
end



%% check energy budget closure
% s(1)*Rs - Iup(1,1) - (1-rg)*(Idn(end,1)+R0(end,1))
% trapz(x,A(:,1))
% % 
% s(2)*Rs - Iup(1,2) - (1-rg)*(Idn(end,2)+R0(end,2))
% trapz(x,A(:,2))

albedo = (Iup(1,1)+Iup(1,2))/Rs;
%% plottings
if nargout==0
subplot(121)
plot([Idn Iup R0],x)
set(gca,'ydir','reverse')
xlabel('radiation')
ylabel('canopy depth (LAI)')
subplot(122)
% plot(cumsum(A)*x(2),x)
plot(A,x)
set(gca,'ydir','reverse')
xlabel('absorbed radiation')
ylabel('canopy depth (LAI)')
end


%% two-stream equations
function dY = Direct(x,Y)
 
k  = spline(x0,k0,x);
dY(1,:) =  -k.*Y;

end

function J = jac1(x,~)
    k  = spline(x0,k0,x);
    J = -k;
end

function dY = TwoStream(x,Y)

	
k = spline(x0,k0,x);                    % Extinction coeficcient for direct radiation
J = spline(x0,J0,x);                    % J-function
gamma = 1/2*(omega+J*delta);            % Backward scattering coefficient for incoming diffuse radiation
sigma_bw = 1/2*omega*k + 1/2*J*delta;   % Backward scattering for direct radiation
sigma_fw = omega*k-sigma_bw;            % Farward scattering for direct radiation

R = deval(temp1,x)*(1-f)*I0;
dY(1,:) =  -(alfa + gamma).*Y(1,:) + gamma.*Y(2,:) + sigma_fw.*R;
dY(2,:) =   (alfa + gamma).*Y(2,:) - gamma.*Y(1,:) - sigma_bw.*R;

end

%% boundary conditions
function res = bcfun(ya,yb)
%     Rbottom = deval(temp1,LAI);
    res = [ya(1)-f*I0; yb(2) - rg*yb(1) - rg*Rbottom];
end   

%% Jacobians
function [dbcdya,dbcdyb] = fjacbc(~,~)
    
dbcdya=[1 0 
        0 0];

dbcdyb=[0   0
        -rg 1];

end
    
    
function Jac = fjac(x,~)

k = spline(x0,k0,x);                    % Extinction coeficcient for direct radiation
J = spline(x0,J0,x);                    % J-function
gamma = 1/2*(omega+J*delta);            % Backward scattering coefficient for incoming diffuse radiation

Jac(1,1) =  -alfa-gamma;
Jac(1,2) =  +gamma;
Jac(2,1) =  -gamma;
Jac(2,2) =   alfa+gamma;

end


end