function [Ab,x,R0,Idn,Iup,albedo] = RT2S_analytical(LAI,Rs,ze,f,RHOL,varargin)


%% default parameters
param.GRoss = 1/2;
param.J = 1/3;

%% user-defined parameters
for i=1:(nargin-5)/2
    param.(varargin{2*(i-1)+1}) = varargin{2*i};
end
% sky conditions
% f = 0.5;
% ze = pi/4;
% I0 =  1; %total incoming radiation
L = LAI;

% G-ross and J functions
% if strcmp(LAD,'spherical')
%     GRoss = 1/2;
%     J = 1/3;
% elseif strcmp(LAD,'planophile')
%     
%     F = @(th) 2/pi*(1+cos(2*th));
%     GRoss = Gfunction(ze,F);
%     J = 3/4;
% elseif strcmp(LAD,'erectophile')
%    
%     F = @(th) 2/pi*(1-cos(2*th));
%     GRoss = Gfunction(ze,F);
%     J = 1/4;
% % else
% %     F = LeafAngleDistribution(LAD,[]);
% %     GRoss = Gfunction(ze,F);     % G-Ross function
% %     J = Jfunction(F);
% end

% spherical LAD
GRoss = param.GRoss;
J = param.J;
%% Visible
% rho = 0.10;    % leaf reflectance
rho = RHOL(1);    % leaf reflectance
tau = 0.05;    % leaf transmittance
rhos  = 0.15;  % soil reflectance
omega = rho+tau;
delta = rho-tau;

k = GRoss/cos(ze);                    % Extinction coeficcient for direct radiation
alfa = 1 - omega;                     % Absorption coefficient for incoming diffuse radiation
gamma = 1/2*(omega+J*delta);          % Backward scattering coefficient for incoming diffuse radiation
sbw = 1/2*omega*k + 1/2*J*delta;      % Backward scattering for direct radiation
sfw = omega*k-sbw;                    % Forward scattering for direct radiation
I0 = 0.50*Rs;


lambda = -sqrt(alfa^2+2*alfa*gamma); % eigenvalue of ODE matrix


%% albedo calculation
G = @(lambda) ((alfa+gamma+lambda-rhos*gamma)*sbw-(rhos*(alfa+gamma-lambda)-gamma)*sfw)/(k+lambda)*(1-f) + ...
    (gamma -rhos*(alfa+gamma-lambda))*f;
Q = (rhos + ((alfa+gamma-k)*sbw+gamma*sfw - rhos*((alfa+gamma+k)*sfw+gamma*sbw))/(k+lambda)/(k-lambda))*2*lambda*(1-f);

rhoc = (G(lambda)*exp(lambda*L)-G(-lambda)*exp(-lambda*L)+Q*exp(-k*L))/...
    ((alfa+gamma+lambda-rhos*gamma)*exp(lambda*L)-(alfa+gamma-lambda-rhos*gamma)*exp(-lambda*L));


%% coefficient of the solution
A = @(lambda) (((alfa+gamma+lambda)*sfw+gamma*sbw)/(k-lambda)*(1-f) + (alfa+gamma+lambda)*f - gamma*rhoc)/(2*lambda);
B = @(lambda) (((alfa+gamma-lambda)*sbw+gamma*sfw)/(k-lambda)*(1-f) + gamma*f - (alfa+gamma-lambda)*rhoc)/(2*lambda);
C = @(k,sfw,sbw) -((alfa+gamma+k)*sfw+gamma*sbw)/((k+lambda)*(k-lambda))*(1-f);


x = linspace(0,L,100);Y=zeros(length(x),3);
Y(:,1) = A(lambda)*exp(-lambda*x) + A(-lambda)*exp(lambda*x) + C(+k,sfw,sbw)*exp(-k*x);
Y(:,2) = B(lambda)*exp(-lambda*x) + B(-lambda)*exp(lambda*x) + C(-k,sbw,sfw)*exp(-k*x);
Y(:,3) = (1-f)*exp(-k*x);
Y = Y*I0;

%% check bottom BC and outputs
err = abs(1-rhos*(Y(100,1)+Y(100,3))./Y(100,2));
if err>1e-6
    display('error bottom boundary condition')
end

Ab(:,1) = alfa.*(Y(:,1) + Y(:,2) + k.*Y(:,3));
R0(:,1) = Y(:,3);
Idn(:,1) = Y(:,1);
Iup(:,1) = Y(:,2);
%% NIR
% rho = 0.40;  % leaf reflectance
rho = RHOL(2);    % leaf reflectance
tau = 0.25;  % leaf transmittance
rhos  = 0.15;  % soil reflectance
omega = rho+tau;
delta = rho-tau;

k = GRoss/cos(ze);                          % Extinction coeficcient for direct radiation
alfa = 1 - omega;                       % Absorption coefficient for incoming diffuse radiation 
gamma = 1/2*(omega+J*delta);            % Backward scattering coefficient for incoming diffuse radiation 
sbw = 1/2*omega*k + 1/2*J*delta;   % Backward scattering for direct radiation
sfw = omega*k-sbw;            % Farward scattering for direct radiation
I0 = 0.50*Rs;

lambda = -sqrt(alfa^2+2*alfa*gamma); % eigenvalue of ODE matrix


%% albedo calculation
G = @(lambda) ((alfa+gamma+lambda-rhos*gamma)*sbw-(rhos*(alfa+gamma-lambda)-gamma)*sfw)/(k+lambda)*(1-f) + ...
    (gamma -rhos*(alfa+gamma-lambda))*f;
Q = (rhos + ((alfa+gamma-k)*sbw+gamma*sfw - rhos*((alfa+gamma+k)*sfw+gamma*sbw))/(k+lambda)/(k-lambda))*2*lambda*(1-f);

rhoc = (G(lambda)*exp(lambda*L)-G(-lambda)*exp(-lambda*L)+Q*exp(-k*L))/...
    ((alfa+gamma+lambda-rhos*gamma)*exp(lambda*L)-(alfa+gamma-lambda-rhos*gamma)*exp(-lambda*L));


%% coefficient of the solution
A = @(lambda) (((alfa+gamma+lambda)*sfw+gamma*sbw)/(k-lambda)*(1-f) + (alfa+gamma+lambda)*f - gamma*rhoc)/(2*lambda);
B = @(lambda) (((alfa+gamma-lambda)*sbw+gamma*sfw)/(k-lambda)*(1-f) + gamma*f - (alfa+gamma-lambda)*rhoc)/(2*lambda);
C = @(k,sfw,sbw) -((alfa+gamma+k)*sfw+gamma*sbw)/((k+lambda)*(k-lambda))*(1-f);


x = linspace(0,L,100);Y=zeros(length(x),3);
Y(:,1) = A(lambda)*exp(-lambda*x) + A(-lambda)*exp(lambda*x) + C(+k,sfw,sbw)*exp(-k*x);
Y(:,2) = B(lambda)*exp(-lambda*x) + B(-lambda)*exp(lambda*x) + C(-k,sbw,sfw)*exp(-k*x);
Y(:,3) = (1-f)*exp(-k*x);
Y = Y*I0;

%% check bottom BC and outputs
err = abs(1-rhos*(Y(100,1)+Y(100,3))./Y(100,2));
if err>1e-6
    display('error bottom boundary condition')
end

Ab(:,2) = alfa.*(Y(:,1) + Y(:,2) + k.*Y(:,3));
R0(:,2) = Y(:,3);
Idn(:,2) = Y(:,1);
Iup(:,2) = Y(:,2);
albedo = (Iup(1,1)+Iup(1,2))/Rs;
%%
if nargout==0
subplot(121)
plot(Y,x);
set(gca,'ydir','reverse')
legend('downward diffuse','upward diffuse','downward direct','location','southeast')
xlabel('fration of incident radiation')
ylabel('canopy depth (cumulative leaf area)')

subplot(122)
Ab =  alfa.*(Y(:,1) + Y(:,2) + k.*Y(:,3));
plot(Ab,x);
set(gca,'ydir','reverse')
legend('absorbed radiation','location','southeast')
xlabel('fration of incident radiation')
ylabel('canopy depth (cumulative leaf area)')

end
% %% check with numerical integration
% [~,x2,~,Idn,Iup,albedo] = RT2S_SW(L,1,ze,f);
% 
% hold all
% plot([Idn(:,1) Iup(:,1)]*2,x2,'--');
% % %% chack integral
% x = linspace(0,L,10000);
% TT=zeros(2);
% for i=1:length(x)
%     TT = TT + expm(Z*(L-x(i)))*exp(-k*x(i))*x(2);
% end
% exp(-k*L)/(k+lambda)/(k-lambda)*(Z-k*eye(2))+...
% exp(lambda*L)/(2*lambda)/(k+lambda)*(Z+lambda*eye(2))-...
% exp(-lambda*L)/(2*lambda)/(k-lambda)*(Z-lambda*eye(2))
    

%% G-function
function GRoss = Gfunction(ze,F)


    if ze==0
        W1 = @(th) F(th).*cos(ze).*cos(th);
        
        GRoss = integral(@(x) W1(x),0,pi/2);
    else
        W1 = @(th) F(th).*cos(ze).*cos(th);
        W2 = @(th) F(th).*cos(ze).*cos(th) ...
            .*2./pi.*(sqrt(tan(th).^2.*tan(ze).^2-1)-asec(tan(th).*tan(ze)));

        GRoss = integral(@(x) W1(x),0,pi/2) + integral(@(x) W2(x),pi/2-ze,pi/2);

    end

end

%% J-Function
function J = Jfunction(F)

    W1 = @(th) F(th).*cos(th).^2;
    J = integral(@(x) W1(x),0,pi/2);
    
end

end