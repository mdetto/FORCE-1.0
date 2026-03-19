%% G-Ross and J-function verical distribution 
function [k,J,kstar,x,zL,zW] = TwoStreamCoeff(MA0,ze,LAI,fw,CI,m)

% Inputs:
%      LAI: leaf area index
%      I0:  total incoming radiation
%      ze:  solar zenit angle (rad)
%      MA0: mean angle at the top of canopy (Deg)
%      fw: fraction of woody elements
%      CI: clumping index (<1)

x = linspace(0,LAI);
n = length(x);

% leaf parameters mu and nu
ME = MA0*(1-0.105*x);
SD =  20*(1-0.058*x);

tbar = ME./90;
st=(SD/90).^2;
s0=tbar.*(1-tbar);
nu = tbar.*(s0./st-1);
mu=(1-tbar).*(s0./st-1);

% wood parameters mu and nu
ME = 80;
SD =  6;

tbar = ME./90;
st=(SD/90).^2;
s0=tbar.*(1-tbar);
nu_w = tbar.*(s0./st-1);
mu_w=(1-tbar).*(s0./st-1);


J = zeros(1,n);
k = zeros(1,n);
kstar = zeros(1,n);
for j=1:n
	%% erectophile woody element angle distribution
	F = @(th) (1-fw)*2./pi./beta(mu(j),nu(j))*(1-2*th/pi).^(mu(j)-1).*(2*th/pi).^(nu(j)-1) + ...
                  fw*2./pi./beta(mu_w,nu_w)*(1-2*th/pi).^(mu_w-1).*(2*th/pi).^(nu_w-1);
    % F = @(th) sin(th); %spherical distribution for testing
    k(j) = CI./OpticalDepth(F,ze);
    mubar = integral(@(x) OpticalDepth(F,x).*sin(x),0,pi/2);
    kstar(j) = CI./mubar;
    
	W2 = @(th) F(th).*cos(th).^2;
	J(j) = integral(@(x) W2(x),0,pi/2);
end


if nargout>4
  
% compute the precentage of sun exposed leaves in each layer a for each
% leaf inclination, according to a leaf angle distribution.
% The inclinations are binned in m classes from 0 to 1 |cos(r.rL)|, 
% where r.rL is the angle between the leaf normal and the sun direction.
% A beta leaf angle distribution that changes with canopy depth (vector x)
% is provided based on observations in BCI (Detto et al., 2015, JGR)
%
%
% OUTPUTS:
% p:    is a n x m matrix with fraction of leaves in full sun, when n is number
%       of layers
% wm:   is the mean abolute value of the cosine between the sun and leaf direction 
%       corresponding to p
    
fi = linspace(0,pi,2^10);
u = linspace(0,1,2^10).';
zL = zeros(100,m);

for i=1:n
    th = betaincinv(u,nu(i),mu(i))*pi/2;
    % th = acos(u);
    z = abs(cos(th).*cos(ze) + sin(th).*sin(ze).*cos(fi));
    q = quantile(z(:),0:1/m:1);
    for j=1:m
        use=z(:)>q(j) & z(:)<=q(j+1);
        zL(i,j) = mean(z(use));
    end

end

    th = betaincinv(u,nu_w,mu_w)*pi/2;
    % th = acos(u);
    z = abs(cos(th).*cos(ze) + sin(th).*sin(ze).*cos(fi));
    q = quantile(z(:),0:1/m:1);
    zW = zeros(1,m);
    for j=1:m
        use=z(:)>q(j) & z(:)<=q(j+1);
        zW(1,j) = mean(z(use));
    end
    
end

%% optical depth
function y = OpticalDepth(F,ze)
    y = zeros(size(ze));
    W0 = @(th) F(th).*cos(th);
    k0 = integral(@(x) W0(x),0,pi/2);
    for ii =1:length(ze)

        if ze(ii)==0
         y(ii)=1./k0;

        elseif ze(ii)==pi/2
             y (ii) = 0;

        else
        	W1 = @(th) F(th).*cos(th).*2./pi.*(sqrt(tan(th).^2.*tan(ze(ii)).^2-1)...
                -asec(tan(th).*tan(ze(ii))));

        	y(ii) = 1./(k0 + integral(@(x) W1(x),pi/2-ze(ii),pi/2));
  
        end
    end

end

end
% %speherical leaf angle distribution used for testing
% x0 = linspace(0,LAI);
% J0 = 1/3*ones(100,1);
% k0 = 1/2/cos(ze)*ones(100,1);

% gamma = 1/2*(omega+J*delta);           % Backward scattering coefficient for incoming diffuse radiation
