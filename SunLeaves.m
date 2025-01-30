function [p,wm] = SunLeaves(x,ze,LAI,ME0)
% compute the precentage of sun exposed leaves in each layer a for each
% leaf inclination, according to a leaf angle distribution.
% The inclinations are binned in m classes from 0 to 1 |cos(r.rL)|, 
% where r.rL is the angle between the leaf normal and the sun direction.
% A beta leaf angle distribution that changes with canopy depth (vector x)
% is provided based on observations in BCI (Detto et al., 2015, JGR)
%
% INPUTS:
% x: vertical coordinate
% ze: solar zenit angle (rad)
% LAI: leaf area index
% ME0: mean angle inclination at the top (Deg)
%
% OUTPUTS:
% p:    is a n x m matrix with fraction of leaves in full sun, when n is number
%       of layers
% wm:   is the mean abolute valu of the cosine between the sun and leaf direction 
%       corresponding to p
    
m = 10;  
n = length(x);
z = 55*(1-x/LAI);
ME=ME0-(ME0-22)*exp(-(z/25).^2);
SD=20-(20-13)*exp(-(z/25).^2);

tbar = ME./90;
st = (SD/90).^2;
s0 = tbar.*(1-tbar);
nu = tbar.*(s0./st-1);
mu = (1-tbar).*(s0./st-1);


w = linspace(0,1,m+1);
fhi = linspace(0,pi,2^10);
u = linspace(0,1,2^10);
p = zeros(n,m);
wm = zeros(1,m);

for i=1:n
    th = betaincinv(u,nu(i),mu(i))*pi/2;
    
    [y,x]=meshgrid(fhi,th);
    
    z = abs(sin(x).*sin(ze).*cos(y)+cos(x).*cos(ze));
    
    for j=1:m
        use=z(:)>=w(j) & abs(z(:))<w(j+1);
        p(i,j) = mean(use);
        if p(i,j)>0
            wm(i,j) = mean(z(use));
        else
            wm(i,j) = w(j)+w(2)/2;
        end
    end
    
end


    
end