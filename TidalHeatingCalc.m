function [H H_total] = TidalHeatingCalc(r,T,por,eta_diss,lat_inp,lon_inp,etaflag)

% will put in a flag so heating can be calculated as an average or
% at a specific lat and long.

% Input arguments required
% temperature and porosity (used to get viscosity)
% passed parameters must be dimensional
% Give a colatitude and longitude if a specific value is requested. DEGREES
% flag = 1 for normal calculation, flag = 2 for asthenoshere heating

% references:
% S&V = Sabadini and Vermeersenm (2004)
% R&N = Roberts and Nimmo (2008)
% Tobie = Tobie et al (2005)
% B&N = Bierson and Nimmo (2016)

if nargin<5
    lat_inp = nan;
    lon_inp = nan;
    etaflag = 1;
end

if isrow(r); r = r'; end
if isrow(T); T = T'; end
if isrow(por); por = por'; end

% downsample arrays if need
fullsize = length(r);
while length(r)>100
    r = downsample(r,10);
    T = downsample(T,10);
    por = downsample(por,10);
end

% pameters
% NOTE the viscosity is a dissipation viscosity, separate to the dynamics
Nrad = length(r); % number of grid points between core and surface
Nlat = 50; % number of grid points in colatitude
Nlon = 100; % number of grid points in longitude
lmax = 2;    % highest spherical harmonic degree
rho0 = 3000; % (kg/m^3) 3000, reference mantle density. Segatz benchmark 3200
rhoc = 7640;     % core density 7640, Segatz benchmark 5150 
omega = 4.11e-5; % (per second) orbital freq
period = 2*pi/omega; % (seconds) orbital period 
ecc = 0.0041; % eccentricity io = 0.004
G = 6.67e-11; % (m^3 kg^-1 s^-2) gravitational constant
R = r(end); % metres, plantary radius
rm = r(1); % metres 700e3, core radius. Segatz benchmark 980
muel0 = 5e10; % reference rigidity (G0 in B&N). Ross fig2a benchmark 5e^10
P1 = [0, 0, 1, 0, 0, 0; 0, 0, 0, 1, 0, 0; 0, 0, 0, 0, 0, 1]; % pull out 3, 4 & 6 entry in vector
P02 = @(x) 0.5*(3*x.^2 -1);
P22 = @(x) 3*(1-x.^2);

if exist('eta_diss','var') == 0
    eta_diss = 7e15; % reference viscosity. Ross fig2a benchmark 2e16
end

if etaflag == 2
    eta_diss_vec = ones(length(r),1);
    crust_ind = find(~por(2:end),1);
    eta_diss_vec(crust_ind:end) = eta_diss*1e6;
    eta_diss_vec(crust_ind-30:crust_ind-1) = eta_diss;
    eta_diss_vec(1:crust_ind-31) = eta_diss*1000;
    eta_diss = eta_diss_vec;
end

if etaflag == 3
    eta_diss_vec = ones(length(r),1);
    crust_ind = find(~por(2:end),1);
    eta_diss_vec(crust_ind:end) = eta_diss*1e6;
    eta_diss_vec(crust_ind-30:crust_ind-1) = eta_diss;
    eta_diss_vec(1:crust_ind-31) = eta_diss*50;
    eta_diss = eta_diss_vec;
end

EA = 3e5; % activation energy (J/mol)
Rg = 8.3145; % gas constant
lamda = 27; % porosity factor for viscosity
Teta = 1500; % Reference temperature for viscosity, take melting point
c_cst = 67/15; % constant for porosity dependent rigidity (eq 27 B&N)

% Get grids
theta = linspace(0.5*pi/Nlat,pi-0.5*pi/Nlat,Nlat)';
phi = linspace(0.5*2*pi/Nlon,2*pi-0.5*2*pi/Nlon,Nlon)';
tvec = linspace(0,7*period/8,8)';

dr = r(2) - r(1); % radial grid spacing

% vector parameters
% Viscosity is a function of temperature and porosity (Katz 2010)
% Rigidity is a function of porosity (B&N)
rho = rho0*ones(Nrad,1);
if etaflag ~= 2
    eta = eta_diss.*exp(EA/Rg *(1./T - 1./Teta) - lamda*por);
else
    eta = eta_diss;
end
muel = muel0; % removed porosity dependence of muel, let it stay constant %./(1+c_cst*por);
mu = 1./(1./muel - i*1./(omega*eta)); % Maxwell complex shear modulus
K = rho*8000^2 - 4*muel/3; % (Pa) bulk modulus, middle value examined in Tobie
lam = K - 2*mu/3; % 2nd Lame parameter 

g = zeros(Nrad,1);
g(1) = (4/3)*pi*G*rhoc*rm;
for j=2:Nrad
    g(j) = (4/3)*pi*G*rho(j)*r(j); % get gravity to current layer assuming constant density
    if j>2
        for ii=j-1:1
            g(j) = g(j) + (4/3)*pi*G*(rho(ii)-rho(ii+1))*r(ii)^3/r(j)^2; % add contribution from different density of each underlying shell
        end
    end
    g(j) = g(j) + (4/3)*pi*G*(rhoc-rho(1))*rm^3/r(j)^2; % add contribution from the core
end

%allocate space for propagator
prop = zeros(6,3,Nrad);
ysoln_2 = zeros(6,Nrad);
ysoln_3 = zeros(6,Nrad);
ysoln_4 = zeros(6,Nrad);

% Main loop
      for l = 2:lmax; % sph degree
          IC = IC_matrix(rm,l,g(1),rhoc); %IC matrix
          prop(:,:,1) = IC;
          % evaluate the propagator up to the surface so as to calculate C.
          for ii = 1:Nrad-1
             Y = Y_matrix(r(ii+1),l,rho(ii+1),g(ii+1),mu(ii+1)); %fundamental matrix for layer
             Yinv = Yinverse_matrix(r(ii),l,rho(ii+1),g(ii+1),mu(ii+1)); %inverse fundamental matrix for layer, evaluated at base of layer
             prop(:,:,ii+1) = Y*(Yinv*prop(:,:,ii)); % product of propagator from below, Yinv of new layer but at base, and Y of new layer at top
          end
          % extract the 3rd, 4th and 6th rows, for which we have surface boundary conditions
          AA = P1*prop(:,:,end);
          % solve 'propagator * C = boundary conds' for C
          C = AA\b_vector(l);
          % solution at the CMB
          ysoln_2(:,1) = prop(:,:,1)*C;
          for ii = 1:Nrad-1
              Y = Y_matrix(r(ii+1),l,rho(ii+1),g(ii+1),mu(ii+1)); %fundamental matrix
              Yinv = Yinverse_matrix(r(ii),l,rho(ii+1),g(ii+1),mu(ii+1)); %inverse fundamental matrix
              ysoln_2(:,ii+1) = (Y*(Yinv*prop(:,:,ii)))*C;
          end
      end

% Use the sum of radial functions with SHD 
y1 = ysoln_2(1,:)';
y2 = ysoln_2(2,:)';
y3 = ysoln_2(3,:)';
y4 = ysoln_2(4,:)';
y5 = ysoln_2(5,:)';
y6 = ysoln_2(6,:)';

% Tidal potential and derivative functions. NOTE this will give an error at the moment if vectors are input, look at sizes of omega, t etc
pot = @(t,theta,phi) R^2*omega^2*ecc*(-3/2*P02(cos(theta)).*cos(omega*t) + 0.25*P22(cos(theta)).*(3*cos(omega*t).*cos(2*phi) + 4*sin(omega*t).*sin(2*phi)));
dpotdth = @(t,theta,phi) 3*R^2*omega^2*ecc*sin(theta).*cos(theta).*cos(phi).*(3*cos(phi).*cos(omega*t) + 4*sin(phi).*sin(omega*t));
dpotdph = @(t,theta,phi) 0.75*R*R*omega*omega*ecc*sin(theta).*sin(theta).*(-7*sin(2*phi - omega*t) + sin(2*phi + omega*t));
d2potd2th = @(t,theta,phi) -1.5*R*R*omega*omega*ecc*cos(2*theta).*cos(phi).*(-7*cos(phi-omega*t)+cos(phi+omega*t));
d2potd2ph = @(t,theta,phi) 1.5*R*R*omega*omega*ecc*sin(theta).*sin(theta).*(-7*cos(2*phi-omega*t)+cos(2*phi+omega*t));
d2potdtp = @(t,theta,phi) -3*R*R*omega*omega*ecc*cos(theta).*sin(theta).*(3*sin(2*phi).*cos(omega*t)-4*cos(2*phi).*sin(omega*t));

dy1dr = -2*lam.*y1./((lam+2*mu).*r) + y3./(lam+2*mu) + lam*l*(l+1).*y2./((lam+2*mu).*r);

% if specific latitudes and longitudes are requested
if isnan(lat_inp) == 0 && isnan(lon_inp) == 0
    lat_inp = lat_inp*pi/180;
    lon_inp = lon_inp*pi/180;
    % Meshgrid Theta and Phi onto different potentials
    % correction to R&N documentation, y3 in there should be y2 and vice versa
    [Dy1dr,tgrid] = ndgrid(dy1dr,tvec);
    e11 = Dy1dr.*pot(tgrid,lat_inp,lon_inp);

    [y2r,tgrid] = ndgrid(y2./r,tvec);
    e22 = y2r.*d2potd2th(tgrid,lat_inp,lon_inp);
    [y1r,tgrid] = ndgrid(y1./r,tvec);
    e22 = e22 + y1r.*pot(tgrid,lat_inp,lon_inp);

    e33 = y1r.*pot(tgrid,lat_inp,lon_inp);
    [y2term,tgrid] = ndgrid(y2./r,tvec);
    e33 = e33 + y2term.*d2potd2ph(tgrid,lat_inp,lon_inp)./(sin(lat_inp).^2);
    [y2term,tgrid] = ndgrid(y2./r,tvec);
    e33 = e33 + y2term.*cot(lat_inp).*dpotdth(tgrid,lat_inp,lon_inp);

    [y4term,tgrid] = ndgrid(y4./mu,tvec);
    e12 = 0.5*y4term.*dpotdth(tgrid,lat_inp,lon_inp);

    [y4term,tgrid] = ndgrid(y4./mu,tvec);
    e13 = 0.5*y4term./sin(lat_inp) .*dpotdph(tgrid,lat_inp,lon_inp);

    [y2term,tgrid] = ndgrid(y2./r,tvec);
    e23 = y2term./sin(lat_inp) .*(d2potdtp(tgrid,lat_inp,lon_inp) - cot(lat_inp).*dpotdph(tgrid,lat_inp,lon_inp));

    [lamgrid,tgrid] = ndgrid(lam,tvec);
    hydro = lamgrid.*(e11 + e22 + e33);
    [mugrid,tgrid] = ndgrid(mu,tvec);
    str11 = 2*mugrid.*e11 + hydro;
    str22 = 2*mugrid.*e22 + hydro;
    str33 = 2*mugrid.*e33 + hydro;
    str12 = 2*mugrid.*e12;
    str13 = 2*mugrid.*e13;
    str23 = 2*mugrid.*e23;

    heat = 0.5*omega*(imag(str11).*real(e11) - real(str11).*imag(e11));
    heat = heat + 0.5*omega*(imag(str22).*real(e22) - real(str22).*imag(e22));
    heat = heat + 0.5*omega*(imag(str33).*real(e33) - real(str33).*imag(e33));
    heat = heat + omega*(imag(str12).*real(e12) - real(str12).*imag(e12));
    heat = heat + omega*(imag(str13).*real(e13) - real(str13).*imag(e13));
    heat = heat + omega*(imag(str23).*real(e23) - real(str23).*imag(e23));

    heatavg = sum(heat,2)/8; % get heat averaged over orbital cycle
    heatint = squeeze(sum(heatavg,1)*dr);

    H = heatavg;
    
    H_total = nan; % can't evulate total heating from 1 lat and lon

    % put back on to full length array if need
    while length(H)<fullsize
        H = interp1(r,H,linspace(700e3,1820e3,fullsize),'linear','extrap')';
        H(H<0) = 0;
    end
    return
end

% Meshgrid Theta and Phi onto different potentials
% correction to R&N documentation, y3 in there should be y2 and vice versa
[Dy1dr,Theta,Phi,tgrid] = ndgrid(dy1dr,theta,phi,tvec);
e11 = Dy1dr.*pot(tgrid,Theta,Phi);

[y2r,Theta,Phi,tgrid] = ndgrid(y2./r,theta,phi,tvec);
e22 = y2r.*d2potd2th(tgrid,Theta,Phi);
[y1r,Theta,Phi,tgrid] = ndgrid(y1./r,theta,phi,tvec);
e22 = e22 + y1r.*pot(tgrid,Theta,Phi);

e33 = y1r.*pot(tgrid,Theta,Phi);
[y2term,Theta,Phi,tgrid] = ndgrid(y2./r,theta,phi,tvec);
e33 = e33 + y2term.*d2potd2ph(tgrid,Theta,Phi)./(sin(Theta).^2);
[y2term,Theta,Phi,tgrid] = ndgrid(y2./r,theta,phi,tvec);
e33 = e33 + y2term.*cot(Theta).*dpotdth(tgrid,Theta,Phi);

[y4term,Theta,Phi,tgrid] = ndgrid(y4./mu,theta,phi,tvec);
e12 = 0.5*y4term.*dpotdth(tgrid,Theta,Phi);

[y4term,Theta,Phi,tgrid] = ndgrid(y4./mu,theta,phi,tvec);
e13 = 0.5*y4term./sin(Theta) .*dpotdph(tgrid,Theta,Phi);

[y2term,Theta,Phi,tgrid] = ndgrid(y2./r,theta,phi,tvec);
e23 = y2term./sin(Theta) .*(d2potdtp(tgrid,Theta,Phi) - cot(Theta).*dpotdph(tgrid,Theta,Phi));

[lamgrid,Theta,Phi,tgrid] = ndgrid(lam,theta,phi,tvec);
hydro = lamgrid.*(e11 + e22 + e33);
[mugrid,Theta,Pih,tgrid] = ndgrid(mu,theta,phi,tvec);
str11 = 2*mugrid.*e11 + hydro;
str22 = 2*mugrid.*e22 + hydro;
str33 = 2*mugrid.*e33 + hydro;
str12 = 2*mugrid.*e12;
str13 = 2*mugrid.*e13;
str23 = 2*mugrid.*e23;

heat = 0.5*omega*(imag(str11).*real(e11) - real(str11).*imag(e11));
heat = heat + 0.5*omega*(imag(str22).*real(e22) - real(str22).*imag(e22));
heat = heat + 0.5*omega*(imag(str33).*real(e33) - real(str33).*imag(e33));
heat = heat + omega*(imag(str12).*real(e12) - real(str12).*imag(e12));
heat = heat + omega*(imag(str13).*real(e13) - real(str13).*imag(e13));
heat = heat + omega*(imag(str23).*real(e23) - real(str23).*imag(e23));

heatavg = sum(heat,4)/8; % get heat averaged over orbital cycle
heatint = squeeze(sum(heatavg,1)*dr);

% Volume of wedges used to calculate average heating
V_wedge = 2*pi*(R^3-r(1)^3)/(3*(Nlon-1)) *(-cos(theta(2:end)) + cos(theta(1:end-1)));
for j=1:Nlat-1
    heatavg_mid_lat(:,j,:) = 0.5*(heatavg(:,j,:)+heatavg(:,j+1,:));
end
for k=1:Nlon-1
    heatavg_mid(:,:,k) = 0.5*(heatavg_mid_lat(:,:,k)+heatavg_mid_lat(:,:,k+1));
end

[V_mesh r_mesh lon_mesh] = meshgrid(V_wedge,r,phi(1:end-1));
heatrad = sum(squeeze(sum(heatavg_mid.*V_mesh/(4/3*pi*(R^3-r(1)^3)),2)),2);

H = heatrad;
H(H<0) = 0;

% need to integrate heating over the whole body
for ii=1:Nrad-1
    heatavg_midr(ii,:,:) = 0.5*(heatavg_mid(ii,:,:)+heatavg_mid(ii+1,:,:));
end

load 'V_mesh.mat';
H_total = sum(heatavg_midr.*V_mesh,'all');

% put back on to full length array if need
while length(H)<fullsize
    H = interp1(r,H,linspace(700e3,1820e3,fullsize),'linear','extrap')';
end

% For benchmarking against Ross fig 2a
% imagesc(phi*180/pi,theta*180/pi,heatint);
% hold on
% [C1,h1] = contour(gca,phi*180/pi,theta*180/pi,heatint,[2 2],'linewidth',2,'Color','white');
% [C1,h1] = contour(gca,phi*180/pi,theta*180/pi,heatint,[1.5 1.5],'linewidth',2,'Color','white');
% [C1,h1] = contour(gca,phi*180/pi,theta*180/pi,heatint,[1 1],'linewidth',2,'Color','white');

 end  % closes main function

      
%%%%%%%%%% subroutines %%%%%%%%

function IC = IC_matrix(r,l,g,rhoc)
   % from S&V eqn. (2.6)
   G = 6.67e-11;
   
   IC = zeros(6,3);
   IC(1,1) = -r^l/g;
   IC(1,3) = 1;
   IC(2,2) = 1;
   IC(3,3) = g*rhoc;
   IC(5,1) = r^l;
   IC(6,1) = 2*(l-1)*r^(l-1);
   IC(6,3) = 4*pi*G*rhoc;
end
   
function Y = Y_matrix(r,l,rho,g,mu);
   % from S&V eqn. (2.42)
   G = 6.67e-11; % (m^3 kg^-1 s^-2) gravitational constant
   
   Y = zeros(6,6); % (row,col)
   Y(1,1) = l*r^(l+1)/(2*(2*l+3));
   Y(1,2) = r^(l-1);
   Y(1,4) = (l+1)*r^(-l)/(2*(2*l-1));
   Y(1,5) = r^(-l-2);
   
   Y(2,1) = (l+3)*r^(l+1)/(2*(2*l+3)*(l+1));
   Y(2,2) = r^(l-1)/l;
   Y(2,4) = (2-l)*r^(-l)/(2*l*(2*l-1));
   Y(2,5) = -r^(-l-2)/(l+1);
   
   Y(3,1) = (l*rho*g*r+ 2*(l^2-l-3)*mu)*r^l/(2*(2*l+3));
   Y(3,2) = (rho*g*r+ 2*(l-1)*mu)*r^(l-2);
   Y(3,3) = rho*r^l; % -ve in R&N code
   Y(3,4) = ((l+1)*rho*g*r- 2*(l^2+3*l-1)*mu)/(2*(2*l-1)*r^(l+1));
   Y(3,5) = (rho*g*r- 2*(l+2)*mu)/r^(l+3);
   Y(3,6) = rho/r^(l+1); % -ve in R&N code 
   
   Y(4,1) = l*(l+2)*mu*r^l/((2*l+3)*(l+1));
   Y(4,2) = 2*(l-1)*mu*r^(l-2)/l;
   Y(4,4) = (l^2-1)*mu/(l*(2*l-1)*r^(l+1));
   Y(4,5) = 2*(l+2)*mu/((l+1)*r^(l+3));
   
   Y(5,3) = r^l; % -ve in R&N code
   Y(5,6) = 1/r^(l+1); % -ve in R&N code
   
   Y(6,1) = 2*pi*G*rho*l*r^(l+1)/(2*l+3);
   Y(6,2) = 4*pi*G*rho*r^(l-1);
   Y(6,3) = (2*l+1)*r^(l-1); % -ve in R&N code
   Y(6,4) = 2*pi*G*rho*(l+1)/((2*l-1)*r^l);
   Y(6,5) = 4*pi*G*rho/r^(l+2);
end
   
function Y = Ybar_matrix(r,l,rho,g,mu)
   % from S&V eqn. (2.47)
   G = 6.67e-11; % (m^3 kg^-1 s^-2) gravitational constant
   
   Y = zeros(6,6); 
   Y(1,1) = (rho*g*r)/mu - 2*(l+2);
   Y(1,2) = 2*l*(l+2);
   Y(1,3) = -r/mu;
   Y(1,4) = l*r/mu;
   Y(1,5) = rho*r/mu;
   
   Y(2,1) = -rho*g*r/mu + 2*(l^2+3*l-1)/(l+1);
   Y(2,2) = -2*(l^2-1);
   Y(2,3) = r/mu;
   Y(2,4) = (2-l)*r/mu;
   Y(2,5) = -rho*r/mu;
   
   Y(3,1) = 4*pi*G*rho;
   Y(3,6) = -1;
   
   Y(4,1) = rho*g*r/mu + 2*(l-1);
   Y(4,2) = 2*(l^2-1);
   Y(4,3) = -r/mu;
   Y(4,4) = -(l+1)*r/mu;
   Y(4,5) = rho*r/mu;
   
   Y(5,1) = -rho*g*r/mu - 2*(l^2-l-3)/l;
   Y(5,2) = -2*l*(l+2);
   Y(5,3) = r/mu;
   Y(5,4) = (l+3)*r/mu;
   Y(5,5) = -rho*r/mu;
   
   Y(6,1) = 4*pi*G*rho*r;
   Y(6,5) = (2*l)+1;
   Y(6,6) = -r;
end
   
function D = D_matrix(r,l)
   % from S&V eqn. (2.46)
   V = [(l+1)/r^(l+1);
       l*(l+1)/(2*(2*l-1)*r^(l-1));
       -1/r^(l-1); % +ve in R&N code
       l*r^l;
       l*(l+1)*r^(l+2)/(2*(2*l+3));
       r^(l+1)]; % -ve in R&N code
   V = V/(2*l+1);
   D = spdiags(V,0,6,6);
end
   
function Y = Yinverse_matrix(r,l,rho,g,mu)
   % from S&V eqn. (2.45)
   Y = D_matrix(r,l) * Ybar_matrix(r,l,rho,g,mu);
end
   
function b = b_vector(l)
   % from S&V eqn. 1.128 & 1.130
   R = 1820e3; % metres, plantary radius (Io value = 1821500)
   
   b = [0; 0; (2*l+1)/R]; % -ve in R&N code
end

function plot_func = plot_func(x,y)
    figure('Units','centimeters','Position',[15 15 12 30],'PaperPositionMode','auto');
    line = plot(x,y/1000);
    hold on
    line.Color = [1 0 0];
    line.LineWidth = 2.000;
    ymin = 0.6e6/1000;
    ymax = 2e6/1000;
    xmin = 0;
    xmax = 2000;

    line2 = plot(linspace(xmin,xmax,2),700*ones(2)); %core line
    line3 = plot(linspace(xmin,xmax,2),1821.5*ones(2)); %surface line
    line4 = plot(zeros(2),linspace(ymin,ymax,2));
    line2(1).Color = [0 0 1];
    line3(1).Color = [0 0 1];
    line2(1).LineWidth = 1.500;
    line3(1).LineWidth = 1.500;
    
    axis([xmin xmax ymin ymax])
    set(gca,'Units','normalized','YTick',ymin:(ymax-ymin)/7:ymax,'XTick',xmin:(xmax-xmin)/4:xmax,'Position',[.2 .2 .75 .7],'FontUnits','points','FontSize',16,'FontName','Times')
    
    ylabel({'Depth $ (km) $'},'FontUnits','points', 'interpreter','latex','FontSize', 20,'FontName','Times')
    xlabel({'$ y_{5} ~(kg/m^{3}) $'},'FontUnits','points','interpreter','latex','FontSize', 20,'FontName','Times')
    
    title('(e) Potential', 'FontUnits','points','Fontweight','normal','FontSize',20,'FontName','Times')
end