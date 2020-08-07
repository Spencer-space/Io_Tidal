function IoTidal_results = IoTidal_script(lam_M,lam_q,eta_flag)

% IoTidal iterator
% Calculates an initial structure using IoTidal with uniform tidal heating.
% Feeds this structure to TidalHeatingCalc to get a new heating profile.
% Heating profile put back into IoTidal, and iterated unti; IoTidal stops changing.

% lam_M = constant emplacement rate, units s^-1
% lam_q = emplacement proportionality to qp, no unit

% If not specified, don't load from file
if nargin < 6
    load_from_file = 0;
end

if eta_flag == 1
    filename = "IoTidal_mant_lamM_" + num2str(lam_M) + "_lamq_" + num2str(lam_q);
else
    filename = "IoTidal_asth_lamM_" + num2str(lam_M) + "_lamq_" + num2str(lam_q);
end

% eta_flag = 1 calculates tidal heating from porosity temperature
% eta_flag = 2 calculates with a high dissipation layer (low visc).
% eta_diss values give total heating ~ 1e14 W
if eta_flag == 1
    eta_diss = 2.7e15;
elseif eta_flag == 2
    eta_diss = 1.7e13;
end

Nlat = 50;
Nlon = 100;

% Dimensional Parameters
K_0 = 10^-7; % reference permeability (m^2)
rho0 = 3000; % density (kg/m^3)
rhoc = 7640; % core density (kg/m^3)
del_rho = 500; % density difference (kg/m^3)
rhol = rho0 - del_rho; % liquid density (kg/m^3)
G = 6.67e-11; % (m^3 kg^-1 s^-2) gravitational constant
g = 1.5; % gravity (m/s^2)
L = 4e5; % latent heat (J/kg)
kappa = 1e-6; % Thermal diffusivity (m^2/s)
c = 1200; % specific heat (J/kg/K)
T_l = 1350; % relative melting point (above surface T) (K)
T_surf = 150; % surface temperature (K)
T0 = T_l + T_surf; % Melting point, reference temp for density (K)
n = 3; % permeability exponent
eta_l = 1; % basalt melt viscosity (Pas)
eta = 1e20; % viscosity (Pas)
alpha = 3e-5; % Coef of thermal expansion (K^-1)

r_s = 1820e3; % Io radius (m)
r_m = 700e3; % core radius

Pc = 0; % not exploring Pc so set to zero

lat_inp = linspace(0,180,Nlat);
lon_inp = linspace(0,360,Nlon);
lat_mids = 0.5*(lat_inp(1:end-1)+lat_inp(2:end));
lon_mids = 0.5*(lon_inp(1:end-1)+lon_inp(2:end));
r = linspace(r_m,r_s,1000);

it = 1; % initialise iteration counter
% try to load from file, if doesn't exist then do the calc
if ~isfile(filename+".mat")
    S = IoTidal(lam_M,lam_q,Pc,1,0,0); % initially a small Pc is used
    [H H_tot] = TidalHeatingCalc(S.r,S.T,S.phi,eta_diss,NaN,NaN,eta_flag);
    Hnew = 0; %initialise Hnew
    
    % don't iterate if doing low viscosity layer
    if eta_flag == 2
        Hnew = H;
        Hnew_tot = H_tot;
        fprintf('Total heating = %.3e \n',Hnew_tot);
        fprintf('Norm(H-Hold) = %.3e \n',norm(Hnew-H));
    end
    
    while norm(Hnew-H)>1e-7
        if it>1
            H = Hnew;
        end
        S = IoTidal(lam_M,lam_q,Pc,H,H_tot,0);
        H_store(:,it) = H;
        S_store(it) = S;
        [Hnew Hnew_tot] = TidalHeatingCalc(S.r,S.T,S.phi,eta_diss,NaN,NaN,eta_flag);
        fprintf('it = %i \n',it);
        fprintf('Total heating = %.3e \n',Hnew_tot);
        fprintf('Norm(H-Hold) = %.3e \n',norm(Hnew-H));
        it = it+1;
    end
    H_avg = H;
    H_tot = Hnew_tot;
    S_avg = IoTidal(lam_M,lam_q,Pc,H_avg,H_tot,1);
    
    heat_matrix = zeros(1000,Nlat-1,Nlon-1);

    % Now go over latitude and longitude
    for i = 1:Nlat-1
        lat = lat_mids(i);
        parfor j = 1:Nlon-1
            lon = lon_mids(j);
            [H H_tot_nan] = TidalHeatingCalc(S_avg.r,S_avg.T,S_avg.phi,eta_diss,lat,lon,eta_flag);
            S = IoTidal(lam_M,lam_q,Pc,H,H_tot,0);
            heat_matrix(:,i,j) = H;
            por_matrix(:,i,j) = S.phi;
            q_matrix(:,i,j) = S.q;
            T_matrix(:,i,j) = S.T;
            crust_thick(i,j) = S.rl;
            P_matrix(:,i,j) = S.P
            erupt_rate(i,j) = S.qe(end)*(100*60*60*24*365);
            cond_heat(i,j) = -kappa*rho0*c* S.dT(end)*T_l/r_s * (4*pi*r_s^2); % global conductive heat flux (W) kappa*rho*c = conductivity
            volc_heat(i,j) = rho0*S.qe(end) *(4*pi*r_s^2)*(c*T_l + L);
            fprintf('Lat = %i / %i, Lon = %i / %i \n',i,Nlat-1,j,Nlon-1);
        end
    end
else
    load(filename);
    Pc = IoTidal_results.Pc;
    density_matrix = IoTidal_results.density_matrix;
    eta_diss = IoTidal_results.eta_diss;
    H_avg = IoTidal_results.avg_heat;
    S_avg = IoTidal_results.avg_soln;
    H_tot = IoTidal_results.tot_heat;
    por_matrix = IoTidal_results.por_matrix;
    T_matrix = IoTidal_results.T_matrix;
    heat_matrix = IoTidal_results.heat_matrix;
    q_matrix = IoTidal_results.q_matrix;
    P_matrix = IoTidal_results.P_matrix;
    crust_thick = IoTidal_results.crust_thick;
    erupt_rate = IoTidal_results.erupt_rate;
    cond_heat = IoTidal_results.cond_heat;
    volc_heat = IoTidal_results.volc_heat;
end

% Isostasy calculation
density_matrix = flip(rho0.*(1 - alpha*(T_matrix - T0))); % gives density as a function down from surface
density_matrix(density_matrix == rho0) = 0; % set values in mantle to zero
z = linspace(0,1820e3-700e3,1000); % get position vector as a function down from surface
B = squeeze(trapz(z,density_matrix/rho0 - 1.*(density_matrix>0))); % get buoyancy matrix
% Find the average buoyancy, weighted by SA
for i = 1:Nlat-1
    SA_func = @(theta,phi) r_s^2*sin(theta);
    SA(i) = integral2(SA_func,lat_inp(i)*pi/180,lat_inp(i+1)*pi/180,lon_inp(1)*pi/180,lon_inp(2)*pi/180);
end
[~,SA_mesh] = meshgrid(lon_mids,SA); % mesh surface area with longitude (uniform in longitude), should use repmat
B_mean = sum(B.*SA_mesh,'all')/(4*pi*r_s^2); % get the SA weighted averaged buoyancy
topography = -(B-B_mean); % get topography, average topography is zero

IoTidal_results.cond_heat = cond_heat;
IoTidal_results.volc_heat = volc_heat;
IoTidal_results.erupt_rate = erupt_rate;
IoTidal_results.crust_thick = crust_thick;
IoTidal_results.por_matrix = por_matrix;
IoTidal_results.q_matrix = q_matrix;
IoTidal_results.P_matrix = P_matrix;
IoTidal_results.T_matrix = T_matrix;
IoTidal_results.heat_matrix = heat_matrix;
IoTidal_results.density_matrix = density_matrix;
IoTidal_results.topography = topography;
IoTidal_results.lat = lat_mids;
IoTidal_results.lon = lon_mids;
IoTidal_results.S_avg = S_avg;
IoTidal_results.Pc = Pc;
IoTidal_results.info = 'Emplacement rate M = lam_M + lam_q*qp';
IoTidal_results.eta_diss = eta_diss;
IoTidal_results.avg_soln = S_avg;
IoTidal_results.avg_heat = H_avg;
IoTidal_results.tot_heat = H_tot;

save(filename,'IoTidal_results');

figure('Units','centimeters','Position',[15 15 70 20],'PaperPositionMode','auto');
fig1 = subplot(1,3,1);
imagesc(fig1,IoTidal_results.lon,IoTidal_results.lat,IoTidal_results.crust_thick);
set(fig1,'Ytick',0:45:180,'Xtick',0:90:360,'Units','normalized','FontUnits','points','FontSize',26,'FontName','Times','ydir','normal');
title(fig1,'a) Crustal Thickness (km)','position',[105 -1], 'FontUnits','points','Fontweight','normal','interpreter','latex','FontSize', 26,'FontName','Times');
xlabel(fig1,{'Longitude'},'FontUnits','points','interpreter','latex','FontSize', 30,'FontName','Times');
ylabel(fig1,{'Colatitude'},'FontUnits','points','interpreter','latex','FontSize', 30,'FontName','Times');
fig1.YDir = 'reverse';
axis([0 360 0 180]);
colorbar
colormap(fig1,bone);

fig2 = subplot(1,3,2);
imagesc(fig2,IoTidal_results.lon,IoTidal_results.lat,IoTidal_results.erupt_rate);
set(fig2,'Ytick',0:45:180,'Xtick',0:90:360,'Units','normalized','FontUnits','points','FontSize',26,'FontName','Times','ydir','normal');
title(fig2,'b) Eruption Rate (cm/yr)','position',[100 -1], 'FontUnits','points','Fontweight','normal','interpreter','latex','FontSize', 26,'FontName','Times');
xlabel(fig2,{'Longitude'},'FontUnits','points','interpreter','latex','FontSize', 30,'FontName','Times');
set(gca,'Yticklabel',[]);
fig2.YDir = 'reverse';
axis([0 360 0 180]);
colorbar
colormap(fig2,pink);

fig3 = subplot(1,3,3);
imagesc(fig3,IoTidal_results.lon,IoTidal_results.lat,IoTidal_results.topography/1000);
hold on
[C1,h1] = contour(fig3,IoTidal_results.lon,IoTidal_results.lat,IoTidal_results.topography/1000,[0 0],'linewidth',2,'Color','white');
[C2,h2] = contour(fig3,IoTidal_results.lon,IoTidal_results.lat,IoTidal_results.topography/1000,[0.1 0.1],'linewidth',2,'Color','black');
[C4,h4] = contour(fig3,IoTidal_results.lon,IoTidal_results.lat,IoTidal_results.topography/1000,[-0.1 -0.1],'linewidth',2,'Color','white');
[C5,h5] = contour(fig3,IoTidal_results.lon,IoTidal_results.lat,IoTidal_results.topography/1000,[-0.25 -0.25],'linewidth',2,'Color','white');
[C6,h6] = contour(fig3,IoTidal_results.lon,IoTidal_results.lat,IoTidal_results.topography/1000,[0.25 0.25],'linewidth',2,'Color','white');
clabel(C1,h1,'Color','white');
clabel(C2,h2,'Color','black');
clabel(C4,h4,'Color','white');
clabel(C5,h5,'Color','white');
clabel(C6,h6,'Color','white');
set(fig3,'Ytick',0:45:180,'Xtick',0:90:360,'Units','normalized','FontUnits','points','FontSize',26,'FontName','Times','ydir','normal');
title(fig3,'c) Topography (km)','position',[80 -1], 'FontUnits','points','interpreter','latex','Fontweight','normal','FontSize', 26,'FontName','Times');
xlabel(fig3,{'Longitude'},'FontUnits','points','interpreter','latex','FontSize', 30,'FontName','Times');
set(gca,'Yticklabel',[]);
fig3.YDir = 'reverse';
axis([0 360 0 180]);
colorbar
colormap(fig3,flipud(brewermap(128,'BuPu')));

num_col = 3;
num_row = 1;
AxesHandle=findobj(gcf,'Type','axes');
AxesHandle=flipud(AxesHandle);
Wstart_offset = 0.05;
Hstart_offset = 0.13;
Wend_offset = 0.07;
Hend_offset = 0.1;
W_separation = 0.06;
H_separation = 0.1;
plot_width  = (1-Wstart_offset-(num_col-1)*W_separation-Wend_offset)/num_col;
plot_height = (1-Hstart_offset-(num_row-1)*H_separation-Hend_offset)/num_row; % taken a bit off from num_row as want last thinner
for (ii=1:num_row)
    for (i=1:num_col)
        set(AxesHandle((ii-1)*num_col + i),'Position',[Wstart_offset+(i-1)*(plot_width+W_separation),Hstart_offset+(num_row-ii)*(plot_height+H_separation),plot_width,plot_height]);
    end
end

end