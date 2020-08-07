function S = IoTidal(lam_M,lam_q,Pc,Psi,Psi_tot,plot_stuff)

% Box model for Io,
% lam_M = constant emplacement rate, units s^-1
% lam_q = emplacement proportionality to qp, no unit
% Pc = Critical overpressure decompacting boundary layer
% Psi is a tidal heating distribution. Can give a vector or an integer
% Psi_tot = total heating rate for non-dimensionalisation, can leave as nan for 1e14 W

% Dimensional Parameters
K_0 = 10^-7; % reference permeability (m^2)
rho = 3000; % density (kg/m^3)
del_rho = 500; % density difference (kg/m^3)
g = 1.5; % gravity (m/s^2)
L = 4e5; % latent heat (J/kg)
kappa = 1e-6; % Thermal diffusivity (m^2/s)
c = 1200; % specific heat (J/kg/K)
T_l = 1350; % relative melting point (above surface T) (K)
n = 3; % permeability exponent
eta_l = 1; % basalt melt viscosity (Pas)
eta = 1e20; % viscosity (Pas)

r_s = 1820e3; % Io radius (m)
r_m = 700e3/1820e3; % normalise core radius

if Psi_tot ~= 0 && ~isnan(Psi_tot)
    Psi_ref = Psi_tot/(4/3 * pi *(r_s^3 - 700e3^3)); % reference tidal heating (W/m^3)
else
    Psi_ref = 1e14/(4/3 * pi *(r_s^3 - 700e3^3)); % reference tidal heating (W/m^3)
end
q_0 = Psi_ref*r_s/rho/L; % reference velocity (m/s)
phi_0 = (q_0*eta_l/(K_0*del_rho*g))^(1/n); % reference porosity

zeta_0 = eta/phi_0;
P_0 = zeta_0*q_0/r_s;
S.P_c = -Pc*P_0;

% Dimensionless Parameters
St = L/c/T_l; % Stefan number
Pe = q_0*r_s/kappa; % Peclet number

% Makes tidal heating an analytical function.
% If a vector of heating is input this fits a function to that data.
% If a scalar is input, tidal heating is that scalar
if length(Psi) > 1
    Psi = Psi/Psi_ref; % normalise input heating to reference
    % if a vector is input, get an analytical function for heating
    x = linspace(r_m,1,length(Psi))';
    psi_func = fit(x,Psi,'linearinterp');
else
    % probably redundant, but makes psi a function
    psi_func = @(x) Psi;
end

% Use shooting method to find lithosphere thickness, see functions for details
opts = odeset('reltol',1e-6);
[l,a,b,iter] = bisect(1e-4,1-r_m-1e-4,1e-6);
% Now with known r_l solve ODEs on correct domain
xint = linspace(r_m,1-l)';
qrc = trapz(xint,xint.^2.*psi_func(xint))/(1-l)^2;
[x,y] = ode15s(@odes1,[1-l 1],[qrc; 1; 0],opts);

if exist('p') == 0
    r = linspace(r_m,1,1000)'; % create position vector
else r = linspace(r_m,1,p.ni)'; % create position vector
end
dr = r(2)-r(1);
r1 = r(r>(1-l)); % Portion of position vector corresponding to lithosphere
r2 = r(r<(1-l)); % Portion of position vector corresponding to mantle

T = 1.*(r<=(1-l)); % Mantle is on the liquidus
T(r>(1-l)) = interp1(x,y(:,2),r1); % Interpolate the ODE temp solution onto lithosphere grid

qe = 0*r.^0; % qe is zero in mantle
qe(r>(1-l)) = interp1(x,y(:,1),r1); % Interpolate the ODE qe solution onto lithosphere grid

dT(r>(1-l)) = interp1(x,y(:,3),r1);

q = zeros(length(r),1);
phi = zeros(length(r),1);
Pbl = zeros(length(r),1);
P = zeros(length(r),1);

rl = 1-l; % radial position of lithosphere boundary

% ****************************************** %
% Decompacting BL
% [phi_bl P_bl Z] = DBL_tidal(Pc,l,psi_func);
% 
% delta = eta*q_0/(phi_0*del_rho*g*r_s^2);
% q(r<rl) = cumtrapz(r(r<rl),r(r<rl).^2.*psi_func(r(r<rl)))./(r(r<rl)).^2;
% phi(r<rl) = q(r<rl).^(1/n) + interp1(Z,phi_bl,(rl-r2)/delta) - qrc^(1/n); % use boundary layer formulation for q
% phi(1) = 0;
% u = -q-qe; % u from conservation of mass
% 
% P(r<rl) = -psi_func(r(r<rl))./(q(r<rl).^(1/n)) + interp1(Z,P_bl,(rl-r2)/delta) + psi_func(rl)/qrc^(1/n); % outer solution + inner solution - P_inf (from BL calc)
% P(isinf(P)) = 0;
% ****************************************** %

% ****************************************** %
% No decompacting BL
q(r<rl) = cumtrapz(r(r<rl),r(r<rl).^2.*psi_func(r(r<rl)))./(r(r<rl)).^2;
phi(r<rl) = q(r<rl).^(1/n);
P(r<rl) = -psi_func(r(r<rl))./(q(r<rl).^(1/n));
u = -q-qe; % u from conservation of mass
% ****************************************** %

% Total heating
Htot = trapz(r*r_s,4*pi*(r*r_s).^2.*psi_func(r)*Psi_ref);

% Store in outputting structure
S.r = r*r_s;
S.T = T*T_l + 150;
S.dT = dT;
S.qe = qe*q_0;
S.q = q*q_0;
S.u = u*q_0;
S.phi = phi*phi_0;
S.P = P*P_0;
S.Psi_tot = Htot;
S.Psi = Psi*Psi_ref;

P = P*P_0/1e6;
T = T*T_l + 150; % Put temperature into dimensional units and shift to absolute T
u = u *Psi_ref*r_s/(rho*L); % dimensionalise velocity
q = q *Psi_ref*r_s/(rho*L); 
qe = qe *Psi_ref*r_s/(rho*L);
r = r * r_s; % dimsonsionalise radial position
rl = rl * r_s;
S.rl = (1820e3-rl)/1000;

if lam_q == 0
    M_dim = lam_M*q_0/r_s; % get dimensional emplacement rate
    qe_lin = Htot/(4*pi*r_s^2*(rho*L + rho*c*T_l)); % eruption rate from global energy balance
    rc_lin = (3/M_dim *qe_lin*r_s^2 + r_s^3 - 3*Htot/(4*pi*M_dim*rho*L))^(1/3);
    qp_lin = qe(end)*r_s^2./r.^2 + M_dim/3 *(r_s^3./r.^2 -r);
    qp_lin = qp_lin*(100*60*60*24*365);
else
    lam_dim = lam_q/r_s; % get dimensional lam_q
    qe_lin = Htot/(4*pi*r_s^2*(rho*L + rho*c*T_l)); % eruption rate from global energy balance
    rc_lin = r_s - 1/lam_dim *log((L+c*T_l)/L);
    qp_lin = qe_lin*r_s^2./r.^2 .*exp(lam_dim*(r_s-r));
    qp_lin = qp_lin*(100*60*60*24*365);
end

if plot_stuff ~= 0
    ss_plot(T,phi,q*(100*60*60*24*365),qe*(100*60*60*24*365),u*(100*60*60*24*365),P,S.Psi,r);
    
    erupt_rate = u(end) *(100*60*60*24*365) % erupt rate in cm/yr
    cond_heat = -kappa*rho*c* dT(end)*T_l/r_s * (4*pi*r_s^2) % global conductive heat flux (W) kappa*rho*c = conductivity
    volc_heat = rho*qe(end) *(4*pi*r_s^2)*(c*T_l + L)
end

function [m,a,b,iter] = bisect(a,b,tol)
    % Looking for the value of lithosphere thickness that gives the correct surface temperature
    % Try two values and look for sign change, if there is one, narrow the range, if not then try another inteval.
    % Two values to try are 'a' and 'b'
    % When the difference between quesses is within a tollerance, success :-)
    iter = 0;
    fa = shoot(a); fb = shoot(b); % return surface temperatures and place them (hopefully) either side of zero.
    if fa*fb > 0
        %warning('bisect: Signs at ends of interval are the same');
        % Try looking for another interval
        for tmp = 0.1:.1:1
            c = a+tmp*(b-a);
            fc = shoot(c); 
            %if fa*fc < 0, b = c; fb = fc; warning(['bisect: Using reduced interval ',num2str([a b])]); break; end
            if fa*fc < 0, b = c; fb = fc; break; end
        end
        % Otherwise give up
        if tmp==1
            warning('bisect: Could not find an interval with a sign change');
            m = NaN;
            return;
        end
    end
    % Goes here when two guesses change sign (i.e. sit either side of surface temperature)
    while b-a>tol
        % Try halfway through the current interval, if it changes sign that's new interval, if not the leftover must be the inteval.
        iter = iter+1;
        m = (a+b)/2;
        fm = shoot(m);
        if fm*fa<0 
           b = m;
        else
           a = m;
           fa = fm;
        end
    end
end

function out = shoot(l)
    % solve ODEs for a guess of the lithosphere thickness, return the predicted surface temperature
    % Guesses for qe and dT take account of whether qe is zero or not. If it's zero in dT case, use q_e BC with q_e = 0.
    xint = linspace(r_m,1-l)';
    qrc = trapz(xint,xint.^2.*psi_func(xint))/(1-l)^2;
    [x1,y1] = ode15s(@odes1,[1-l 1],[qrc; 1; 0],opts);
    out = y1(end,2);
end

function dydx = odes1(x,y)
    % stiff ODE solver for the 3 coupled 1st order equations
    qe = y(1,:);
    T = y(2,:);
    dT = y(3,:);
    M = (lam_M + lam_q.*qe).*(qe>0);
    dydx(1,:) = - M -2*qe/x; % M = 0 if qe<0
    dydx(2,:) = dT;
    dydx(3,:) = -Pe*qe*dT - Pe*St*psi_func(x) - 2*dT/x - Pe*M.*(St + (1-T)).*(qe>0); % Energy equation rewritten for d2T/dr2. M = 0 if qe<0
end

function lines = ss_plot(T,phi,q,qe,u,P,Psi,r)
    for ii=2:length(phi)
        if phi(ii)<=1e-5 && phi(ii-1)>1e-5
            lith = r(ii);
        end
    end
    
    figure('Units','centimeters','Position',[15 15 50 25],'PaperPositionMode','auto');
    fig1 = subplot(1,5,1);
    patch(fig1,[-10 -10 10 10],[lith/1000 r_s/1000 r_s/1000 lith/1000],[0.852 0.852 0.852])
    hold on
    patch(fig1,[-10 -10 10 10],[600 700 700 600],[0.632 0.5 0.5])
    gradient_patch = zeros(1,4,3);
    gradient_patch(1,1,:) = [0.696 0.696 0.696];
    gradient_patch(1,2,:) = [0.9 0.516 0.168];
    gradient_patch(1,3,:) = [0.9 0.516 0.168];
    gradient_patch(1,4,:) = [0.696 0.696 0.696];
    patch(fig1,[-10 -10 10 10],[700 lith/1000 lith/1000 700],gradient_patch)
    line = plot(fig1,q+qe,r/1000,'linewidth',2);
    %line = plot(fig1,[-10 10],[1820-elas_thick 1820-elas_thick],'k--','linewidth',1.5);
    line = plot(fig1,u,r/1000,'linewidth',2,'color',[0 0.416 0.1]);
    line = plot(fig1,[-10 10],[700 700]);
    line = plot(fig1,[0 0], [700 1820],':','linewidth',2);
    plot(fig1,[-10 10],[rc_lin/1000 rc_lin/1000],'--k','linewidth',1.5);
    %plot(fig1,qp_lin,r/1000,'linewidth',1.5);
    
    
    axis([-10 10 600 1820]);
    set(fig1,'Units','normalized','FontUnits','points','FontSize',24,'FontName','Times', 'Layer', 'top')
    ylabel(fig1,{'Radial position $r$ (km)'},'FontUnits','points', 'interpreter','latex','FontSize', 24,'FontName','Times')
    xlabel(fig1,{'Upwelling flux (cm/yr)'},'FontUnits','points','interpreter','latex','FontSize', 24,'FontName','Times')
    title(fig1,'a)','position',[-10 1825], 'FontUnits','points','Fontweight','normal','FontSize',24,'FontName','Times')
    %legend(fig3,'Liquid','Solid','location','southwest','FontSize',14)
    
    fig2 = subplot(1,5,2);
    patch(fig2,[0 0 2000 2000],[lith/1000 r_s/1000 r_s/1000 lith/1000],[0.852 0.852 0.852])
    hold on
    patch(fig2,[0 0 2000 2000],[600 700 700 600],[0.632 0.5 0.5])
    patch(fig2,[0 0 2000 2000],[700 lith/1000 lith/1000 700],gradient_patch)
    line = plot(fig2,T,r/1000,'linewidth',2);
    line = plot(fig2,[0 2000],[700 700]);
    %line = plot(fig2,[0 2000],[1820-elas_thick 1820-elas_thick],'k--','linewidth',1.5);
    axis([0 2000 600 1820]);
    set(fig2,'Units','normalized','FontUnits','points','FontSize',24,'FontName','Times', 'Layer', 'top')
    xlabel(fig2,{'Temperature (K)'},'FontUnits','points','interpreter','latex','FontSize', 24,'FontName','Times')
    title(fig2,'b)','position',[0 1825],'FontUnits','points','Fontweight','normal','FontSize',24,'FontName','Times')
    set(gca,'Yticklabel',[])

    fig3 = subplot(1,5,3);
    patch(fig3,[0 0 20 20],[lith/1000 r_s/1000 r_s/1000 lith/1000],[0.852 0.852 0.852])
    hold on
    patch(fig3,[0 0 20 20],[600 700 700 600],[0.632 0.5 0.5])
    patch(fig3,[0 0 20 20],[700 lith/1000 lith/1000 700],gradient_patch)
    line = plot(fig3,phi*phi_0*100,r/1000,'linewidth',2);
    %line = plot(fig3,[0 20],[1820-elas_thick 1820-elas_thick],'k--','linewidth',1.5);
    line = plot(fig3,[0 5],[700 700]);
    axis([0 5 600 1820]);
    set(fig3,'Units','normalized','XTick',0:5:20,'FontUnits','points','FontSize',24,'FontName','Times', 'Layer', 'top')
    xlabel(fig3,{'Porosity (\%)'},'FontUnits','points','interpreter','latex','FontSize', 24,'FontName','Times')
    title(fig3,'c)','position',[0 1825], 'FontUnits','points','Fontweight','normal','FontSize',24,'FontName','Times')
    set(gca,'Yticklabel',[])
    
    fig4 = subplot(1,5,4);
    patch(fig4,[-150 -150 30 30],[lith/1000 r_s/1000 r_s/1000 lith/1000],[0.852 0.852 0.852])
    hold on
    patch(fig4,[-150 -150 30 30],[600 700 700 600],[0.632 0.5 0.5])
    patch(fig4,[-150 -150 30 30],[700 lith/1000 lith/1000 700],gradient_patch)
    line = plot(fig4,[0 0], [700 1820],'-k');
    line = plot(fig4,P(2:end),r(2:end)/1000,'linewidth',2);
    %line = plot(fig4,[-150 30],[1820-elas_thick 1820-elas_thick],'k--','linewidth',1.5);
    line = plot(fig4,[-150 30],[700 700]);
    axis([-150 30 600 1820]);
    set(fig4,'Units','normalized','XTick',-150:50:0,'FontUnits','points','FontSize',24,'FontName','Times', 'Layer', 'top')
    xlabel(fig4,{'Pressure (MPa)'},'FontUnits','points','interpreter','latex','FontSize', 24,'FontName','Times')
    title(fig4,'d)','position',[-150 1825], 'FontUnits','points','Fontweight','normal','FontSize',24,'FontName','Times')
    set(gca,'Yticklabel',[])
    
    fig5 = subplot(1,5,5);
    patch(fig5,[0 0 20 20],[lith/1000 r_s/1000 r_s/1000 lith/1000],[0.852 0.852 0.852])
    hold on
    patch(fig5,[0 0 20 20],[600 700 700 600],[0.632 0.5 0.5])
    patch(fig5,[0 0 20 20],[700 lith/1000 lith/1000 700],gradient_patch)
    line = plot(fig5,Psi*(1000^3)/1000,r/1000,'linewidth',2);
    line = plot(fig5,[0 20],[700 700]);
    %line = plot(fig5,[0 20],[1820-elas_thick 1820-elas_thick],'k--','linewidth',1.5);
    axis([0 20 600 1820]);
    set(fig5,'Units','normalized','FontUnits','points','FontSize',24,'FontName','Times', 'Layer', 'top')
    xlabel(fig5,{'Tidal dissipation (kW/km$^3$)'},'FontUnits','points','interpreter','latex','FontSize', 24,'FontName','Times')
    title(fig5,'e)','position',[0 1825],'FontUnits','points','Fontweight','normal','FontSize',24,'FontName','Times')
    set(gca,'Yticklabel',[])
    
    linkaxes([fig1 fig2 fig3 fig4 fig5],'y');
    
    num_col = 5;
    num_row = 1;
    AxesHandle=findobj(gcf,'Type','axes');
    AxesHandle=flipud(AxesHandle);
    Wstart_offset = 0.07;
    Hstart_offset = 0.1;
    Wend_offset = 0.03;
    Hend_offset = 0.04;
    W_separation = 0.04;
    H_separation = 0.12;
    plot_width  = (1-Wstart_offset-(num_col-1)*W_separation-Wend_offset)/num_col;
    plot_height = (1-Hstart_offset-(num_row-1)*H_separation-Hend_offset)/num_row; % taken a bit off from num_row as want last thinner
    for (ii=1:num_row)
        for (i=1:num_col)
            set(AxesHandle((ii-1)*num_col + i),'Position',[Wstart_offset+(i-1)*(plot_width+W_separation),Hstart_offset+(num_row-ii)*(plot_height+H_separation),plot_width,plot_height]);
        end
    end
    
end

end