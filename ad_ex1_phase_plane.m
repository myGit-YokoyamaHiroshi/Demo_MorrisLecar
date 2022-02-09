clear 
close all
clc
filepath = fileparts(mfilename('fullpath'));
addpath(genpath(filepath));
%% code font settings
%%%% Set "Arial" as the Default font
set(0,'defaultAxesFontSize',16);
set(0,'defaultAxesFontName','Arial');
set(0,'defaultTextFontSize',16);
set(0,'defaultTextFontName','Arial');

set(0,'defaultUipanelFontName','Arial');
set(0,'defaultUicontrolFontName','Arial');
%%
Nt     = 50000;  % Num. of sample
dt     = 0.01;   % time step for numerical integration; unit : msec
time   = linspace(0, Nt-1, Nt) * dt; % time vector; unit : msec
%%%%% parameter settings
%%% typical parameter setting for Type I mode
C    =  5;
gL   =  2;
gK   =  8;
gCa  =  4;
VL   = -60;
VK   = -80;
VCa  =  120;
V1   = -1.2;
V2   =  18;
V3   =  12;
V4   =  17.4;
Iext =  39.8;
phi  =  1/15; %unit: 1/msec 

%%% typical parameter setting for Type II mode
% C    =  5;
% gL   =  2;
% gK   =  8;
% gCa  =  4.4;
% VL   = -60;
% VK   = -80;
% VCa  =  120;
% V1   = -1.2;
% V2   =  18;
% V3   =  2;
% V4   =  30;
% Iext =  82;
% phi  =  1/25; %unit: 1/msec 

X0     = [0, 0]; % initial value of state variables
                 % X0(1): membrane potential, v
                 % X0(2): recovery variable,  w
%%%%% parameter settings
%% Solve differential equation
X_trj      = zeros(Nt, length(X0));
X_trj(1,:) = X0;

for i = 2:Nt
    X_now  = X_trj(i-1,:);
    %%%%% Numerical integral scheme with 4th order Runge Kutta method
    X_trj(i,:) = runge_kutta(X_now, dt, @MorrisLecar, ...
                                         C, gL, gK, gCa,...
                                         VL, VK, VCa,...
                                         V1, V2, V3, V4,...
                                         Iext, phi);
end
%% Get Nullcline
[V_null, N_null] = get_nullcline_MorrisLecar(gL, gK, gCa,...
                                             VL, VK, VCa,...
                                             V1, V2, V3, V4,...
                                             Iext, -65, 65);
%% Determine the stability of each equilibrium points

x_init = -60:1:60;
n_init = zeros(size(x_init));

X_eq   = zeros(length(x_init), 2);

%%%%% Calculate equilibrium points
for i = 1:length(x_init)
    x0   = [x_init(i), n_init(i)];
    Xtmp = newton_method(x0, @MorrisLecar, 1E-16, 5000, C, gL, gK, gCa,...
                                                        VL, VK, VCa,...
                                                        V1, V2, V3, V4,...
                                                        Iext, phi);
    X_eq(i,:) = Xtmp;
end

%%% Extract unique solution within torelance
tmpval = uniquetol(X_eq(:,2), 0.05);
v_eq   = zeros(1, length(tmpval));
n_eq   = zeros(1, length(tmpval));
for i = 1:length(tmpval)
    idx = find(X_eq(:,2)==tmpval(i),1);
    
    %%%%% Get solution of equilibrium point
    v_eq(i) = X_eq(idx,1);
    n_eq(i) = X_eq(idx,2);
end

eqpt_labels  = cell(1, length(v_eq));
for i = 1:length(v_eq)
    X = [v_eq(i), n_eq(i)];
    %%%%%%% Get jacobian matrix at equilibrium point [v_eq(i), w_eq(i)]
    J = jacobian_matrix_MorrisLecar(X, C, gL, gK, gCa, ...
                                       VL, VK, VCa,...
                                       V1, V2, V3, V4,...
                                       Iext, phi);
    [eigvec, eigvalue] = eig(J);
    eigvalue = diag(eigvalue);
    
    %%%%%%% Determine its stability
    if all(imag(eigvalue)==0)
        if all(real(eigvalue)>0)
            eqpt_labels{i} = 'Unstable node';
        elseif all(real(eigvalue)<0)
            eqpt_labels{i} = 'Stable node';
        else
            eqpt_labels{i} = 'Saddle node';
        end
    else
        if all(real(eigvalue)<0)
            eqpt_labels{i} = 'Stable focus';
        elseif any(real(eigvalue)>0)
            eqpt_labels{i} = 'Unstable focus';
        else
            eqpt_labels{i} = 'Center (Hopf)';
        end
    end
end
%% plot results
fig = figure(1);
figure_setting(60, 60, fig)

sfh1 = subplot(2,1,1,'parent', fig);
plot(time, X_trj(:,1), 'LineWidth', 3);

xlabel('time (ms)')
ylabel('membrane potential \it v')

%%% plot trajectory
sfh2 = subplot(2,1,2,'parent', fig);
plot(X_trj(:,1), X_trj(:,2), 'r', 'LineWidth', 3);
trj_labels = {'trajectory'};

hold on

%%% plot nullcline
plot(V_null, N_null, 'k-', 'LineWidth', 1,'HandleVisibility','off'); 

for i = 1:length(v_eq)
    %%% plot equilibrium point
    h_eq(i) = scatter(v_eq(i), n_eq(i), 60,'filled');
end
xlabel('membrane potential \it V')
ylabel('recovery variable \it N')
legend([trj_labels, eqpt_labels], 'location', 'northeastoutside')


title('phase plane')
xlim([-60, 60])
ylim([-0.2, 1.2])

alpha(0.8)
sfh2.Position = sfh2.Position - [0.05, 0, 0, 0];
axis square
fname = [filepath, filesep, 'figures', filesep, 'ad_ex1', filesep, 'stability'];
figure_save(fig, fname)
hold off
