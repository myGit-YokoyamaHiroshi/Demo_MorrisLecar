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
% C    =  5;
% gL   =  2;
% gK   =  8;
% gCa  =  4;
% VL   = -60;
% VK   = -80;
% VCa  =  120;
% V1   = -1.2;
% V2   =  18;
% V3   =  12;
% V4   =  17.4;
% Iext =  -10:2:150;
% phi  =  1/15; %unit: 1/msec 

%%% typical parameter setting for Type II mode
C    =  5;
gL   =  2;
gK   =  8;
gCa  =  4.4;
VL   = -60;
VK   = -80;
VCa  =  120;
V1   = -1.2;
V2   =  18;
V3   =  2;
V4   =  30;
Iext =  0:5:240;
phi  =  1/25; %unit: 1/msec 

X0     = [0, 0]; % initial value of state variables
                 % X0(1): membrane potential, v
                 % X0(2): recovery variable,  w
%%%%% parameter settings
%%
color_list  = turbo(6);
eqpt_labels = {};
eqpt_idx    = [];
V_list      = [];
I_list      = [];

cnt         = 1;
h = waitbar(0,'running');
for i = 1:length(Iext)
    %% %%%%% Calculate equilibrium points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_init = -60:1:60;
    n_init = zeros(size(x_init));
    
    X_eq   = zeros(length(x_init), 2);
    
    %%%%% Calculate equilibrium points
    for k1 = 1:length(x_init)
        x0   = [x_init(k1), n_init(k1)];
        Xtmp = newton_method(x0, @MorrisLecar, 1E-16, 5000, C, gL, gK, gCa,...
                                                            VL, VK, VCa,...
                                                            V1, V2, V3, V4,...
                                                            Iext(i), phi);
        X_eq(k1,:) = Xtmp;
    end
    %%% Extract unique solution within torelance
    tmpval = uniquetol(X_eq(:,2), 0.05);
    v_eq   = zeros(1, length(tmpval));
    n_eq   = zeros(1, length(tmpval));
    for k2 = 1:length(tmpval)
        idx = find(X_eq(:,2)==tmpval(k2),1);
        
        %%%%% Get solution of equilibrium point
        v_eq(k2) = X_eq(idx,1);
        n_eq(k2) = X_eq(idx,2);
    end
    %%%%%%% Calculate equilibrium points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 
    for j = 1:length(v_eq)
        X = [v_eq(j), n_eq(j)];
        %%%%%%% Get jacobian matrix at equilibrium point [v_eq(i), n_eq(i)]
        J = jacobian_matrix_MorrisLecar(X, C, gL, gK, gCa,...
                                           VL, VK, VCa,...
                                           V1, V2, V3, V4,...
                                           Iext(i), phi);
        [eigvec, eigvalue] = eig(J);
        eigvalue           = diag(eigvalue);
        
        V_list(cnt) = v_eq(j);
        I_list(cnt) = Iext(i);
        %%%%%%% Determine its stability
        if all(imag(eigvalue)==0)
            if all(real(eigvalue)>0)
                labels        = 'Unstable node';
                eqpt_idx(cnt) = 1;
            elseif all(real(eigvalue)<0)
                labels        = 'Stable node';
                eqpt_idx(cnt) = 2;
            else
                labels        = 'Saddle node';
                eqpt_idx(cnt) = 3;
            end
        else
            if all(real(eigvalue)<0)
                labels        = 'Stable focus';
                eqpt_idx(cnt) = 4;
            elseif any(real(eigvalue)>0)
                labels        = 'Unstable focus';
                eqpt_idx(cnt) = 5;
            else
                labels        = 'Center (Hopf)';
                eqpt_idx(cnt) = 6;
            end
        end
        
        if cnt == 1
            eqpt_labels = {labels};
        else
            eqpt_labels = [eqpt_labels, {labels}];
        end
        cnt = cnt + 1;
    end

    %%%%%% Progress bar
    if mod(floor(i/length(Iext)*100), 1) == 0
        waitbar(i/length(Iext), h, ['Progress...', num2str(floor(i/length(Iext)*100)) , '%'])
    end
end
close(h)
%% Determine regions where the system has a limit cycle
periodic = zeros(size(V_list));
for i = 1:length(Iext)
    idx = find(I_list==Iext(i));
    tmp = zeros(size(idx));
    for j = 1:length(idx)
        if contains(eqpt_labels{idx(j)}, 'Stable') %&& length(idx)~=1
            tmp(j) = 1;
        end
    end

    if sum(tmp)==0
        periodic(i) = 1;
    end
end

Iperi   = Iext(periodic==1);
Vmaxmin = zeros(length(Iperi), 2);
for i = 1:length(Iperi)
    X      = zeros(Nt, length(X0));
    X(1,:) = X0;
    
    for j = 2:Nt
        X_now  = X(j-1,:);
        %%%%% Numerical integral scheme with 4th order Runge Kutta method
        X(j,:) = runge_kutta(X_now, dt, @MorrisLecar, C, gL, gK, gCa,...
                                                      VL, VK, VCa,...
                                                      V1, V2, V3, V4,...
                                                      Iperi(i), phi);
    end
    
    Vmaxmin(i,:) = [min(X(:,1)), max(X(:,1))];
end

%%
fig = figure(1);
figure_setting(40, 20, fig)
%%%%%% Show the region where the system has a limit cycle
plot(Iperi, Vmaxmin, 'b', 'linewidth', 4,'HandleVisibility','off')
hold on
X  =[Iperi, fliplr(Iperi)];
Y = [Vmaxmin(:,1).', fliplr(Vmaxmin(:,2).')];
h = fill(X, Y, 'k','DisplayName','limit cycle');

set(h,'facealpha',.1)
%%%%%% Show the I-V curve 
for i = 1:6
    if sum(eqpt_idx==i) ~=0
        idx = find(eqpt_idx==i);
        scatter(I_list(idx), V_list(idx), ...
            'MarkerEdgeColor', color_list(i,:),...
            'MarkerFaceColor', color_list(i,:),...
            'DisplayName', eqpt_labels{idx(1)})

        legend('-DynamicLegend');
    end
end
xlabel('parameter \it I')
ylabel('Membrane potential V_*')
legend('location', 'southeastoutside')

fname = [filepath, filesep, 'figures', filesep, 'ad_ex2', filesep, 'bifurcation'];
figure_save(fig, fname)
%%
