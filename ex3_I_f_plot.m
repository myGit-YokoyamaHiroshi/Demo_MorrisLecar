clear 
close all
clc
filepath = fileparts(mfilename('fullpath'));
addpath(genpath(filepath));
%% code font settings
%%%% Set "Arial" as the Default font
set(0,'defaultAxesFontSize',14);
set(0,'defaultAxesFontName','Arial');
set(0,'defaultTextFontSize',14);
set(0,'defaultTextFontName','Arial');

set(0,'defaultUipanelFontName','Arial');
set(0,'defaultUicontrolFontName','Arial');
%%
Nt     = 50000;  % Num. of sample
dt     = 0.01;   % time step for numerical integration; unit : msec
time   = linspace(0, Nt-1, Nt) * dt; % time vector; unit : msec

mode   = 'type1';

%%%%% parameter settings
if strcmp(mode, 'type1')
    %%% typical parameter setting for Type I mode
    C      =  5;
    gL     =  2;
    gK     =  8;
    gCa    =  4;
    VL     = -60;
    VK     = -80;
    VCa    =  120;
    V1     = -1.2;
    V2     =  18;
    V3     =  12;
    V4     =  17.4;
    I_list =  30:5:120;
    phi    =  1/15; %unit: 1/msec 

    I_plt  = [30, 40, 50, 60, 70, 80];
elseif strcmp(mode, 'type2')    
    %%% typical parameter setting for Type II mode
    C      =  5;
    gL     =  2;
    gK     =  8;
    gCa    =  4.4;
    VL     = -60;
    VK     = -80;
    VCa    =  120;
    V1     = -1.2;
    V2     =  18;
    V3     =  2;
    V4     =  30;
    I_list =  70:5:160;
    phi    =  1/25; %unit: 1/msec 

    I_plt = [80, 90, 100, 110, 120, 130];
end

X0     = [0, 0]; % initial value of state variables
                 % X0(1): membrane potential, v
                 % X0(2): recovery variable,  w
%%%%% parameter settings
%% Solve differential equation


freqs  = zeros(size(I_list));
X      = zeros(Nt, length(X0), length(I_list));

h = waitbar(0,'running');
for j = 1:length(I_list)
    X(1,:,j) = X0;
    Iext     = I_list(j);
    for i = 2:Nt
        X_now  = X(i-1,:,j);
        %%%%% Numerical integral scheme with 4th order Runge Kutta method
        X(i,:,j) = runge_kutta(X_now, dt, @MorrisLecar, ...
                                          C, gL, gK, gCa,...
                                             VL, VK, VCa,...
                                             V1, V2, V3, V4,...
                                             Iext, phi);
    end
    
    pks      = findpeaks(X(:,1,j));
    Npks     = sum(pks>20);
    T        = (Nt * dt) * 10^-3;
    freqs(j) = Npks/T; 
    
    %%%%%% Progress bar
    if mod(floor(j/length(I_list)*100), 5) == 0
        waitbar(j/length(I_list), h, ['Progress...', num2str(floor(j/length(I_list)*100)) , '%'])
    end
end

close(h)
%%
fig = figure(1);
figure_setting(60, 70, fig);

sfh1 = subplot(2,4,2:3);
plot(I_list, freqs, '-o', 'linewidth', 1);
xlabel('\it I_{ext}')
ylabel('frequency (Hz)')
sfh1.Position = sfh1.Position + [0, 0, 0, -0.05];

xlim([I_list(1), I_list(end)])
ylim([0, 45])

for i = 1:length(I_plt)
    sfh2 = subplot(4,3,6+i);
    plot(time, X(:,1, I_list == I_plt(i)), 'LineWidth', 2);
    xlabel('time (ms)')
    ylabel(' \it V')
    xticks(0:200:Nt*dt)
    xlim([0, Nt*dt])
    title(['\it I_{ext} = ', num2str(I_plt(i))])

    sfh2.Position = sfh2.Position + [0, 0, 0, -0.02];
end

fname = [filepath, filesep, 'figures', filesep, 'ex3', filesep, 'result'];
figure_save(fig, fname)
