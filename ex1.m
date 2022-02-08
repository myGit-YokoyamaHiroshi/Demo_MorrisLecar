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
% Iext =  100.5;
% phi  =  1/25; %unit: 1/msec 

X0     = [0, 0]; % initial value of state variables
                 % X0(1): membrane potential, v
                 % X0(2): recovery variable,  w
%%%%% parameter settings
%% Solve differential equation
X      = zeros(Nt, length(X0));
X(1,:) = X0;

for i = 2:Nt
    X_now  = X(i-1,:);
    %%%%% Numerical integral scheme with 4th order Runge Kutta method
    X(i,:) = runge_kutta(X_now, dt, @MorrisLecar, ...
                                    C, gL, gK, gCa,...
                                       VL, VK, VCa,...
                                       V1, V2, V3, V4,...
                                       Iext, phi);
end
%%
fig = figure(1);
figure_setting(60, 40, fig);

sfh1 = subplot(2,1,1,'parent', fig);
plot(time, X(:,1), 'LineWidth', 3);
hold on
plot(time, X(:,2), 'LineWidth', 3);
hold off

xlabel('time (ms)')
ylabel('V, N')
lgnd = legend({'membrane potential \it V', 'recovery variable \it N'}, 'location', 'northeastoutside');
%%%%%%%
sfh2 = subplot(2,1,2,'parent', fig);
plot(X(:,1), X(:,2), 'r', 'LineWidth', 3);
xlabel('membrane potential \it V')
ylabel('recovery variable \it N')
title('phase space')
axis square
sfh2.Position = sfh2.Position - [0.1, 0, 0, 0];

fname = [filepath, filesep, 'figures', filesep, 'ex1', filesep, 'result'];
figure_save(fig, fname)