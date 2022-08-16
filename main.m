clc; clear all; close all;

%% Defining constants
lambda    = 1;
epsilonr  = 4;
k1        = 2*pi/lambda;
k2        = k1*sqrt(epsilonr);
phi0      = 0; %incident angle
theta0    = pi/2;

%% Building the scatterer 
% plot_flag = 1;
%[X, Y] = structure(N,lambda,radius,plot_flag);
N         = 10;
radius    = 0.3*lambda;
l         = 2*pi*radius/N; %length of each segment
step      = l;
test_pts  = 0:step:2*pi-step/2;
scatter(radius*cos(test_pts),radius*sin(test_pts),'k','filled');
grid on; axis('equal');




%% Formulating the problem


%% Plotting the result
