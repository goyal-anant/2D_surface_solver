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
N         = 50;
radius    = 0.3*lambda;
circumfer = 2*pi*radius;
l         = circumfer/N; %length of each segment
test_pts  = linspace(0,)




%% Formulating the problem


%% Plotting the result