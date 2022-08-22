clc; clear all; close all;

%% Defining constants
lambda    = 1;
epsilonr  = 4;
k1        = 2*pi/lambda;
k2        = k1*sqrt(epsilonr);
phi0      = 0; %incident angle
theta0    = pi/2;
E0        = 1; %incident field amplitude
gamma     = 0.577; %for hankel fun approximation

%% Building the scatterer 
% plot_flag = 1;
%[X, Y] = structure(N,lambda,radius,plot_flag);
N          = 50; %number of segments on the surface(boundary)
radius     = 0.3*lambda;
circum     = 2*pi*radius;
l          = circum/N; %length of each segment
step       = 2*pi/N; %angle measure of each segment
nodes      = 0:step:2*pi-step; %nodes of each segment
test_pts   = step/2 : step : 2*pi-step/2; %testing points which are choosen as the center of each segment
% scatter(radius*cos(nodes),radius*sin(nodes),'k','filled');
% hold on; grid on; axis('equal');pause(2)
% scatter(radius*cos(test_pts),radius*sin(test_pts),'r','filled');
% grid on; axis('equal');


%% Formulating the problem
%defining the known side of Ax = c i.e. c
%for phi_i(or Ei) vector
X1 = radius*cos(test_pts);
Y1 = radius*sin(test_pts);
alpha = (X1 * sin(theta0) * cos(phi0)) + (Y1 * sin(theta0) * sin(phi0));
Ei = E0 * exp(-1i * k1 * alpha);     %incident electric field 
%defining 'c' vector
c = [Ei zeros(1,N)];
%creating the 'A' matrix
gdiag = @(k) -1j/4*(l - k^2*l^3/48 - 1j*(2*l/pi*(log(l*k/(4*exp(1)))+gamma)));
for i = 1:N
    rp = [radius*cos(test_pts(i)) radius*sin(test_pts(i))];
    for j = 1:N
        r    = [radius*cos(test_pts(j)) radius*sin(test_pts(j))];
        nhat = [cos(test_pts(j)) sin(test_pts(j))];
        if i == j
            A(i,j)     = gdiag(k1);
            A(i,j+N)   = 1/2;

            A(i+N,j)   = gdiag(k2);
            A(i+N,j+N) = -1/2;
        else
            A(i,j)     = integral(@(r) green(k1)) 
        end
    end
end

%% Plotting the result

%% Functions

%% green function
function g = green(k,rho)
    g = (-1j/4)*besselh(0,2,k*rho);
end
%% grad_green
function grad_g = gradg(k,rho)
    grad_g = (1j*k/4)*besselh(1,2,k*rho);
end