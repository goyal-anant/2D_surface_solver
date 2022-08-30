clear;

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
Ei = E0 * 1.25 * exp(-1i * k1 * alpha);     %incident electric field 

%defining 'c' vector
c = [Ei zeros(1,N)];

%creating the 'A' matrix, i.e., A = [U V; W X]
A = zeros(2*N,2*N);
gdiag = @(k) -1j/4*(l - k^2*l^3/48 - 1j*(2*l/pi*(log(l*k/(4*exp(1)))+gamma)));
% syms m
for i = 1:N
    rpnl = [radius*cos(test_pts(i) - (step/2)) radius*sin(test_pts(i) - (step/2))];
    rpnu = [radius*cos(test_pts(i) + (step/2)) radius*sin(test_pts(i) + (step/2))];
    for j = 1:N
        rnl = [radius*cos(test_pts(j) - (step/2)) radius*sin(test_pts(j) - (step/2))];
        rnu = [radius*cos(test_pts(j) + (step/2)) radius*sin(test_pts(j) + (step/2))];
        nhat   = [cos(test_pts(j)) sin(test_pts(j))];
        if i == j
            A(i,j)     = gdiag(k1);
            A(i,j+N)   = 1/2;
            A(i+N,j)   = gdiag(k2);
            A(i+N,j+N) = -1/2;
        else
            g          = @(m) green(rpnl,rpnu,0.5,rnl,rnu,m,k1);
            A(i,j)     = l * integral(g,0,1,'ArrayValued',true);

            gg         = @(m) gradg(rpnl,rpnu,0.5,rnl,rnu,m,nhat,k1);
            A(i,j+N)   = -l * integral(gg,0,1,'ArrayValued',true);

            g          = @(m) green(rpnl,rpnu,0.5,rnl,rnu,m,k2);
            A(i+N,j)   = l * integral(g,0,1,'ArrayValued',true);

            gg         = @(m) gradg(rpnl,rpnu,0.5,rnl,rnu,m,nhat,k2);
            A(i+N,j+N) = -l * integral(gg,0,1,'ArrayValued',true);
        end
    end
end

x = A\c';

%finding the total field now using huygen principle
oradius   = 200 * lambda;
oN        = 100;
ostep     = 2*pi/oN;
onodes    = 0:ostep:2*pi-ostep;
otest_pts = ostep/2 : ostep : 2*pi-ostep/2;
farfield  = zeros(1,oN); %scattered far field

for i = 1:oN
    rpnl = [oradius*cos(otest_pts(i) - (ostep/2)) oradius*sin(otest_pts(i) - (ostep/2))];
    rpnu = [oradius*cos(otest_pts(i) + (ostep/2)) oradius*sin(otest_pts(i) + (ostep/2))];
    for j = 1:N
        rnl  = [radius*cos(test_pts(j) - (step/2)) radius*sin(test_pts(j) - (step/2))];
        rnu  = [radius*cos(test_pts(j) + (step/2)) radius*sin(test_pts(j) + (step/2))];
        nhat = [cos(test_pts(j)) sin(test_pts(j))];


        g  = @(m) green(rpnl,rpnu,0.5,rnl,rnu,m,k1);
        gg = @(m) gradg(rpnl,rpnu,0.5,rnl,rnu,m,nhat,k1);

        farfield(i) = farfield(i) - l * (x(j) * integral(g,0,1,"ArrayValued",true) - ...
            x(j+N) * integral(gg,0,1,"ArrayValued",true));
    end
%     farfield(i) = farfield(i) + Ei;
end

%% Plotting the result
% phi = linspace(0,200,length(farfield));
% plot(phi,abs(farfield));
% grid on; 

%calling the volume_disk script for comparison
volume_disk

s = polarplot(onodes, -20*log10(2*pi*oradius*abs(farfield)),'blue');
set(s,'LineWidth',3);


%% function definitions
%%green function
function g = green(rpnl,rpnu,n,rnl,rnu,m,k)
    rps = ((1-n).*rpnl) + (n.*rpnu);
    rt  = ((1-m).*rnl)  + (m.*rnu);
    rho = norm(rt-rps);
    g   = (-1j/4) * besselh(0,2,k.*rho);
end

%%grad_green
function gg = gradg(rpnl,rpnu,n,rnl,rnu,m,nhat,k)
    rps  = ((1-n).*rpnl) + (m.*rpnu);
    rt   = ((1-n).*rnl)  + (m.*rnu);
    rho  = norm(rt-rps);
    rhat = (rt-rps)./rho;
    gg   = (1j*k/4) * besselh(1,2,k.*rho) .* dot(rhat,nhat);
end
