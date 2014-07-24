
n=100; % number of loops
R=0.05; % radius in meters
mu=4*pi*10e-7; %permeability in Tm/A
x=0.05; % distance along the axis from the center in meters
I=1; % current in A

R=(1:10)/100;

B=mu*n*I*R.^2./(2*(R.^2+x^2).^1.5)*10e4;
plot(R, B)