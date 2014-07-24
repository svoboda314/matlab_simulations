% Script to compute contour of whisker 
%
% Whisker starts along the x axis without curvature
% Deflected upward by a pole at position s_pole providing a force.
% Methods similar to Solomon and Hartmann, Nature
%
% r - vector from current length element and point of force application
% beta - angle between vertical and r
% force_angle - angle between force and vertical axis
%

% whisker parameters
conical=1;
ym=3.5e9; % young's modulus in pascal
base_radius=1.00e-4; % whisker radius at the base, in meters
L=0.060; % whisker length in meters
s_pole=0.012; % position along the whisker contour where the force is applied, in meters

n_segments=1000;
segment=L/n_segments; % segment length in meters; total length is n_segments * segment

% initialize arrays
phi=zeros(1, n_segments); % angles of whisker segments wrt to x-axis
x=zeros(1, n_segments+1);% x-position of segments
y=zeros(1, n_segments+1); % y-position of segmetns

% Whisker geometry
if ~ conical
    I(1:n_segments)=(pi/4)*base_radius^4; %uniform
else
    I(1:n_segments)=(pi/4)*(base_radius*(1-(1:n_segments)*segment/L)).^4;  %conical
end

force_index=fix(s_pole/segment); %segment index where the pole applies force
force=1/300; % force in Newtons

x(force_index)=0;
x(force_index-1)=-segment;
y(force_index)=0;
y(force_index-1)=0;
phi(force_index)=0;
phi(force_index-1)=0;

for k=2:(force_index-1) 
    r=sqrt(((y(force_index-k+1))^2+(x(force_index-k+1))^2)); % U
    beta=acot((x(force_index-k+1))/(y(force_index-k+1)))+pi/2; % U
    moment=r*force*sin(beta); % U
    phi(force_index-k)=phi(force_index-k+1)-moment*segment/(ym*I(force_index-k+1)); % going backwards
    x(force_index-k)=x(force_index-k+1)-cos(phi(force_index-k))*segment;
    y(force_index-k)=y(force_index-k+1)-sin(phi(force_index-k))*segment;
    %plot( force_index-k, sin(phi(force_index-k)))
end
figure(19); cla; hold
phi=phi-phi(1); 
%phi=-phi;
%phi=phi+pi*30/180;

x=zeros(1, length(phi)+1);
y=zeros(1, length(phi)+1);
x(2:end)=cumsum(cos(phi)*segment);
y(2:end)=cumsum(sin(phi)*segment);

plot(x*1000, y*1000, 'b', 'LineWidth', 1.5);
plot(x(force_index)*1000, y(force_index)*1000, 'rx'); % point of force application
x(force_index)
%Plot small angle approximations
x1=(0:n_segments-1)*segment;
y1=zeros(1, (n_segments));
temp1=y1;
%x(1:force_index); %(0:n_segments)/(n_segments+1)*s_pole; %distance in meters
x_poleL=x1(force_index);
a=x_poleL;
if ~ conical
    y1(1:force_index)=(force./(6*ym*I(1:force_index))).*(3*x_poleL*x1(1:force_index).^2-x1(1:force_index).^3); % uniform thin beam model
    y1((force_index+1):end)=(force./(6*ym*I((force_index+1):end))).*(3*x_poleL^2*x1((force_index+1):end)-x_poleL^3);
else
    for i=1:force_index
    y1(i)=(2*force*L*x1(i).^2/(3*ym*pi*(base_radius^4)))*((3*L*x_poleL-L*x1(i)-2*x_poleL*x1(i))/(L-x1(i)).^2); % tapered thin beam model
    temp1(i)=(L-x1(i)).^2;
    end
    y12(1:force_index)=((2*force*L*x1(1:force_index).^2)/(3*ym*pi*(base_radius^4))).*((3*L*x_poleL-L*x1(1:force_index)-2*x_poleL*x1(1:force_index))./(L-x1(1:force_index)).^2); 
    temp12=(L-x1(1:force_index)).^2;
%    y1(1:force_index)=(2*force*L*x1(1:force_index).^2/(3*ym*pi*base_radius^4))*((3*L*s_pole-L*x1(1:force_index)-2*s_pole*x1(1:force_index))/(L-x1(1:force_index)).^2); % tapered thin beam model
    y1((force_index+1):end)=(2*force*L*(x_poleL^2)/(3*ym*pi*base_radius^4))*(((3*L*x1((force_index+1):end))-L*x_poleL-2*x_poleL*x1((force_index+1):end))/(L-x_poleL).^2); % tapered thin beam model    
    %            b=base_radius/L;
    %    y12=(4*force/(pi*ym))*(((b*(L-3*x1)+2*base_radius)./(6*b^3*(base_radius-b*x1).^2))-(r+2*b*L)*x1/(6*b^2*base_radius^3)-(b*L+2*base_radius)/6*b^3*base_radius^2);
    
end
% plot(x1*1000, y1*1000, 'm')
% b=base_radius/L;
% y12=(4*force/(pi*ym))*(((b*(L-3*x1)+2*base_radius)./(6*b^3*(base_radius-b*x1).^2))-(r+2*b*L)*x1/(6*b^2*base_radius^3)-(b*L+2*base_radius)/6*b^3*base_radius^2);
plot(x1*1000, y1*1000, '--k', 'LineWidth', 1.0)

hold off
