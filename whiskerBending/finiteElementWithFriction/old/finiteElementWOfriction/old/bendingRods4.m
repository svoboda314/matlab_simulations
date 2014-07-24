

function polePosError=bendingRods4(x, a)
%
% INPUTS: 
% position along contour in meters; x_pole
% force applied normal to the whisker; force 
% targetPolePos = [xTarget yTarget] 
%
% OUTPUTS:
% distance between pole position and targetPolePos

% target pole position
targetPolePos= a;

% whisker parameters
conical=1;
ym=3.5e9; % young's modulus in pascal
base_radius=1.00e-4; % whisker radius at the base, in meters
L=0.060; % whisker length in meters

% external force
x_pole=x(1);
force=x(2);
%x_pole=0.012; % position where the force is applied, in meters; ALONG THE contour
%force=1/300; % force in Newtons

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

force_index=fix(x_pole/segment); %segment index where the pole applies force

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
%figure(19); cla; hold
phi=phi-phi(1); 
%phi=-phi;
%phi=phi+pi*30/180;

x=zeros(1, length(phi)+1);
y=zeros(1, length(phi)+1);
x(2:end)=cumsum(cos(phi)*segment);
y(2:end)=cumsum(sin(phi)*segment);

%plot(x*1000, y*1000, 'b', 'LineWidth', 1.5);

polePos=[x(force_index), y(force_index)];
%plot(polePos(1)*1000, polePos(2)*1000, 'rx'); % point of force application

polePosError=norm(polePos-targetPolePos);

hold off


