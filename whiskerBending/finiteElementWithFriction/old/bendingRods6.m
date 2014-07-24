
function polePosError=bendingRods6(temp, targetPolePos, wh)
% Differs from bendingRods5 by addition of friction
% INPUTS: 
% position along contour in meters; temp(1)
% force applied normal to the whisker; temp(2)
% targetPolePos = [xTarget yTarget] 
% wh - structure containing whisker info
% wh.I moment of inertia along the whisker; the number of elements in this
%       array determines the discretization of the simulation
% wh.L whisker length in meters
% wh.ym young's modulus
%
% OUTPUTS:
% distance between pole position and targetPolePos

% external force
s_pole=temp(1); % contour location where force f_pole is applied
f_pole=temp(2)*[wh.friction 1 0]; % force that is applied at s_pole - 3-element vector
% fc=temp(3); % friction coefficient < 1 

% whisker parameters

% wh.ym=3.5e9; % young's modulus in pascal
%base_radius=1.00e-4; % whisker radius at the base, in meters

%s_pole=0.012; % position where the force is applied, in meters; ALONG THE contour
%f_pole=1/300; % force in Newtons

n_segments=length(wh.I);
segment=wh.L/n_segments; % segment length in meters; total length is n_segments * segment

% initialize arrays
phi=zeros(1, n_segments); % angles of whisker segments wrt to x-axis
x=zeros(1, n_segments+1);% x-position of segments
y=zeros(1, n_segments+1); % y-position of segmetns

force_index=fix(s_pole/segment); %segment index where the pole applies force

x(force_index)=0;
x(force_index-1)=-segment;
y(force_index)=0;
y(force_index-1)=0;
phi(force_index)=0;
phi(force_index-1)=0;

for k=2:(force_index-1) 
    %r=sqrt(((y(force_index-k+1))^2+(x(force_index-k+1))^2)); % U
    %beta=acot((x(force_index-k+1))/(y(force_index-k+1)))+pi/2; % U
    r=[x(force_index-k+1) y(force_index-k+1) 0];
    moment=norm(cross(r, f_pole));
    %moment=r*f_pole*sin(beta); % U
    phi(force_index-k)=phi(force_index-k+1)-moment*segment/(wh.ym*wh.I(force_index-k+1)); % going backwards
    x(force_index-k)=x(force_index-k+1)-cos(phi(force_index-k))*segment;
    y(force_index-k)=y(force_index-k+1)-sin(phi(force_index-k))*segment;
    %plot( force_index-k, sin(phi(force_index-k)))
end

% CHANGE

phi=phi-phi(1); 
%phi=-phi;
%phi=phi+pi*30/180;

x=zeros(1, length(phi)+1);
y=zeros(1, length(phi)+1);
x(2:end)=cumsum(cos(phi)*segment);
y(2:end)=cumsum(sin(phi)*segment);

polePos=[x(force_index), y(force_index)];
polePosError=norm(polePos-targetPolePos);

% figure(17); cla; hold
% plot(x*1000, y*1000, 'b', 'LineWidth', 1.5);
% plot(polePos(1)*1000, polePos(2)*1000, 'rx'); % point of force application
% hold off



