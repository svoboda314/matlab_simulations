
function wh=bendingRods6wPlot(temp, wh)
%
% INPUTS: 
% position along contour in meters; x(1)
% force applied normal to the whisker; x(2)
% wh - structure containing whisker info
% wh.I moment of inertia along the whisker; the number of elements in this
%       array determines the discretization of the simulation
% wh.L whisker length in meters
% wh.ym young's modulus
%
% OUTPUTS:
% contour(1, n_segments) - x values
% contour(2, n_segments) - y values

% external force
s_pole=temp(1); % contour location where force f_pole is applied
f_pole=temp(2)*[wh.friction 1 0]; % force that is applied at s_pole - 3-element vector

% whisker parameters

%wh.ym=3.5e9; % young's modulus in pascal
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

wh.phi=-phi+phi(1);%+pi*30/180; 
%phi=-phi;
%phi=phi+pi*30/180;

wh.x=zeros(1, length(wh.phi)+1);
wh.y=zeros(1, length(wh.phi)+1);
wh.x(2:end)=cumsum(cos(wh.phi)*segment);
wh.y(2:end)=cumsum(sin(wh.phi)*segment);
wh.force_index=force_index;
wh.segment=segment;

wh.f_axial=f_pole*cos(wh.phi(1)-wh.phi(force_index-1)-pi/2);

%  compute the moment
gamma=atan((wh.y(force_index)-wh.y(1))/(wh.x(force_index)- wh.x(1)));
r=sqrt(((wh.y(force_index)-wh.y(1))^2+(wh.x(force_index)-wh.x(1))^2));
beta=abs(wh.phi(force_index)-pi/2-gamma);
wh.moment=r*f_pole*sin(beta);





