
function [polePosError, varargout]=bendingRods5(temp, wh, varargin) % targetPolePos)
%
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

if nargin > 2
 targetPolePos= varargin{1};
else
 targetPolePos=[];
end

% external force
s_pole=temp(1); % contour location where force f_pole is applied
% f_poleDUMMY=temp(2); % force that is applied at s_pole
f_pole=[wh.friction*temp(2) -temp(2) 0]; % <<<<<<<<<<<NEW

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
%    rDUMMY=sqrt(((y(force_index-k+1))^2+(x(force_index-k+1))^2)); % U
%    beta=acot((x(force_index-k+1))/(y(force_index-k+1)))+pi/2; % U
    r=[x(force_index-k+1) y(force_index-k+1) 0]; % <<<<<<<<<<<NEW
    moment=(r(1)*f_pole(2)-r(2)*f_pole(1));% <<<<<<<<<<<NEW
%    momentDUMMY=rDUMMY*f_poleDUMMY*sin(beta) % U
    phi(force_index-k)=phi(force_index-k+1)-moment*segment/(wh.ym*wh.I(force_index-k+1)); % going backwards
    x(force_index-k)=x(force_index-k+1)-cos(phi(force_index-k))*segment;
    y(force_index-k)=y(force_index-k+1)-sin(phi(force_index-k))*segment;
    %plot( force_index-k, sin(phi(force_index-k)))
end

wh.phi=phi-phi(1); 
%phi=-phi;
%phi=phi+pi*30/180;

wh.x=zeros(1, length(wh.phi)+1);
wh.y=zeros(1, length(wh.phi)+1);
wh.x(2:end)=cumsum(cos(phi)*segment);
wh.y(2:end)=cumsum(sin(phi)*segment);

polePos=[wh.x(force_index), wh.y(force_index)];
if ~ isempty(targetPolePos)
    polePosError=norm(polePos-targetPolePos);
else
    polePosError=NaN;
end

wh.force_index=force_index;
wh.segment=segment;

% wh.f_axial=f_pole*cos(wh.phi(1)-wh.phi(force_index-1)-pi/2); % FIX

wh.f_axial= dot([sin(wh.phi(1)) sin(wh.phi(1)) 0]/sqrt(2), f_pole);
r=[(wh.x(force_index)-wh.x(1)) (wh.y(force_index)-wh.y(1)) 0];
wh.moment = (r(1)*f_pole(2)-r(2)*f_pole(1));
varargout{1}=wh;
% figure(17); cla; hold
% plot(x*1000, y*1000, 'b', 'LineWidth', 1.5);
% plot(polePos(1)*1000, polePos(2)*1000, 'rx'); % point of force application
% hold off



