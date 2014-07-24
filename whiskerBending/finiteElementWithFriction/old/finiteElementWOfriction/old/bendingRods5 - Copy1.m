
function [polePosError, varargout]=bendingRods5(temp, wh, varargin) %targetPolePos
% INPUTS: 
% position along contour in meters; temp(1)
% force applied normal to the whisker; temp(2)
% wh - structure containing whisker info
% wh.I moment of inertia along the whisker; the number of elements in this
%       array determines the discretization of the simulation
% wh.L whisker length in meters
% wh.ym young's modulus
% varargin targetPolePos = [xTarget yTarget] 
%
% OUTPUTS:
% distance between pole position and targetPolePos
% varargout wh - structure containing whisker data

if nargin > 2
 targetPolePos= varargin{1};
else
 targetPolePos=[];
end

% external force
s_pole=temp(1); % contour location where force f_pole is applied
f_pole=[wh.friction*temp(2) temp(2) 0]; % <<<<<<<<<<<NEW

% whisker parameters
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
    r=[x(force_index-k+1) y(force_index-k+1) 0]; 
    moment=-(r(1)*f_pole(2)-r(2)*f_pole(1));
    phi(force_index-k)=phi(force_index-k+1)-moment*segment/(wh.ym*wh.I(force_index-k+1)); % going backwards
    x(force_index-k)=x(force_index-k+1)-cos(phi(force_index-k))*segment;
    y(force_index-k)=y(force_index-k+1)-sin(phi(force_index-k))*segment;
end


wh.f_axial=dot([sin(phi(1)) cos(phi(1)) 0], f_pole); % ALREADY NEW COORDINATE SYSTEM
wh.moment=0; % AlREADY NEW COORDINATE SYSTEM
r=[(x(force_index)-x(1)) (y(force_index)-y(1)) 0];
wh.moment = (r(1)*f_pole(2)-r(2)*f_pole(1));

phi=phi-phi(1);

x=zeros(1, length(phi)+1);
y=zeros(1, length(phi)+1);
x(2:end)=cumsum(cos(phi)*segment);
y(2:end)=cumsum(sin(phi)*segment);

polePos=[x(force_index), y(force_index)];
if ~ isempty(targetPolePos)
    polePosError=norm(polePos-targetPolePos);
else
    polePosError=NaN;
end

wh.phi=-phi;
wh.x=zeros(1, length(wh.phi)+1);
wh.y=zeros(1, length(wh.phi)+1);
wh.x(2:end)=cumsum(cos(wh.phi)*segment);
wh.y(2:end)=cumsum(sin(wh.phi)*segment);
wh.force_index=force_index;
wh.segment=segment;

varargout{1}=wh;


