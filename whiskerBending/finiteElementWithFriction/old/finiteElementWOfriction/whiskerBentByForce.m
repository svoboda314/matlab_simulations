
function [polePosError, varargout]=whiskerBentByForce(params, wh, varargin) 
%
% Force (parametrized by f_pole) is applied at position params(1) along the whisker 
% Calculates the contour of the whisker 
%
% INPUTS: 
% all units are MKS (meter, Newton, kilogram, etc)
% position along contour in meters; params(1)
% force applied normal to the whisker; params(2) 
%
% wh - structure containing whisker info
% wh.I - moment of inertia along the whisker; the number of elements in this
%       array determines the discretization of the simulation
% wh.L - whisker length in meters
% wh.ym - young's modulus
% wh.friction - friction coefficient
%
% varargin targetPolePos = [xTarget yTarget] 
%
% OUTPUTS:
% polePosError - distance between pole position and targetPolePos; used for whiskerAgainstObject
% varargout wh - structure containing whisker data
% wh.segment - segment length
% wh.x - x coordinate of the whisker (n+1 points)
% wh.y - y coordinate of the whisker (n+1 points)
% wh.theta - theta angle of each contour segment (n segments)
% wh.f_tot - total force acting on whisker
% wh.f_axial - axial force pushing the whisker into the face along theta base
% wh.f_norm - force acting normal to the whisker at the contact site 
% wh.f_friction - frictional force   
% wh.moment - torque acting on the whisker base
% wh.force_index - index where force acts -- i.e. point of contact along the contour

if nargin > 2
 targetPolePos= varargin{1};
else
 targetPolePos=[];
end

% external force
s_pole=params(1); % contour location where force f_pole is applied
f_pole=[-params(2) wh.friction*params(2) 0]; % force vector in x-y plane CHANGE 011011 - cleaned up signs; still a bit confused by this

% whisker parameters
n_segments=length(wh.I);
segment=wh.L/n_segments; % segment length; total length is n_segments * segment

% initialize arrays
theta=zeros(1, n_segments); % angles of whisker segments wrt to x-axis
x=zeros(1, n_segments+1);% x-position of segments
y=zeros(1, n_segments+1); % y-position of segmetns

force_index=fix(s_pole/segment); %segment index where the pole applies force

x(force_index)=0;
x(force_index-1)=0;
y(force_index)=0;
y(force_index-1)=-segment;
theta(force_index)=0;
theta(force_index-1)=0;

for k=2:(force_index-1) 
    r=[x(force_index-k+1) y(force_index-k+1) 0]; % Note - force is applied at origin 
    moment=(r(1)*f_pole(2)-r(2)*f_pole(1)); % changed sign with change signe in previous line % check sign!!
    theta(force_index-k)=theta(force_index-k+1)+moment*segment/(wh.ym*wh.I(force_index-k+1)); % going backwards
    x(force_index-k)=x(force_index-k+1)-sin(theta(force_index-k))*segment;
    y(force_index-k)=y(force_index-k+1)-cos(theta(force_index-k))*segment;
end

wh.f_axial=dot([sin(theta(1)) -cos(theta(1)) 0], f_pole); % Note: minus sign in fron of sin is because friction pulls 

r=[-(x(force_index)-x(1)) -(y(force_index)-y(1)) 0];
wh.moment = (r(1)*f_pole(2)-r(2)*f_pole(1)); % check sign!!

wh.f_tot=norm(f_pole);
wh.f_friction=abs(f_pole(2));
wh.f_norm=abs(f_pole(1));

theta=theta-theta(1);

x=zeros(1, length(theta)+1);
y=zeros(1, length(theta)+1);
x(2:end)=cumsum(sin(theta)*segment);
y(2:end)=cumsum(cos(theta)*segment);

polePos=[x(force_index), y(force_index)];
if ~ isempty(targetPolePos)
    polePosError=norm(polePos-targetPolePos);
else
    polePosError=NaN;
end

wh.theta=-theta;
wh.x=zeros(1, length(wh.theta)+1);
wh.y=zeros(1, length(wh.theta)+1);
wh.x(2:end)=cumsum(sin(wh.theta)*segment);
wh.y(2:end)=cumsum(cos(wh.theta)*segment);
wh.force_index=force_index;
wh.segment=segment;

varargout{1}=wh;


