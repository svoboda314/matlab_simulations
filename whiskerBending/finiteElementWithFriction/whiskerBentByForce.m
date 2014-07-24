
function [polePosError, varargout]=whiskerBentByForce(params, wh, varargin) 
%
% Force (i.e. f_pole) is applied at position params(1) along the whisker. 
% The function calculates the contour of the whisker using the repeated application of the relationship 
% kappa(i) ym I(i)  =  r(i) X F. Here kappa(i) is the curvature at node i, r(i) is the vector connecting node i 
% to the point of contact, I(i) is the second moment of inertia at i, ym is young's modulus. 
% See supplemental methods in Solomon and Hartmann, Nature, 2006.
% Note: For the calculation, the coordinate system (x', y') is such that
% the whisker is parallel to the y' direction at the point of force. The
% normal force points along -x', and the frictional force along y'.
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
% wh.x - x coordinate of the whisker (n+1 points) - 
% wh.y - y coordinate of the whisker (n+1 points) - 
% wh.theta - theta angle of each contour segment (n segments) - theta at base = 0; 
% wh.f_tot - total force acting on whisker
% wh.f_axial - axial force pushing the whisker into the face along theta base
% wh.f_norm - force acting normal to the whisker at the contact site 
% wh.f_friction - frictional force   
% wh.moment - torque acting on the whisker base
% wh.force_index - index where force acts -- i.e. point of contact along the contour
%
% KS 011011
%
if nargin > 2
 targetPolePos= varargin{1};
else
 targetPolePos=[];
end

% force paramters
s_pole=params(1); % contour location along the whisker where force f_pole is applied
f_pole=[-params(2) wh.friction*params(2) 0]; % force vector in x', y' plane 

% whisker parameters
n_segments=length(wh.I);
segment=wh.L/n_segments; % segment length; total length is n_segments * segment
force_index=fix(s_pole/segment); %segment index where the pole applies force

% initialize arrays
theta=zeros(1, n_segments); % angles of whisker segments wrt to x-axis
x=zeros(1, n_segments+1);% x-position of segments
y=zeros(1, n_segments+1); % y-position of segmetns
x(force_index)=0;
x(force_index-1)=0;
y(force_index)=0;
y(force_index-1)=-segment;
theta(force_index)=0;
theta(force_index-1)=0;

% Apply the relationship kappa(i) ym I(i)  =  r(i) X F iteratively,
% starting at the point of contact
for k=2:(force_index-1) 
    r=[x(force_index-k+1) y(force_index-k+1) 0]; 
    moment=(r(1)*f_pole(2)-r(2)*f_pole(1)); 
    theta(force_index-k)=theta(force_index-k+1)+moment*segment/(wh.ym*wh.I(force_index-k+1)); % going backwards towards the follicle
    x(force_index-k)=x(force_index-k+1)-sin(theta(force_index-k))*segment;
    y(force_index-k)=y(force_index-k+1)-cos(theta(force_index-k))*segment;
end

wh.f_axial=dot([sin(theta(1)) -cos(theta(1)) 0], f_pole); % Note: signs are because friction pulls, normal force pushes 

r=[-(x(force_index)-x(1)) -(y(force_index)-y(1)) 0];
wh.moment = (r(1)*f_pole(2)-r(2)*f_pole(1)); 
wh.f_tot=norm(f_pole);
wh.f_friction=abs(f_pole(2));
wh.f_norm=abs(f_pole(1));

theta=theta-theta(1); % rotate whisker so that the base is along the y axis

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


