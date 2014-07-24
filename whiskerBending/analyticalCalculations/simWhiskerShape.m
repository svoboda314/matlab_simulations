
function w_shape=simWhiskerShape(a, F, varargin) % (x, wh)

% Computes whisker shape with a force applied at a point a along the contour of the whisker
% Uses small-angle approcimation (see Birdwell JA et al JNP 2007)
%
% 
% INPUTS  
% a; position at which the force is exerted in mm; a < L
% F; force applied at a in uN
%
% parameter-value pairs - ALL OPTIONAL
% ym; Young's modulus - Young's modulus GPa (1 Pa == 1 N/m^2)
% L - length of the whisker in mm
% rb - radius at the base im um
% conical - 1 == conical model; 0 = cylindrical model
%
% OUTPUTS
% contour of the whisker shape - in millimeters

p = inputParser;   % Create instance of inputParser class.
p.addRequired('a', @(x)x>1 && x<100); % mouse 5
p.addRequired('F', @(x)x>0 && x<10000); % mouse 50; 
p.addParamValue('L', 16, @(x)x>0 && x<100); % mouse 16; rat 60
p.addParamValue('ym', 5, @(x)x>0 && x<10); % mouse 5 GPa; rat 2.75 GPa
p.addParamValue('rb', 33.5, @(x)x>0 && x<1000); % mouse 33.5; rat 100
p.addParamValue('conical', 1, @(x) ismember(x, [0 1]));
p.addParamValue('n_segments', 10000, @(x) x > 10 && x < 1e7);
p.parse(a, F, varargin{:});

% convert to SI units, meters, Newtons
rb=p.Results.rb/1e6;
L=p.Results.L/1000;
a=p.Results.a/1000;
ym=p.Results.ym*1e9;
F=p.Results.F/1e6;
n_segments=p.Results.n_segments;

segment=L/n_segments;

% some other checks
if a > L
    error('a > L; whisker does not contact object')
end

%calculate moment of inertia I
if ~ p.Results.conical
    I(1:n_segments)=(pi/4)*rb^4; %uniform
else
    I(1:n_segments)=(pi/4)*(rb*(1-(1:n_segments)*segment/L)).^4;  %conical
end

% external force
s_pole=a; % contour location where force f_pole is applied
f_pole=F; % force that is applied at s_pole

% whisker parameters

%ym=3.5e9; % young's modulus in pascal
%rb=1.00e-4; % whisker radius at the base, in meters

%s_pole=0.012; % position where the force is applied, in meters; ALONG THE contour
%f_pole=1/300; % force in Newtons

%n_segments=length(wh.I);sim
 % segment length in meters; total length is n_segments * segment

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
    r=sqrt(((y(force_index-k+1))^2+(x(force_index-k+1))^2)); % U
    beta=acot((x(force_index-k+1))/(y(force_index-k+1)))+pi/2; % U
    moment=r*f_pole*sin(beta); % U
    phi(force_index-k)=phi(force_index-k+1)-moment*segment/(ym*I(force_index-k+1)); % going backwards
    x(force_index-k)=x(force_index-k+1)-cos(phi(force_index-k))*segment;
    y(force_index-k)=y(force_index-k+1)-sin(phi(force_index-k))*segment;
    %plot( force_index-k, sin(phi(force_index-k)))
end

phi=-phi+phi(1);%+pi*30/180; 
%phi=-phi;
%phi=phi+pi*30/180;

x=zeros(1, length(phi)+1);
y=zeros(1, length(phi)+1);
x(2:end)=cumsum(cos(phi)*segment);
y(2:end)=cumsum(sin(phi)*segment);
w_shape=[x*1000; -y*1000];
% %wh.force_index=force_index;
% %wh.segment=segment;
% 
% f_axial=f_pole*cos(phi(1)-phi(force_index-1)-pi/2);
% 
% %  compute the moment
% gamma=atan((wh.y(force_index)-wh.y(1))/(wh.x(force_index)- wh.x(1)));
% r=sqrt(((wh.y(force_index)-wh.y(1))^2+(wh.x(force_index)-wh.x(1))^2));
% beta=abs(wh.phi(force_index)-pi/2-gamma);
% wh.moment=r*f_pole*sin(beta);





