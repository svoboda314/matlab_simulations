
function w_shape=calcWhiskerShape(a, F, varargin)
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
p.addParamValue('npts', 10000, @(x) x > 10 && x < 1e7);
p.parse(a, F, varargin{:});

% convert to SI units, meters, Newtons
rb=p.Results.rb/1e6;
L=p.Results.L/1000;
a=p.Results.a/1000;
ym=p.Results.ym*1e9;
F=p.Results.F/1e6;
npts=p.Results.npts;

% some other checks
if a > L
    error('a > L; whisker does not contact object')
end



% calculation parameters
% npts=10000;
x=(0:npts)*L/npts;
n_a=round((npts+1)*a/L);   % index corresponding to application of force
y=zeros(1, npts+1);

if ~ p.Results.conical
    % cylndrical model; check if npts > n_a
    I = pi*rb^4/4;
    y(1:n_a)=(F/(6*ym*I))*(3*x(1:n_a).^2*a-x(1:n_a).^3);
    y((n_a+1):end)=(F/(6*ym*I))*(3*x((n_a+1):end)*a^2-a^3);
else
    % conical model
    y(1:n_a)= ((2*F*L*x(1:n_a).^2)/(3*ym*pi*rb^4)).*(3*L*a-L*x(1:n_a)-2*a*x(1:n_a))./((L-x(1:n_a)).^2);
    y((n_a+1):end)=((2*F*L*a^2)/(3*ym*pi*rb^4)).*(3*L*x((n_a+1):end)-L*a-2*a*x((n_a+1):end))/((L-a).^2);
end

w_shape=[x*1000; y*1000];