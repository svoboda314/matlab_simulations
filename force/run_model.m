function [Z, R, X] = run_model(x0, t, net, varargin)
% function [Z, R, X] = run_model(x0, t, net, varargin)
%
% Return the outputs Z, firing rates R and currents X of simulating the model net for t seconds.  This is an 1xT or
% NxT matrix.
%
% If you change anything in this file, you should also change it in run_compressed_model.m, in the same directory.  They are
% intended to be kept in identical shape. -DCS:2011/01/19
%
optargin = size(varargin,2);

has_input = 0;
TF = 0.0;
F = [];
I = [];
do_net_noise = 1;
for i = 1:2:optargin
    switch varargin{i}
     case 'input'
      I = varargin{i+1}; % Should be an Ixtime matrix
     case 'target'
      F = varargin{i+1};
     case 'TF'
      TF = varargin{i+1};		% TF %, from 0 to 1.
     case 'donetnoise'
      do_net_noise = varargin{i+1};
    end
end
has_target = ~isempty(F);
has_input = ~isempty(I);

% Simple Euler integration of one of the feedback networks trained to run the Romo task.
% -DCS:2010/09/13
N = net.N;				% network size
B = net.B;				% output size (usually M in C implementations)
niters = round(t/net.DT);
gJ = net.g * net.J;
wIn = net.wIn;
wfb = net.wFB;
wout = net.wOut;
dt = net.DT;
tau = net.tau;

dt_div_tau = dt/tau;

currentBias = net.currentBias;

do_use_bias_neuron = net.doUseBiasNeuron;
bias = net.biasValue;
nbidx = net.biasIdx;

X = zeros(N,niters);
R = zeros(N,niters);
Z = zeros(B,niters);



% Set up the initial conditions for the experiment. 
inp = zeros(N,1);
f = zeros(B,1);

x = x0;
if ( do_use_bias_neuron )
    x(nbidx) = bias;
end

net_noise_sigma = net.netNoiseSigma;
act_fun = net.actFun;
r = act_fun(x);

z = wout'*r;
if ( has_target )
    f = F(1:B,1);
end    
zfb = (1-TF)*z + TF*f;			% should be offset by one, but it's taking the first value twice, probably 

%R(:,1) = r;  % These should be saved in the last time step of the data from the last simulation.
%Z(:,1) = z;

for i = 1:niters

    if ( has_input )
	inp = wIn*I(:,i);
    end    
    
    dxdt = -x + gJ*r + wfb*zfb + inp + currentBias;          
    
    if ( do_net_noise & net_noise_sigma > 0.0 )
	dxdt = dxdt + net_noise_sigma * randn(N,1);
    end    

    x = x + dxdt*dt_div_tau;
       
    if ( do_use_bias_neuron )		% Comes last after any other operations to x! 
	x(nbidx) = bias;
    end

    
    
    r = act_fun(x);
    
    z = wout'*r;    
    
    % Handle teacher forcing, if there is any target. 
    if ( has_target )
	f = F(1:B,i);
    end    
    zfb = (1-TF)*z + TF*f;
    
    X(:,i) = x;
    R(:,i) = r;
    Z(:,i) = z;
end

