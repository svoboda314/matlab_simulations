function [net, varargout] = learn_model(x0, net, F, varargin)
% function net = learn_model(net, F, varargin)
% function [net, learn] = learn_model(net, F, varargin)
% function [net, learn, simdata] = learn_model(net, F, varargin)

% Variable output considerations.
nout = max(nargout,1)-1;
do_save_network_data = 0;
if ( nout > 0 )
    do_save_network_data = 1;
end

% Variable input considerations.
optargin = size(varargin,2);
linewidth = 3;
fontsize = 14;
fontweight = 'bold';

alpha = 1.0;				% regularization constant for RLS, it diagonally loads the inverse
                                        % correlation matrix.
has_input = 0;
TF = [];
DL = [];
Inp = [];
learn_struct = [];
do_plot = 0;
npasses = 1;
learn_every = 1;
for i = 1:2:optargin
    switch varargin{i}
     case 'learnstruct'
      learn_struct = varargin{i+1};		% Reuse the same struct, since the P has meaningful info.
     case 'learnevery'
      learn_every = varargin{i+1};
     case 'alpha'
      alpha = varargin{i+1};
     case 'npasses'
      npasses = varargin{i+1};		% N passes through the training. 
     case 'input'
      Inp = varargin{i+1}; % Should be an Ixtime matrix
     case 'target'			% Can put the target function in either way.
      F = varargin{i+1};
     case 'TF'
      TF = varargin{i+1};		% Soft teacher forcing so TF is from 0 to 1, for each time step, or a simple
                                        % 0 or 1 
     case 'DL'
      DL = varargin{i+1};		% Do we learn on this particular time step, [0,1] for each time step, or
                                        % a simple 0 or 1.
     case 'doplot'
      do_plot = varargin{i+1};
    end
end

N = net.N;				% network size
B = net.B;				% output size (usually M in C implementations)
I = net.I;
niters = length(F);

% BULLET PROOF THAT ALL THE INPUTS, TARGETS, DL, AND TF HAVE THE SAME LENGTH. 
% TF can be either a single row, empty (so always 0), or specified at each time point by the user.
if ( isempty(TF) )
    TF = zeros(B,niters);
end
if ( size(TF,2) == 1 )
    TF = repmat(TF, 1, niters);
end
assert ( length(TF) == niters, 'All learning input vectors should have the same length.');
assert ( size(TF,1) == B, 'There must be a TF row for each output');
if ( isempty(DL) )
    DL = ones(B,niters);
end
% DL can be either a single row, empty (so always 1), or specified at each time point by the user.
if ( size(DL,2) == 1 )
    DL = repmat(DL, 1, niters);
end
assert ( length(DL) == niters, 'All learning input vectors should have the same length.');
assert ( size(DL,1) == B, 'There must be a DL row for each output');
% The input can be either empty, or specified at each time point by the user.
has_input = ~isempty(Inp);
if ( has_input )
    assert ( length(Inp) == niters, 'All learning input vectors should have the same length.');
    assert ( size(Inp,1) == I, 'There must be an input entry for each input vector.');
end


% Unpack the structure for clearness and a little speed. 
gJ = net.g * net.J;
win = net.wIn;
wfb = net.wFB;
wout = net.wOut;
dt = net.DT;
tau = net.tau;

dt_div_tau = dt/tau;

net_noise_sigma = net.netNoiseSigma;
currentBias = net.currentBias;

act_fun = net.actFun;

do_use_bias_neuron = net.doUseBiasNeuron;
bias = net.biasValue;
nbidx = net.biasIdx;


% Setup the learning structure. 
if ( isempty(learn_struct) )
    P = (1.0/alpha)*eye(N);
else
    P = learn_struct.P;
end

% Initialize the memeory used in the learning. 
if ( do_save_network_data )
    simdata.X = zeros(N,niters);
    simdata.R = zeros(N,niters);
    simdata.Z = zeros(B,niters);
    simdata.W = zeros(N,B,niters);
end


% Set up the initial conditions for the experiment. 
inp = zeros(N,1);
f = zeros(B,1);
x = x0;
if ( do_use_bias_neuron )
    x(nbidx) = bias;
end
r = act_fun(x);
z = wout'*r;
f = F(1:B,1);
tf = TF(:,1);
zfb = (1-tf).*z + tf.*f;			% should be offset by one, but it's taking the first value twice, probably 

%R(:,1) = r;  % These should be saved in the last time step of the data from the last simulation.
%Z(:,1) = z;

for n = 1:npasses
    for i = 1:niters
	if ( has_input )
	    inp = win*Inp(:,i);
	end    
	
	% x
	dxdt = -x + gJ*r + wfb*zfb + inp + currentBias;          
	if ( net_noise_sigma > 0.0 )
	    dxdt = dxdt + net_noise_sigma*randn(N,1);
	end
	
	x = x + dxdt*dt_div_tau;
	
	if ( do_use_bias_neuron )		% Always goes after the integration.
	    x(nbidx) = bias;
	end
	% r & z
	r = act_fun(x);    
	z = wout'*r;    
	
	% Recursive Least Squares learning implements FORCE
	% See page 443 in Haykin's adaptive filter theory (3rd Ed). 
	f = F(:,i);
	do_learn = sum(DL(:,i)) > 0;	% Does any output need to learn this time step?
	do_learn = do_learn & mod(i, learn_every) == 0; % Is this a time step that we learn at all?
	if  do_learn
	    % update inverse correlation matrix
	    k = P*r;
	    rPr = r'*k;
	    c = 1.0/(1.0 + rPr);
	    k = k*c;
	    P = P - k*(k'/c);
	    
	    % update the error for the linear readout
	    err = z-f;
	    
	    % update the output weights
	    kall = repmat(k,1,B);	% only works under the assumption that all units get all n outputs, which is
                                        % true in create_random_model.
	    % DL simply makes sure all w vectors "want" to learn on this time step, by multiplying by 0 if not. 
	    dw = -bsxfun(@times, kall, (DL(:,i).*err)');
	    wout = wout + dw;	
	end      
	
	% Handle teacher forcing, if there is any target. 
	tf = TF(:,i);
	zfb = (1-tf).*z + tf.*f;
	
	if ( do_save_network_data )
	    simdata.X(:,i) = x;
	    simdata.R(:,i) = r;
	    simdata.Z(:,i) = z;
	    simdata.W(:,:,i) = wout;
	end
    end
end

% Put the learned part back in the network structure.
net.wOut = wout;
learn_struct.P = P;


if ( nout >= 1 )
    varargout{1} = learn_struct;
end
if ( nout >=2 ) 
    varargout{2} = simdata;
end