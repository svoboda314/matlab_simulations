function net = create_random_model(N, B, I, p, g, dt, tau, varargin)
% N - the number of recurrent neurons in network
% B - the number of outputs
% I - the number of inputs
% p - the sparseness of the J (connectivity) matrix, e.g. 0.1
% g - the spectral scaling of J
% dt - the integration time constant
% tau - the time constant of each neuron
% net_noise_sigma - random gaussian noise to the dxdt at each dt
% current_bias_sigma - magnitude of current offset to each neuron
% act_fun_type - 0 for tanh, or 5 for the PRLE 2010 function.
net_noise_sigma = 0.0;
current_bias_sigma = 0.0; 
act_fun_type = 0;

optargin = size(varargin,2);
for i = 1:2:optargin
    switch varargin{i}
     case 'netnoisesigma'
      net_noise_sigma = varargin{i+1};
     case 'currentbiassigma'
      current_bias_sigma = varargin{i+1};
     case 'actfuntype'
      act_fun_type = varargin{i+1};
    end
end

scale = 1.0/sqrt(p*N);
J = sprandn(N,N,p)*scale;
J = full(J);

net.I = I;
net.B = B;
net.N = N;
net.sep1 = '-------------';

net.g = g;
net.J = J;
net.h = 0.0;				% lateral eigenshift, zero no shift.
net.wOut = zeros(N,B);			% doesn't have to be zero.
net.wFB = 2.0*(rand(N,B)-0.5);
net.wIn = 2.0*(rand(N,I)-0.5);

net.sep2 = '-------------';
net.DT = dt;				
net.tau = tau;
net.doUseBiasNeuron = false;		% Should one neuron always give a constant value?
net.biasIdx = -1;
net.biasValue = 0.0;
net.netNoiseSigma = net_noise_sigma;
net.currentBias = current_bias_sigma * randn(N,1);

net.sep3 = '-------------';
net.actFunType = act_fun_type;

% Reference forcelib/network/network_utils.h for the values and meanings of this enumerated type. 
switch net.actFunType
 case 0
  net.R0 = repmat([1.0], N, 1);				% Should be a parameter. -DCS:2011/03/17  
  net.actFun = @tanh;			% Good old iron-sides.
  net.actFunDerivX = @(x) 1 - tanh(x).^2;
  net.actFunDerivR = @(r) 1 - r.^2;
 case 5					% Kanaka, Larry & Haim's function. 
  net.R0 = repmat([0.1], N, 1);				% Should be a parameter. -DCS:2011/03/17a
  net.actFun = @(x) LHFI(x, net.R0);	% this could be problematic if we want to use the function in other ways. -DCS:2011/03/12
  net.actFunDerivX = @(x) LHFI_deriv(x, net.R0); % I think I may have to redefine these functions if R0 changes. Not sure.
  net.actFunDerivR = @(r) LHFI_deriv_r(r, net.R0);

 otherwise
  assert ( false, 'Case not implemented yet.');
end
