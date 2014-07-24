function r = LHFI(x, R0)
% See Rajan, Abbott, Sompolinsky, PRLE 2010 
%
% For r0 = 1, we recover the often-used tanh function, but we use a smaller
% value of r0 = 0.1, which is more biologically reasonable.

% The tanh function has the disadvantage of having the “resting” rate φ(0) halfway between the minimum and maximum
% rates.  This generalization allows us to adjust the value of φ(0) to be closer to the minimum of this range, while
% retaining the desirable feature that the maximum of the derivative of φ is at x = 0.

nx = length(x);
r = zeros(nx,1);
for i = 1:nx
    xi = x(i);
    R0i = R0(i);
    if  xi <= 0.0    
	r(i) =  R0i*tanh(xi./R0i); % + R0;
    else   
	r(i) = (2.0-R0i)*tanh(xi./(2.0-R0i));  % + R0;
    end
end

