
function f = twoPdepth(z, ls)
% Equation from Theer and Denk 2003 Opt Lett
% set NA = n = lambda = 1

f=ls/(z*z*2*pi) - exp(-2*z/ls);
