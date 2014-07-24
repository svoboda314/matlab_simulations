


function out=psf1d(x, x_p, s) 

out=(1/(s*sqrt(2*pi)))*exp(-(x-x_p).^2/(2*s^2));

