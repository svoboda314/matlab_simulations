


function out=psf1d_shell(x, x_p, s, r, d) 
% x_p center of psf
% s - half width
% r - radius of cell
% d - thickness of shell (< r)



%x=-100:100;
%d=2;
%r=20;
out=zeros(1, length(x));
n1=find(x<-r);
n2=find(x>-r+d & x < r-d);
n3=find(x>r);
n_zero=[n1, n2, n3];
out(n_zero)=0;
n=setdiff(1:length(x), n_zero);

out(n)=(1/(s*sqrt(2*pi)))*exp(-(x(n)-x_p).^2/(2*s^2));

