
Fo_cell = 100; 
A=0.2 % DF/F ...
DF = A*Fo_cell;

alpha = 1.00 % neuropil bleedthrough
Fo_neuropil = Fo_cell;

Fo = Fo_cell + alpha*Fo_neuropil; 

samples=100000;

bl=poissrnd(Fo, samples, 1);
resp=poissrnd((Fo+DF), samples, 1);

x=min(bl):1:max(resp);

n_bl=hist(bl, x);
n_bl=n_bl/sum(n_bl);
n_resp=hist(resp, x);
n_resp=n_resp/sum(n_resp);

if (exist('hf') & hf~=1)
hf=figure(1)
set(hf, 'Position', [100 100 1100 500] );
hal=axes('Position', [0.08 0.1 0.4 0.8], 'FontSize', 16);
har=axes('Position', [0.58 0.1 0.4 0.8], 'FontSize', 16 );
else
    cla(hal)
    cla(har)
end
    

line(x, n_bl, 'Parent', hal, 'LineWidth', 1.5, 'Color', [0 0 0]);
line(x, n_resp, 'Parent', hal, 'LineWidth', 1.5);


%plot(x, n_bl, 'k')
%hold
%plot(x, n_resp, 'b')
%hold off

% far=zeros(length(x)-1, 1);
% hr=far;
% for i=1:(length(x)-1)
%    far(i)=sum(n_bl(i:end));
%    hr(i)=sum(n_resp(i:end));
% end
% 
% plot(far, hr, 'o')
% hold
% plot(far, far, 'k')
% 
% auc=0;
% for i=1:(length(far)-1)
%     auc=auc+(far(i+1)-far(i))*(hr(i+1)+hr(i))/2;
% end
% auc = abs(auc)-1/2;

sc=[bl',resp'];
la=[zeros(1,samples), ones(1,samples)];

[X, Y, T, auc]=perfcurve(la, sc, 1);
auc

line(X, Y, 'Parent', har, 'LineWidth', 1.5, 'Color', [0 0 0]);
line(0:.1:1,0:.1:1,  'Parent', har, 'LineWidth', 1.0, 'LineStyle', '--','Color', [0 0 0]);
% [X,Y] = perfcurve(resp,bl)

% plot(bl);
% hold
% plot(resp, 'g')

