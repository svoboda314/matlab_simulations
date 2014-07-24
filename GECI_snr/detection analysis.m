
Fo_cell = 10; 
%A=0.2 % DF/F ...
DF = 4;
broadening = 1; % neuropil noise is bigger than shotnoise
alpha = 0 % neuropil bleedthrough
Fo_neuropil = Fo_cell;

Fo = Fo_cell + alpha*Fo_neuropil; 

samples=100000;

%bl=poissrnd(Fo, samples, 1);

bl=poissrnd(broadening*Fo, samples, 1);
bl=bl-(broadening-1)*Fo;

resp=poissrnd(broadening*(Fo+DF), samples, 1);
resp=resp-(broadening-1)*(Fo+DF);

x=min(bl):1:max(resp);

n_bl=hist(bl, x);
n_bl=n_bl%/sum(n_bl);
n_resp=hist(resp, x);
n_resp=n_resp%/sum(n_resp);

if ~exist('hf')
hf=figure(11);
set(hf, 'Position', [100 100 1100 500]);
hal=axes('Position', [0.08 0.1 0.4 0.8], 'FontSize', 16, 'Linewidth', 1, 'YtickLabel', [], 'XLim', [0 300]);
xlabel('Signal photons');
ylabel('Measurements');
har=axes('Position', [0.58 0.1 0.4 0.8], 'FontSize', 16, 'Linewidth', 1);
xlabel('False alarm rate');
ylabel('Hit rate');
else
    cla(hal)
    cla(har)
end
    
line(x, n_bl, 'Parent', hal, 'LineWidth', 1.5, 'Color', [0 0 0]);
line(x, n_resp, 'Parent', hal, 'LineWidth', 1.5);


sc=[bl',resp'];
la=[zeros(1,samples), ones(1,samples)];

[X, Y, T, auc]=perfcurve(la, sc, 1);
auc

line(X, Y, 'Parent', har, 'LineWidth', 5, 'Color', [1 0 0]);
line(0:.1:1,0:.1:1,  'Parent', har, 'LineWidth', 1, 'LineStyle', '--','Color', [0 0 0]);
% [X,Y] = perfcurve(resp,bl)

% plot(bl);
% hold
% plot(resp, 'g')

