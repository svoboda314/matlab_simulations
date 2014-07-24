
figure; hold
mu=.1:.1:25;
n=0:100; % num
p_see=zeros(10, length(mu));
for i=1:length(mu) % i is the average number of photons per flash
    pn_mu = poisspdf(n,mu(i));
    for j=1:size(p_see, 1)
        p_see(j, i)=sum(pn_mu((j+1):end));
    end
end
semilogx(mu, p_see, 'k', 'LineWidth', 1.5)
xlabel('Log avg number of absorbed photons', 'fontsize',16)
ylabel('P(n>K)', 'fontsize',16)
title('Probability of "seeing"', 'fontsize',16)
set(gca, 'fontsize',16)


%% standard threshold theory
x=-5:0.01:5;
psychoC=zeros(1, length(x));
th=1.0*randn(1,10);
figure; hold

ind1=find(x>0);
psychoC(ind1)=1;
plot(x, psychoC, 'r', 'LineWidth', 4)
set(gca,'YLim', [-0.01, 1.01], 'fontsize',16, 'LineWidth', 1.5);
xlabel('Stimulus relative to average threshold', 'fontsize',16)
ylabel('Probability of "seeing"', 'fontsize',16)
title('Psychometric function', 'fontsize',16)

for i=1:length(th)
    
    psychoC(1:end)=0;
    ind1=find(x>th(i));
    psychoC(ind1)=1;
    plot(x, psychoC, 'k', 'LineWidth', 1.5)
end

psychoC=normcdf(x,0,1);
plot(x, psychoC, 'g', 'LineWidth',4);

%% shot noise experiment wo noise (Hecht et al)
x=.1:0.1:20;
n=0:200;
% draw threshold at 10


%psychoC=zeros(1, length(x));
%th=1.0*randn(1,10);
figure('Position', [600 200 300 300]); 
hb=gca;
set(gca,'XLim', [min(x) max(x)], 'YLim', [0 1], 'fontsize',16, 'LineWidth', 1.5);
xlabel('mean # photons absorbed', 'fontsize',16)
ylabel('Probability seeing', 'fontsize',16)


figure('Position', [200 200 300 300]); 
ha=gca;
set(gca,'XLim', [0, 20], 'fontsize',16, 'LineWidth', 1.5);
xlabel('# photons absorbed', 'fontsize',16)
ylabel('Probability, P(n | <n>) ', 'fontsize',16)
line([10 10], [0 1], 'Color', 'r', 'LineWidth', 3) % threshold line

% draw threshold at 10 - mean at 10
nMean=10;
lm=line([nMean nMean], [0 1], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5)
y=poisspdf(n,nMean);
ld=line(n, y, 'LineWidth', 1.5)
ylabel(strcat('Probability, P(n |', num2str(nMean), ') '), 'fontsize',16)

axes(hb)
line(nMean, sum(y(10:end)), 'Color', 'g', 'Marker', 'o', 'MarkerFaceColor', 'g', 'MarkerSize', 10)

% draw threshold at 10 - mean at 8
axes(ha)
delete([lm ld])
nMean=8;
lm=line([nMean nMean], [0 1], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5)
y=poisspdf(n,nMean);
ld=line(n, y, 'LineWidth', 1.5)
ylabel(strcat('Probability, P(n |', num2str(nMean), ') '), 'fontsize',16)

axes(hb)
line(nMean, sum(y(10:end)), 'Color', 'g', 'Marker', 'o', 'MarkerFaceColor', 'g', 'MarkerSize', 10)

% draw threshold at 10 - mean at 12
axes(ha)
delete([lm ld])
nMean=12;
lm=line([nMean nMean], [0 1], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5)
y=poisspdf(n,nMean);
ld=line(n, y, 'LineWidth', 1.5)
ylabel(strcat('Probability, P(n |', num2str(nMean), ') '), 'fontsize',16)

axes(hb)
line(nMean, sum(y(10:end)), 'Color', 'g', 'Marker', 'o', 'MarkerFaceColor', 'g', 'MarkerSize', 10)

cla
pC=[];
for i=1:length(x)
    y=poisspdf(n,x(i));
    pC=[pC sum(y(10:end))];
end
line(x, pC, 'Color', 'g', 'LineWidth', 2)

%% Shot noise experiment with noise (Barlow)

x=.1:0.1:20;
n=0:200;
theta=10;

figure('Position', [600 200 300 300]); 
hb=gca;
set(gca,'XLim', [min(x) max(x)], 'YLim', [0 1], 'fontsize',16, 'LineWidth', 1.5);
xlabel('mean # photons absorbed', 'fontsize',16)
ylabel('Probability seeing', 'fontsize',16)

figure('Position', [200 200 300 300]); 
ha=gca;
set(gca,'XLim', [0, 20], 'fontsize',16, 'LineWidth', 1.5);
xlabel('# photons absorbed', 'fontsize',16)
ylabel('Probability', 'fontsize',16)
% line([theta theta], [0 1], 'Color', 'r', 'LineWidth', 3) % threshold line

% draw threshold at 10 - mean at 10
nMean=10;
nDark=1;
%lm=line([nMean nMean], [0 1], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5)
ym=poisspdf(n,nMean);
lm=line(n, ym, 'LineWidth', 1.5)
yd=poisspdf(n,nDark);
ld=line(n, yd, 'LineWidth', 1.5, 'Color', 'k')
ylabel(strcat('Probability'), 'fontsize',16)

% calculate probability of seeing for differnt thresholds

pC=zeros(1,3);

for k=1:3
    theta=4+2*k-1 % need to have count larger than this to detect photons
    for i=1:200
        for j=1:200
            if (i+j-2) > theta % theta is already detected; i,j =1 corresponds to p(n=)
                pC(k)=pC(k)+yd(i)*ym(j);
            end
        end
    end
end

axes(hb);
line([10 10 10], pC, 'Color', 'g', 'Marker', 'o', 'MarkerFaceColor', 'g', 'MarkerSize', 10, 'LineStyle', 'none')
line([10], sum(ym(11:end)), 'Color', 'k', 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 10, 'LineStyle', 'none')



for k=1:3
    theta=4+2*k-1 % need to have count larger than this to detect photons
    for i=1:200
        for j=1:200
            if (i+j-2) > theta % theta is already detected; i,j =1 corresponds to p(n=)
                pC(k)=pC(k)+yd(i)*ym(j);
            end
        end
    end
end

% Compute psychometric curve for theta=9; 5
numPts=length(x);
pC9=zeros(1, length(x));
pC5=zeros(1, length(x));
theta=9;
for k=1:numPts
    ym=poisspdf(n,x(k));
    for i=1:200
        for j=1:200
            if (i+j-2) > 9 
                pC9(k)=pC9(k)+yd(i)*ym(j);
            end
            if (i+j-2) > 5 
                pC5(k)=pC5(k)+yd(i)*ym(j);
            end
        end
    end
end

line(x, pC9, 'Color', 'g','LineWidth', 2)
line(x, pC5, 'Color', 'g','LineWidth', 2)
            
% Compute psychometric curve for theta=5





% simulate for threshold 2-10; nMean 0-20; claculate 

psychoC=zeros(9, 21);
simNum=10000;
for k=1:21
    countsM=poissrnd(k-1, [1 simNum]);
    countsD=poissrnd(nDark, [1 simNum]);
    countsTotal=countsD+countsM;
    for threshold=2:10
        ind=find(threshold<countsTotal);
        pychoC(threshold, k)=length(ind)/simNum;
    end
end



