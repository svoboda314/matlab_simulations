
rate=5; % Hz
time=100000; % in ms
synStrength=0.1;
synDelay=2; % delay of onset of impact in ms
synSpread=3; %duration of impact in ms

% generate presynaptic spike-train
temp=rand(1,time);
presyn=zeros(1, time);
indpre=find(temp < (rate/ 1000));
presyn(indpre)=1;

% generate unperturbed postsynaptic spike-train
temp=rand(1,time);
postsyn=zeros(1, time);
ind=find(temp < (rate/ 1000));
postsyn(ind)=1;

% for each presynaptic AP -- a prob of generating a postsynaptic AP
temp=rand(1,length(indpre));
postsyn2=zeros(1, length(indpre));
postsyn2(temp < synStrength)=1;

% jitter corresponding to synSpread
temp=rand(1,length(indpre));
ind0=find(temp<(1/3));
ind2=find(temp>(2/3));
ind1=setdiff(1:length(indpre), [ind0, ind2]);
temp(ind0)=0;
temp(ind1)=1;
temp(ind2)=2;

indpre=indpre+temp;
IndecesLargerThanTime=find(indpre > time);
num=length(IndecesLargerThanTime);
indpre=indpre(1:(end-num)); 
postsyn2=postsyn2(1:(end-num));

% if indpre(end) > time; 
%     indpre=indpre(1:(end-1)); 
%     postsyn2=postsyn2(1:(end-1));
% end

postsyn(indpre)=postsyn(indpre)+postsyn2;
c=xcorr(postsyn, presyn);
timewindow=25;
plot(c((time-timewindow):(time+timewindow)), '-*');
