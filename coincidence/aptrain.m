
rate=5; % Hz
time=100000; % in ms
synStrength=0.1;
synDelay=2;

temp=rand(1,time);
presyn=zeros(1, time);
indpre=find(temp < (rate/ 1000));
presyn(indpre)=1;

temp=rand(1,time);
postsyn=zeros(1, time);
ind=find(temp < (rate/ 1000));
postsyn(ind)=1;

% for eaxh presynaptic AP -- a prob of generating an AP
temp=rand(1,length(indpre));
postsyn2=zeros(1, length(indpre));
postsyn2(temp < synStrength)=1;

indpre=indpre+2;
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
plot(c(time:(time+timewindow)));