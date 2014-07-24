
figure; hold

nTrials=100;
stimPerTrial=7;

% compute ISIs
meanISI=0.001;
minISI=0.030;
isi = random('exp',meanISI, 1,1000);
isi(find(isi < minISI))=minISI;
%hist(isi, 0:0.01:10)

stimTimes=zeros(nTrials, stimPerTrial);
for i=1:nTrials
    stimTimes(i, :)=cumsum(isi(((i-1)*stimPerTrial+1):i*stimPerTrial))-isi((i-1)*stimPerTrial+1);
    plot(stimTimes(i, :), i*ones(1, stimPerTrial), 'o');
end    
