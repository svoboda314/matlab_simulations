
trials=10;

bins=30;
n=1:bins;
nStim=5;
stim=zeros(trials,bins);
stimInd=ones(trials, nStim); % index of stimuli; range 1 - bins
stimPeriod=1; % in seconds
timeStep=stimPeriod/bins; % min time between stims, in seconds
time=n*timeStep; % 

% generate random isi's with min delay; use exponential distribution



for j=1:trials
    for i=1:nStim
        
    end
end













figure; hold
for j=1:trials
    for i=1:nStim
        temp=randi(30);
        while (stim(j, temp) == 1)
            temp=randi(30);
        end
        stimInd(j, :)=temp;
        plot(stimInd(j, :), j*ones(1, nStim), 'o')
        stim(j, temp)=1;
    end
end

plot(sum(stim), 'o');

figure; hold
stimInd=ones(trials, nStim);
meanIsi=1;
for j=1:trials
    isi=round(rand*(meanIsi)+1);
    stimInd(j, 1)=isi;
    for i=2:nStim
        isi=round(rand*(meanIsi)+1);
        stimInd(j, i)=stimInd(j, i-1)+isi;
    end
    plot(stimInd(j, :), j*ones(1, nStim), 'o')
        stim(j, temp)=1;
end


