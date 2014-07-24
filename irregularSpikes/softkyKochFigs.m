


%% n Poisson neuron onto one IF neuron - Fig 6
% for leaky integrator make tau = 10

n_neur = 10;

duration = 10;
tstep=1; %in ms
t=0:1:(duration*1000); %in ms
tsteps=length(t);
sp_rate=40; % in Hz
rate=ones(1, tsteps)*sp_rate;
duration=duration+0.001; %in seconds
epsp = 0.05; % in Vth units
tau = 1000; % in milliseconds

sp_times = cell(1,n_neur); 
for i=1:n_neur
temp=inhPoisson(rate, duration);
%temp=round(temp*1000);
sp_times{i}=temp; % in ms
end

vm=zeros(1, tsteps);
sp_out=[];
for i=2:tsteps
%    dv=vm(i-1)*(1-exp(-tstep/tau));
%    vm(i)=vm(i-1)-dv;
     vm(i)=vm(i-1)*exp(-tstep/tau);
    for j=1:n_neur
    if ismember(i-1, sp_times{j})
        vm(i)=vm(i)+epsp;
    end
    end
    if vm(i) > 1
    sp_out=[sp_out, i];
    vm(i)=0;
    end
end

plot(vm)
hold
plot(sp_out, ones(1, length(sp_out)), 'o')
hold off
figure;
hist(diff(sp_out), 5:10:500)

cv= std(diff(sp_out))/mean(diff(sp_out))


spike_rate = length(sp_out)/duration