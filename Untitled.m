
t=0:1000;

st=[100 120 130 300 400 410];

g=zeros(1, length(t));

tau=10;

for i=1:length(st)
    tt=st(i):1000;
    g((st(i)+1):end)=g((st(i)+1):end)+(t((st(i)+1):end)/tau).*exp(1-t((st(i)+1):end)/tau);
end