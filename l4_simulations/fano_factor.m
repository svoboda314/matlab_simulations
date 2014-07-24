

tlen=1000;
t=1:tlen;
sr=10*(1.5+sin(2*pi*t/1000));

nn=100;
min_tbin=10;

spike_count=zeros(nn, round(tlen/min_tbin));

tt=1:min_tbin:tlen;

for i=1:nn
    st=inhPoisson(sr);
    
    for j=1:length(tt)
        spike_count(i,j)= length(find(st > (tt(j)-1) & st < (tt(j)+10)));
    end
    
end


ff=zeros(1, length(tt));


for j=1:length(tt)
   i=1;
   while i < length(tt)-j
       temp=sum(spike_count(:, i:(i+j-1)), 2);
       ff(j)=ff(j)+std(temp)^2/mean(temp);
       i=i+j;
   end
   ff(j)=ff(j)/((i-1)/j); 
    
end



% plot(st, ones(1, length(st)), 'o');
% hold
% plot(tt, spike_count(1, :), 'x');
% hold off