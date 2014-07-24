th=0.05;
corr=0.02;
t=100000;

ci=[];
for i=1:100

x1=rand([t 1]);
ind=find(x1>th);
x1(ind)=1;
x1=fix(x1);
% plot(x1);

x2=rand([t 1]);
ind2=find(x2>th);
x2(ind2)=1;
x2=fix(x2);

xHi=rand([length(ind) 1]);
xHi(find(xHi>(th-corr)))=1;
xLo=rand([length(ind) 1]);
xLo(find(xLo>(th+corr)))=1;

x2(ind)=xHi;
x2(ind+1)=xLo;
x2=x2(1:t);

ci=[ci, (mean(x1.*x2)-mean(x1)*mean(x2))/sqrt(mean(x1)*mean(x2))];

end