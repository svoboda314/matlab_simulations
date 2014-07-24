
tr=3;
td=50;

vv=-80:10:0;
tt=0:.1:500;

ginf=nmdaV(vv);

epsc=zeros(length(tt), length(vv));

%temp=(exp(-tt/td)-exp(-tt/tr))';

%temp1=temp*ginf;

epsc=(exp(-tt/td)-exp(-tt/tr))'*ginf;

plot(epsc)