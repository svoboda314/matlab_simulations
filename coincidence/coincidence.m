
f1=rand(100000,1);
f2=rand(100000,1);

t1=0.8;
t2=0.8;

f1(find(f1>t1))=1;
f2(find(f2>t2))=1;
f1=fix(f1);
%sum(f1)
f2=fix(f2);
%sum(f2)
plot(f1, 'ro');
hold
plot(f2, 'go');

co=f1.*f2;
sum(co)
plot(co, 'bx');

