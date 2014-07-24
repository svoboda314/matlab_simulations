
t=1:.01:1000;
tf=100;
tr=10;
gm=10;

c=((tr/(tr+tf))^(tr/tf)-(tr/(tr+tf))^((tr+tf)/tf))^-1;

g_fnmda=gm*c*(exp(-t/tf)-exp(-t/(tr*tf/(tr+tf))));

%g=exp(-t/tf)-exp(-t/t2);

