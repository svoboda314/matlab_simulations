
x=rand(40,1);
y=x+2*rand(40,1)-1;
plot(x, y, 'o');
p = polyfit(x,y,1)
%%

an=10;
angles=10;

slps=zeros(an, angles);

for i=1:an
    
    for j=1:angles
        x=rand(40,1);
    y=x+2*rand(40,1)-1;
    p = polyfit(x,y,1);
    slps(i, j)=p(1);
    end
end
    
for i=1:an
    [temp m_ind]= max(slps(i, :));
    slps(i, 1:(angles-m_ind+1))=slps(i, m_ind:angles);
    slps(i, (angles-m_ind+2):angles)=nan;
    plot(nanmean(slps)) 
end