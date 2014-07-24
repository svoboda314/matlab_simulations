

ls=50:10:500;
depth=[];
for ls=50:10:500
    z = fzero(@(z) twoPdepth(z,ls), 5*ls)
    depth=[depth, z];
end

plot(50:10:500, depth)


%%

ls=50;
z=50:1000;
y=zeros(1, length(z));
for i=1:length(z)
    y(i)=twoPdepth(z(i), ls);
end
plot(z, y)
min(y)
