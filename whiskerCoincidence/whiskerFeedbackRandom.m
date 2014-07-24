
nIterations=10000;

hit=zeros(nIterations, 1);
threshold = 11;

for i=1:nIterations

dominant=rand(1,22);
principal=rand(1,22);

dominant(find(dominant < (1/3)))=1;
dominant(find(dominant < (2/3)))=2;
dominant(find(dominant < 1))=3;

principal(find(principal < (1/3)))=1;
principal(find(principal < (2/3)))=2;
principal(find(principal < 1))=3;

pd=abs(dominant - principal);

ind=find(pd==0);
if length(ind)> threshold
hit(i)=1;
end
end

p=sum(hit)/nIterations
