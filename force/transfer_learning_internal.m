function new_net = transfer_learning_internal(net, learn_struct)

I = net.I;
B = net.B;
N = net.N;

P = learn_struct.P;
P = (P + P')/2.0;			% Make certain that P really is symmetric.
C = inv(P);				% Since P is diagonally loaded, this inverse should go through.

wfb = net.wFB;
wout = net.wOut;
Cw = C*wout;
g = net.g;
J = net.J;			% Have to subsume into J, otherwise the delta gets multiplied by g.
newJ = zeros(N,N);
parfor n = 1:N				% If you don't have the latest matlab, you should do 'for'.
    idxs = find(J(n,:));
    SCS = C(idxs,idxs);
    SCw = Cw(idxs,:);
    
    if ( B == 1 )			% separate for clarity
	delta_j{n} = (wfb(n) * pinv(SCS) * SCw)';  % row 	
    else
	delta_j_nwfb = pinv(SCS) * SCw; 	   % columns
	delta_j_all = bsxfun(@times, delta_j_nwfb, wfb(n,:)); % rows
	delta_j{n} = sum(delta_j_nwfb'); 	   % row
    end

    jidxs{n} = idxs;
end
for n = 1:N
    idxs = jidxs{n};
    newJ(n,idxs) = J(n,idxs) + delta_j{n}/g;
end

new_net = net;
%new_net.g = 1.0;			% Have to subsume into J, otherwise the delta gets multiplied by g.
new_net.J = newJ;
new_net.wFB = zeros(N,B);














































