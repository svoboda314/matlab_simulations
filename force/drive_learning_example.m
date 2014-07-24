
% Setup an entire learning example. 

% Learn a simple sine wave in the presence of a lot of noise.
if ( 1 ) 
    N = 300;
    B = 1;
    I = 0;
    p = 0.5;
    g = 1.5;
    dt = 0.1;
    tau = 1.0;
    
    net1 = create_random_model(N, B, I, p, g, dt, tau, 'netnoisesigma', 0.0, 'actfuntype', 0);
    
    amp = 1.0;
    freq = 1.0/100.0;
    nsecs = 1.0/freq * 20.0;		% n iterations of the basic frequency
    simtime = 0:dt:nsecs-dt;
    ft = (amp/1.0)*sin(2.0*pi*freq*simtime);   
    
    [net1, learn_struct, simdata] = learn_model(0.1*randn(N,1), net1, ft, 'doplot', 1, 'learnevery', 10);   
    
    time = 2000.0;
    [Z, R, X] = run_model(0.1*randn(net1.N,1), time, net1);
    figure; 
    plot(Z');
    hold on;
    title('Example 1');
    
    net1_nfb = transfer_learning_internal(net1, learn_struct);
    time = 2000.0;
    [Znfb, Rnfb, Xnfb] = run_model(0.1*randn(N,1), time, net1_nfb);
    figure; 
    plot(Znfb', 'r')
    title('Example 1 (Internal learning)');    
    
end

% Here's a more sophisticated example that uses inputs to determine the output (again sine waves).  In addition, to
% flex the code, there are two outputs now, each feeding back to the network.
if ( 1 )
    N = 300;
    B = 2;
    I = 2;
    p = 0.5;
    g = 1.5;
    dt = 0.1;
    tau = 1.0;
    
    net2 = create_random_model(N, B, I, p, g, dt, tau, 'actfuntype', 0);
    
    amp = 1.0;
    freq = 1.0/100.0;
    nsecs = 1.0/freq * 10.0;		% n iterations of the basic frequency
    simtime = 0:dt:nsecs-dt;
    ft11 = (amp/1.0)*sin(2.0*pi*freq*simtime);
    ft12 = (amp/1.0)*cos(2.0*pi*freq*simtime);
    ft1 = [ft11; ft12];
    ft21 = (amp/1.0)*sin(3.0*pi*freq*simtime);
    ft22 = (amp/1.0)*cos(3.0*pi*freq*simtime);
    ft2 = [ft21; ft22];
    inp1 = ones(I,length(simtime));
    inp1(1,:) = 0;
    inp2 = ones(I,length(simtime));
    inp2(2,:) = 0;
    
    learn_every = 5;
    
    time = 100.0;			% Just get the network running a little bit. 
    [Z, R, X] = run_model(0.1*randn(net2.N,1), time, net2);
    sd.X = X;				% hack for the loop
    learn_struct = [];
    figure;
    dl_yes = [1;1];
    dl_no = [0;0];
    tf_yes = [1;1];
    tf_no = [1;1];
    for i = 1:3
	[net2, learn_struct, sd] = learn_model(sd.X(:,end), net2, ft1, 'input', inp1, ... % just used to run the network some
					      'learnstruct', learn_struct, 'learnevery', learn_every, 'DL', dl_no, 'TF', tf_yes);
	
	[net2, learn_struct, sd] = learn_model(sd.X(:,end), net2, ft1, 'input', inp1, ...
					      'learnstruct', learn_struct, 'learnevery', learn_every, 'DL', dl_yes, 'TF', tf_no);
	plot(sd.Z'); hold on; pause(1.0);
	[net2, learn_struct, sd] = learn_model(sd.X(:,end), net2, ft2, 'input', inp2, ... % just used to run the network some
					      'learnstruct', learn_struct, 'learnevery', learn_every, 'DL', dl_no, 'TF', tf_yes);
	
	[net2, learn_struct, sd] = learn_model(sd.X(:,end), net2, ft2, 'input', inp2, ...
					      'learnstruct', learn_struct, 'learnevery', learn_every, 'DL', dl_yes, 'TF', tf_no);
	plot(sd.Z', 'r'); pause(1.0);
    end
    title('Example 2 Learning');
    
    % Note that the two sets of sine waves run at different times.
    time = 1000.0;
    [Z1, R1, X1] = run_model(0.1*randn(net2.N,1), time, net2, 'input', inp1);
    [Z2, R2, X2] = run_model(0.1*randn(net2.N,1), time, net2, 'input', inp2);
    
    figure; 
    plot(Z1', 'b'); 
    hold on; 
    plot(Z2', 'r');
    title('Example 2');
       
    
    net2_nfb = transfer_learning_internal(net2, learn_struct);
    time = 1000.0;
    [Z1nfb, R1nfb, X1nfb] = run_model(0.1*randn(net2.N,1), time, net2_nfb, 'input', inp1);
    [Z2nfb, R2nfb, X2nfb] = run_model(0.1*randn(net2.N,1), time, net2_nfb, 'input', inp2);
    figure; 
    plot(Z1', 'b'); 
    hold on; 
    plot(Z2', 'r');
    title('Example 2 (Internal learning)');
    
    
end    
    
