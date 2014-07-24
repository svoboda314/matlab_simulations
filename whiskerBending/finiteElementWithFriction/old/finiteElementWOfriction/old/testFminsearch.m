%% test fminsearch from manual

options = optimset('TolFun',0.00001);
x=fminsearch(@funcToMinimize, [1 1 1], options);


%% bendingRods3
contourLoc=0.01;
f=1/300;
initialCondition = [contourLoc f];

x = fminsearch(@bendingRods3,[contourLoc f])

%% bendingRods4 - whisker paramters are defined in function
a=[0.01 0.005];
contourLoc=0.01;
f=1/300;
initialCondition = [contourLoc f];

x = fminsearch(@(x) bendingRods4(x, a),initialCondition);

%% bendingRods5 - whisker parameters are passed as parameter

% --------------setting whisker parameters 
conical=1;
% rat whiskers
wh.ym=3.5e9; % young's modulus in pascal
wh.L=0.060; % whisker length in meters
base_radius=1.00e-4;

n_segments=1000;
segment=wh.L/n_segments; % in meters; total length is n_segments * segment
 
if ~ conical
    I(1:n_segments)=(pi/4)*base_radius^4; %uniform
else
    I(1:n_segments)=(pi/4)*(base_radius*(1-(1:n_segments)*segment/wh.L)).^4;  %conical
end
wh.I=I;

% --------------setting simulation parameters 
a=[0.02 0.005];
contourLoc=0.01;
f=1/300;
initialCondition = [contourLoc f];

options = optimset('TolFun',0.0000001);
x = fminsearch(@(x) bendingRods5(x, a, wh),initialCondition);

wh=bendingRods5wPlot(x, wh);

figure(19); cla; hold
plot(wh.x*1000, wh.y*1000, 'b', 'LineWidth', 1.5);
plot(a(1)*1000, -a(2)*1000, 'rx'); % point of force application
axis equal
hold off

%% bendingRods5 with rotation
% whisker angle at based - theta0
% object position - (x y)
% rotate frame of reference by -theta0 by applying a rotation matrix
%
% to do
% animation
% 

% --------------setting whisker parameters 

% rat whiskers
% wh.ym=3.5e9; % young's modulus in pascal
% wh.L=0.060; % whisker length in meters
% base_radius=1.00e-4;

% mouse pad elasticitiy
% assume that moment = kTheta * dTheta 
% moment = [N m]; dTheta = [rad]; kTheta = [N m / rad]

kTheta = 1e-6; % in rad/(Nm); larger means stiffer

% mouse whisker parameters
conical=1;
wh.ym=5e9;
base_radius= 33.5e-6; % in meters
wh.L = 16e-3; % in meters

pole_pos_in_mm=5;

n_segments=1000;
segment=wh.L/n_segments; % in meters; total length is n_segments * segment

if ~ conical
    I(1:n_segments)=(pi/4)*base_radius^4; %uniform
else
    I(1:n_segments)=(pi/4)*(base_radius*(1-(1:n_segments)*segment/wh.L)).^4;  %conical
end
wh.I=I;

% --------------setting up plots
r0=[pole_pos_in_mm/1000 0]'; % pole pos in meters
figure(21); cla; hold
plot(r0(1)*1000, r0(2)*1000, 'rx'); % point of force application
axis equal

% --------------simulation 
numPosForces=20;
incr=2;
thetaMin=incr;
thetaMax=numPosForces*incr;

numZeroForces=5;
anglesAtBases=(-(numZeroForces-1))*incr:incr:thetaMax; 

tolerance=0.0000001;
options = optimset('TolFun',tolerance);

contourLoc=0.004;
f=1/30000;
initialCondition = [contourLoc f];

f_axial=[];
moment=[];

contours=zeros(length(anglesAtBases), 2, n_segments+1); % to store trajectories

for j=1:numPosForces
    
    theta0=pi*anglesAtBases(j+numZeroForces)/180;
    rotTheta0=[cos(-theta0) -sin(-theta0); sin(-theta0) cos(-theta0)];
    r1=rotTheta0*r0;
    rotTheta0=[cos(-theta0) -sin(-theta0); sin(-theta0) cos(-theta0)];
    a=r1';
    
    [x,fval,exitflag,output] = fminsearch(@(x) bendingRods5(x, a, wh),initialCondition, options); % x = [s_pole force]
    
    %fval
    %exitflag
    %output
    %x
    %pause
    
    if fval > tolerance*100 % Kluge 
        wh.phi(1:n_segments)=theta0;
        wh.x=zeros(1, n_segments+1);
        wh.y=zeros(1, n_segments+1);
        wh.x(2:end)=cumsum(segment*cos(wh.phi));
        wh.y(2:end)=-cumsum(segment*sin(wh.phi));
        wh.f_axial=0;
        wh.moment=0
    else
        wh=bendingRods5wPlot(x, wh);
        for i=1:size(wh.x, 2)
            temp=[wh.x(i) wh.y(i)]';
            temp=rotTheta0*temp;
            wh.x(i)=temp(1);
            wh.y(i)=temp(2);
        end
    end
    f_axial=[f_axial wh.f_axial];
    moment=[moment wh.moment];
    
    %   plot(wh.x*1000, -wh.y*1000, 'b', 'LineWidth', 1.0); % minus is to keep with convention
    %   plot(r1(1)*1000, r1(2)*1000, 'rx')
    
    plot(wh.x*1000, -wh.y*1000, 'r', 'LineWidth', 1.0); % minus is to keep with convention
    contours(j+numZeroForces, 1, :)=wh.x*1000;
    contours(j+numZeroForces, 2, :)=wh.y*1000;
end
hold off

moment=abs(moment);
f_axial=abs(f_axial);

% construct some undeflected contours

for i=1:numZeroForces
    contours(i, 1, :)=segment*((1:(n_segments+1))-1)*cosd(anglesAtBases(i))*1000;
    contours(i, 2, :)=-segment*((1:(n_segments+1))-1)*sind(anglesAtBases(i))*1000;
end

moment=[zeros(1,numZeroForces)  moment];
f_axial=[zeros(1,numZeroForces)  f_axial];

%intendedAnglesAtBases=anglesAtBases+(180/pi)*moment/kTheta;

figure(22); 
subplot(2,2,1); cla
xlabel('whisker angle at follicle (degrees)')
ylabel('moment (Nm)')
line(anglesAtBases, moment)

subplot(2,2,2); cla
xlabel('whisker angle at follicle (degrees)')
ylabel('axial_force (N)')
line(anglesAtBases, f_axial)

subplot(2, 2, 3); cla
ylabel('axial_force (N)')
xlabel('moment (Nm)')
line(moment, f_axial, 'LineStyle','none', 'Marker', '.')


%% animate whisking
% assume that contours have been tracked
% to do: 

h23=figure(23); % arbitrary figure number
set(h23, 'Position', [50 50 1000 500], 'Color', 'w')

% axes on the left
x_lim=[-wh.L*1000*0.1 wh.L*1000*1.1]; y_lim=[-wh.L*1000*0.80 wh.L*1000*0.25];
hl=axes('Units', 'normalized', 'Position', [.1 .15 .35 .75], 'XLim', x_lim, 'YLim', ...
    y_lim, 'NextPlot', 'replacechildren');
xlabel('Distance (mm)'); ylabel('Distance (mm)')

% plot circle to indicate pole
h_p=line(r0(1)*1000,r0(2)*1000+0.3, 'Marker', 'o', 'Color', 'r', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');

h_w=line(squeeze(contours(1, 1, :)), -squeeze(contours(1, 2, :)), 'Linewidth', 1.5);
%define and plot the follicle
follX=[contours(1, 1, 1), -50*(contours(1, 1, 2)- contours(1, 1, 1))];
follY=[contours(1, 2, 1), 50*(contours(1, 2, 2)- contours(1, 2, 1))];
h_f=line(follX, follY, 'Color', 'k', 'Linewidth', 3);

% axes on the right
xr=abs(max(anglesAtBases)-min(anglesAtBases)); x_lim=[min(anglesAtBases)-0.05*xr max(anglesAtBases)+0.05*xr]; 
yr=abs(max(moment)-min(moment)); y_lim=[min(moment)-0.05*yr max(moment)+0.05*yr];
hr=axes('Units', 'normalized', 'Position', [.6 .15 .35 .75], 'NextPlot', 'replacechildren', ...
    'XLim', x_lim,'YLim', y_lim);
xlabel('Angle at follicle (degrees)'); ylabel('Moment (Nm)')

% plot the moment vs angle line with green dot on top
h_m=line(anglesAtBases, moment, 'Linewidth', 1.5);
h_mp=line(anglesAtBases(1), moment(1), 'Marker', 'o', 'Color', 'g', 'LineWidth', 1.5, 'MarkerFaceColor', 'g', 'MarkerSize', 10); % green point - to move

aviobj = avifile('testCinepak2.avi', 'fps', 5, 'compression', 'Cinepak', 'quality', 100);

for i=1:size(contours,1)
    set(h_w, 'XData', squeeze(contours(i, 1, :)), 'YData', -squeeze(contours(i, 2, :))) 
    follX=[contours(i, 1, 1), -50*(contours(i, 1, 2)- contours(i, 1, 1))];
    follY=[contours(i, 2, 1), 50*(contours(i, 2, 2)- contours(i, 2, 1))];
    set(h_f, 'XData', follX, 'YData', follY);
    set(h_mp, 'XData', anglesAtBases(i),'YData',moment(i));
    F(i)=getframe(h23);
    aviobj=addframe(aviobj,F(i));
    pause(0.2)
end

for i=size(contours,1):-1:1
    set(h_w, 'XData', squeeze(contours(i, 1, :)), 'YData', -squeeze(contours(i, 2, :))) 
    follX=[contours(i, 1, 1), -50*(contours(i, 1, 2)- contours(i, 1, 1))];
    follY=[contours(i, 2, 1), 50*(contours(i, 2, 2)- contours(i, 2, 1))];
    set(h_f, 'XData', follX, 'YData', follY);
    set(h_mp, 'XData', anglesAtBases(i),'YData',moment(i));
    F(2*size(contours,1)-i+1)=getframe(h23);
    aviobj=addframe(aviobj,F(i));
    pause(.2)
end

aviobj=close(aviobj);
%movie(h23, F, 3)
%movie2avi(F, 'testRLE1.avi', 'fps', 5, 'compression', 'RLE', 'quality', 100) % no MSVC, RLE, 
