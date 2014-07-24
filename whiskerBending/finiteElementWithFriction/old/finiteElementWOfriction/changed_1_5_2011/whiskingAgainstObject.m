

%% bendingRods5 with rotation
% whisker angle at based - theta0
% object position - (x y)
% rotate frame of reference by -theta0 by applying a rotation matrix
%
% 
% --------------misc simulation paramters
clear all; close all

n_segments=1000; % number of elements in the whisker
pole_pos_in_mm=8; % distance of pole form the face; along x- axis <-- CHANGE

numPosForces=19; % number of simulations with non-zero forces
incr=1.0; % increment of intended angle at base, degrees
thetaMin=incr; % min angle aqafter contact (i.e. non-zero force)
thetaMax=numPosForces*incr; % max angle with positive force

anglesAtBases=incr:incr:thetaMax; % intended angles at base
numZeroForces=15; % number of simulations with zero forces (i.e. before whisker strikes object)

% --------------setting up plots
r0=[pole_pos_in_mm/1000 0]'; % pole pos in meters (x, y)
figure(21); cla; hold
plot(r0(1)*1000, r0(2)*1000, 'kx', 'MarkerSize', 10); % point of force application
axis equal
title('Whisker contours with follicle at (0,0)');
xlabel('Distance (mm)');
ylabel('Distance (mm)');

% --------------mouse pad elasticitiy
% assume that moment = kTheta * dTheta 
% moment = [N m]; dTheta = [rad]; kTheta = [N m / rad]
kTheta = 5e-6; % in rad/(Nm); larger means stiffer <-- CHANGE

% --------------setting whisker parameters 
% rat whiskers
% wh.ym=3.5e9; % young's modulus in pascal
% wh.L=0.060; % whisker length in meters
% base_radius=1.00e-4;

% mouse whisker parameters
conical=1; % if 1 then conical whisker with linear taper
wh.ym=5e9; % young's modulus in pascal
base_radius= 33.5e-6; % in meters
wh.L = 16e-3; % whisker lenght in meters
wh.friction = 0.0;

segment=wh.L/n_segments; % in meters; total length is n_segments * segment

if ~ conical
    I(1:n_segments)=(pi/4)*base_radius^4; %uniform
else
    I(1:n_segments)=(pi/4)*(base_radius*(1-(1:n_segments)*segment/wh.L)).^4;  %conical
end
wh.I=I; % 2. moment of intertia

 % --------------simulation 
% numPosForces=20;
% incr=1.0;
% thetaMin=incr;
% thetaMax=numPosForces*incr;
% 
% anglesAtBases=incr:incr:thetaMax; 
% numZeroForces=15;
%C anglesAtBases=(-(numZeroForces-1))*incr:incr:thetaMax; 

tolerance=0.0000001; % this is a weird parameter - if smaller the fminsearch doesn't work; if larger, results are noisier as expected
options = optimset('TolFun',tolerance);

contourLoc=0.008; % in meters
f=1/30000;
initialCondition = [contourLoc f];

f_axial=[];
moment=[];

contours=zeros(length(anglesAtBases)+numZeroForces, 2, n_segments+1); % to store trajectories

for j=1:numPosForces
    
    theta0=pi*anglesAtBases(j)/180;
    rotTheta0=[cos(-theta0) -sin(-theta0); sin(-theta0) cos(-theta0)];
    r1=rotTheta0*r0;
    rotTheta0=[cos(-theta0) -sin(-theta0); sin(-theta0) cos(-theta0)];
    a=r1';
    
    [x,fval,exitflag,output] = fminsearch(@(x) bendingRods5(x, wh, a),initialCondition, options); % x = [s_pole force]
    
    if fval < tolerance*100 % Kluge; I don't understand why fminsearch 'converges' even though the 'tolerance' is never reached
        [temp, wh]=bendingRods5(x, wh);
        for i=1:size(wh.x, 2)
            temp=[wh.x(i) wh.y(i)]';
            temp=rotTheta0*temp;
            wh.x(i)=temp(1);
            wh.y(i)=temp(2);
        end
        %plot(wh.x*1000, -wh.y*1000, 'r', 'LineWidth', 1.0); % minus is to keep with convention
        contours(j+numZeroForces, 1, :)=wh.x*1000;
        contours(j+numZeroForces, 2, :)=wh.y*1000;
        plot(squeeze(contours(j+numZeroForces, 1, :)), -squeeze(contours(j+numZeroForces, 2, :)), 'r', 'LineWidth', 1.0);
    else
        wh.f_axial=0;
        wh.moment=0;
    end
    f_axial=[f_axial wh.f_axial];
    moment=[moment wh.moment];
    
end

%%
%intendedAnglesAtBases=anglesAtBases-(180/pi)*moment/kTheta; % correct for compliance in the whisker pad
intendedAnglesAtBases=anglesAtBases; % torn off correction for compliance in the whisker pad

%hold off
%figure; hold
ind=find(f_axial==0);
if ~isempty(ind)
    startind=min(ind)-1;
    for i=1:length(ind)
        intendedAnglesAtBases(startind+i)=intendedAnglesAtBases(startind+i-1)+incr;
        contours(startind+i+numZeroForces, 1, :)=segment*((1:(n_segments+1))-1)*cosd(intendedAnglesAtBases(startind+i))*1000;
        contours(startind+i+numZeroForces, 2, :)=-segment*((1:(n_segments+1))-1)*sind(intendedAnglesAtBases(startind+i))*1000;
        plot(squeeze(contours(startind+i+numZeroForces, 1, :)), -squeeze(contours(startind+i+numZeroForces, 2, :)), 'b', 'LineWidth', 1.0);
        %startind+i
    end
end

moment=abs(moment);
f_axial=abs(f_axial);
moment=[zeros(1, numZeroForces) moment];
f_axial=[zeros(1, numZeroForces)  f_axial];

intendedAnglesAtBases=[(-(numZeroForces-1)*incr):incr:0 intendedAnglesAtBases];

for i=1:numZeroForces
    contours(i, 1, :)=segment*((1:(n_segments+1))-1)*cosd(intendedAnglesAtBases(i))*1000;
    contours(i, 2, :)=-segment*((1:(n_segments+1))-1)*sind(intendedAnglesAtBases(i))*1000;
    plot(squeeze(contours(i, 1, :)), -squeeze(contours(i, 2, :)), 'b', 'LineWidth', 1.0);
end

hold off

figure(22); 
subplot(2,2,1); cla
xlabel('whisker angle at follicle (degrees)')
ylabel('moment (Nm)')
line(intendedAnglesAtBases, moment)

subplot(2,2,2); cla
xlabel('whisker angle at follicle (degrees)')
ylabel('axial_force (N)')
line(intendedAnglesAtBases, f_axial)

subplot(2, 2, 3); cla
ylabel('axial_force (N)')
xlabel('moment (Nm)')
line(moment, f_axial, 'LineStyle','none', 'Marker', '.')

simResults.f_axial=f_axial;
simResults.moment=moment;
simResults.contours=contours;
simResults.intendedAnglesAtBases=intendedAnglesAtBases;

%% ______________________animate whisking - two locations
% assume that contours have been tracked
% contours, f_axial, moment all need to be defined
% to do: 

load('simResults5');
f_axialGO=simResults5.f_axial;
contoursGO=simResults5.contours;
momentGO=simResults5.moment;

load('simResults8');
f_axialNG=simResults8.f_axial;
contoursNG=simResults8.contours;
momentNG=simResults8.moment;
intendedAnglesAtBases=simResults8.intendedAnglesAtBases

% load('intendedAnglesAtBases');
% load('contoursGO');load('contoursNG');
% load('momentGO');load('momentNG');
% load('f_axialGO');load('f_axialNG');

h23=figure(23); % arbitrary figure number
set(h23, 'Position', [50 50 1000 500], 'Color', 'w')

% axes on the left
x_lim=[-wh.L*1000*0.1 wh.L*1000*1.1]; y_lim=[-wh.L*1000*0.60 wh.L*1000*0.45];
hl=axes('Units', 'normalized', 'Position', [.07 .20 .25 .5], 'XLim', x_lim, 'YLim', ...
    y_lim, 'NextPlot', 'replacechildren', 'FontSize', 16);
axis square
xlabel('Distance (mm)'); ylabel('Distance (mm)')
ht=title('Go - object near', 'Color','b', 'FontSize', 20)

% plot circle to indicate pole
%h_pNG=line(8,0, 'Marker', 'o', 'Color', 'r', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
h_p=line(4,0, 'Marker', 'o', 'Color', 'b', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');

h_w=line(squeeze(contoursGO(1, 1, :)), -squeeze(contoursGO(1, 2, :)), 'Linewidth', 1.5, 'Color', 'k');
%define and plot the follicle
follX=[contours(1, 1, 1), -50*(contours(1, 1, 2)- contours(1, 1, 1))];
follY=[contours(1, 2, 1), 50*(contours(1, 2, 2)- contours(1, 2, 1))];
h_f=line(follX, follY, 'Color', 'm', 'Linewidth', 3);

% axes in the middle
xm=abs(max(intendedAnglesAtBases)-min(intendedAnglesAtBases)); x_lim=[min(intendedAnglesAtBases)-0.05*xm max(intendedAnglesAtBases)+0.05*xm]; 
ym=abs(max(momentGO)-min(momentGO)); y_lim=[min(momentGO)-0.05*ym max(momentGO)+0.05*ym];
hm=axes('Units', 'normalized', 'Position', [.38 .15 .25 .6], 'NextPlot', 'replacechildren', ...
    'XLim', x_lim,'YLim', y_lim, 'FontSize', 16);
xlabel('Angle at follicle (degrees)'); ylabel('Moment (Nm)')

% plot the moment vs angle line with green dot on top
h_mGO=line(intendedAnglesAtBases, momentGO, 'Linewidth', 1.5, 'Color', 'b');
h_mNG=line(intendedAnglesAtBases, momentNG, 'Linewidth', 1.5, 'Color', 'r');
h_mp=line(intendedAnglesAtBases(1), momentGO(1), 'Marker', 'o', 'Color', 'b', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'MarkerSize', 10); % green point - to move

% axes on the right
yr=abs(max(f_axialGO)-min(f_axialGO)); y_lim=[min(f_axialGO)-0.05*yr max(f_axialGO)+0.05*yr];
hr=axes('Units', 'normalized', 'Position', [.70 .15 .25 .6], 'NextPlot', 'replacechildren', ...
    'XLim', x_lim,'YLim', y_lim, 'FontSize', 16);
xlabel('Angle at follicle (degrees)'); ylabel('Axial Force (N)')


% plot the moment vs angle line with green dot on top
h_rGO=line(intendedAnglesAtBases, f_axialGO, 'Linewidth', 1.5, 'Color', 'b');
h_rNG=line(intendedAnglesAtBases, f_axialNG, 'Linewidth', 1.5, 'Color', 'r');
h_rp=line(intendedAnglesAtBases(1), f_axialGO(1), 'Marker', 'o', 'Color', 'b', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'MarkerSize', 10); % green point - to move

aviobj = avifile('testNoCompr1.avi', 'fps', 5, 'compression', 'None', 'quality', 100);

frameNumber=0;
for i=1:size(contoursGO,1)
    frameNumber=frameNumber+1;
    set(h_w, 'XData', squeeze(contoursGO(i, 1, :)), 'YData', -squeeze(contoursGO(i, 2, :))) 
    follX=[contoursGO(i, 1, 1), -50*(contoursGO(i, 1, 2)- contoursGO(i, 1, 1))];
    follY=[contoursGO(i, 2, 1), 50*(contoursGO(i, 2, 2)- contoursGO(i, 2, 1))];
    set(h_f, 'XData', follX, 'YData', follY);
    set(h_mp, 'XData', intendedAnglesAtBases(i),'YData',momentGO(i));
    set(h_rp, 'XData', intendedAnglesAtBases(i),'YData',f_axialGO(i));
    F(frameNumber)=getframe(h23);
%    aviobj=addframe(aviobj,F(frameNumber));
    %pause
end

for i=size(contoursGO,1):-1:1
    frameNumber=frameNumber+1;
    set(h_w, 'XData', squeeze(contoursGO(i, 1, :)), 'YData', -squeeze(contoursGO(i, 2, :))) 
    follX=[contoursGO(i, 1, 1), -50*(contoursGO(i, 1, 2)- contoursGO(i, 1, 1))];
    follY=[contoursGO(i, 2, 1), 50*(contoursGO(i, 2, 2)- contoursGO(i, 2, 1))];
    set(h_f, 'XData', follX, 'YData', follY);
    set(h_mp, 'XData', intendedAnglesAtBases(i),'YData',momentGO(i));
    set(h_rp, 'XData', intendedAnglesAtBases(i),'YData',f_axialGO(i));
    F(frameNumber)=getframe(h23);
%    aviobj=addframe(aviobj,F(frameNumber));
    %pause
end

% set up NG conditions
% remove green dot; add red dot

axes(hl);
set(h_p, 'XData', 8, 'YData', 0, 'Color', 'r', 'MarkerFaceColor', 'r');
%h_p=line(8,0, 'Marker', 'o', 'Color', 'r', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
delete(ht);
ht=title('NoGo - object far', 'Color','r', 'FontSize', 20)

axes(hm);
set(h_mp, 'XData', intendedAnglesAtBases(1), 'YData', momentNG(1), 'Color', 'r', 'MarkerFaceColor', 'r');
%h_mp=line(intendedAnglesAtBases(1), momentNG(1), 'Marker', 'o', 'Color', 'r', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'MarkerSize', 10); 

axes(hr);
set(h_rp, 'XData', intendedAnglesAtBases(1), 'YData', f_axialNG(1), 'Color', 'r', 'MarkerFaceColor', 'r');
%h_rp=line(intendedAnglesAtBases(1), f_axialNG(1), 'Marker', 'o', 'Color', 'r', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'MarkerSize', 10); 

for i=1:size(contoursNG,1)
    frameNumber=frameNumber+1;
    set(h_w, 'XData', squeeze(contoursNG(i, 1, :)), 'YData', -squeeze(contoursNG(i, 2, :))) 
    follX=[contoursNG(i, 1, 1), -50*(contoursNG(i, 1, 2)- contoursNG(i, 1, 1))];
    follY=[contoursNG(i, 2, 1), 50*(contoursNG(i, 2, 2)- contoursNG(i, 2, 1))];
    set(h_f, 'XData', follX, 'YData', follY);
    set(h_mp, 'XData', intendedAnglesAtBases(i),'YData',momentNG(i));
    set(h_rp, 'XData', intendedAnglesAtBases(i),'YData',f_axialNG(i));
    F(frameNumber)=getframe(h23);
%    aviobj=addframe(aviobj,F(frameNumber));
    %pause
end

for i=size(contoursNG,1):-1:1
    frameNumber=frameNumber+1;
    set(h_w, 'XData', squeeze(contoursNG(i, 1, :)), 'YData', -squeeze(contoursNG(i, 2, :))) 
    follX=[contoursNG(i, 1, 1), -50*(contoursNG(i, 1, 2)- contoursNG(i, 1, 1))];
    follY=[contoursNG(i, 2, 1), 50*(contoursNG(i, 2, 2)- contoursNG(i, 2, 1))];
    set(h_f, 'XData', follX, 'YData', follY);
    set(h_mp, 'XData', intendedAnglesAtBases(i),'YData',momentNG(i));
    set(h_rp, 'XData', intendedAnglesAtBases(i),'YData',f_axialNG(i));
    F(frameNumber)=getframe(h23);
%    aviobj=addframe(aviobj,F(frameNumber));
    %pause
end

aviobj=close(aviobj);
movie(h23, F, 3)
% movie2avi(F, 'testRLE1.avi', 'fps', 5, 'compression', 'RLE', 'quality', 100) % no MSVC, RLE, 

%% animate whisking - one location
% assume that contours have been tracked
% contours, f_axial, moment all need to be defined
% to do: 

h23=figure(23); % arbitrary figure number
set(h23, 'Position', [50 50 1000 500], 'Color', 'w')

% axes on the left
x_lim=[-wh.L*1000*0.1 wh.L*1000*1.1]; y_lim=[-wh.L*1000*0.60 wh.L*1000*0.45];
hl=axes('Units', 'normalized', 'Position', [.05 .20 .25 .5], 'XLim', x_lim, 'YLim', ...
    y_lim, 'NextPlot', 'replacechildren');
axis square
xlabel('Distance (mm)'); ylabel('Distance (mm)')

% plot circle to indicate pole
h_p=line(r0(1)*1000,r0(2)*1000+0.3, 'Marker', 'o', 'Color', 'r', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');

h_w=line(squeeze(contours(1, 1, :)), -squeeze(contours(1, 2, :)), 'Linewidth', 1.5);
%define and plot the follicle
follX=[contours(1, 1, 1), -50*(contours(1, 1, 2)- contours(1, 1, 1))];
follY=[contours(1, 2, 1), 50*(contours(1, 2, 2)- contours(1, 2, 1))];
h_f=line(follX, follY, 'Color', 'k', 'Linewidth', 3);

% axes in the middle
xm=abs(max(intendedAnglesAtBases)-min(intendedAnglesAtBases)); x_lim=[min(intendedAnglesAtBases)-0.05*xm max(intendedAnglesAtBases)+0.05*xm]; 
ym=abs(max(moment)-min(moment)); y_lim=[min(moment)-0.05*ym max(moment)+0.05*ym];
hm=axes('Units', 'normalized', 'Position', [.4 .15 .25 .6], 'NextPlot', 'replacechildren', ...
    'XLim', x_lim,'YLim', y_lim);
xlabel('Angle at follicle (degrees)'); ylabel('Moment (Nm)')

% plot the moment vs angle line with green dot on top
h_m=line(intendedAnglesAtBases, moment, 'Linewidth', 1.5);
h_mp=line(intendedAnglesAtBases(1), moment(1), 'Marker', 'o', 'Color', 'b', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'MarkerSize', 10); % green point - to move

% axes on the right
yr=abs(max(f_axial)-min(f_axial)); y_lim=[min(f_axial)-0.05*yr max(f_axial)+0.05*yr];
hr=axes('Units', 'normalized', 'Position', [.70 .15 .25 .6], 'NextPlot', 'replacechildren', ...
    'XLim', x_lim,'YLim', y_lim);
xlabel('Angle at follicle (degrees)'); ylabel('Axial Force (N)')

% plot the moment vs angle line with green dot on top
h_r=line(intendedAnglesAtBases, f_axial, 'Linewidth', 1.5);
h_rp=line(intendedAnglesAtBases(1), f_axial(1), 'Marker', 'o', 'Color', 'b', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'MarkerSize', 10); % green point - to move

aviobj = avifile('testCinepak2.avi', 'fps', 5, 'compression', 'Cinepak', 'quality', 100);

for i=1:size(contours,1)
    set(h_w, 'XData', squeeze(contours(i, 1, :)), 'YData', -squeeze(contours(i, 2, :))) 
    follX=[contours(i, 1, 1), -50*(contours(i, 1, 2)- contours(i, 1, 1))];
    follY=[contours(i, 2, 1), 50*(contours(i, 2, 2)- contours(i, 2, 1))];
    set(h_f, 'XData', follX, 'YData', follY);
    set(h_mp, 'XData', intendedAnglesAtBases(i),'YData',moment(i));
    set(h_rp, 'XData', intendedAnglesAtBases(i),'YData',f_axial(i));
    F(i)=getframe(h23);
    aviobj=addframe(aviobj,F(i));
    pause(0.2)
end

for i=size(contours,1):-1:1
    set(h_w, 'XData', squeeze(contours(i, 1, :)), 'YData', -squeeze(contours(i, 2, :))) 
    follX=[contours(i, 1, 1), -50*(contours(i, 1, 2)- contours(i, 1, 1))];
    follY=[contours(i, 2, 1), 50*(contours(i, 2, 2)- contours(i, 2, 1))];
    set(h_f, 'XData', follX, 'YData', follY);
    set(h_mp, 'XData', intendedAnglesAtBases(i),'YData',moment(i));
    set(h_rp, 'XData', intendedAnglesAtBases(i),'YData',f_axial(i));
    F(2*size(contours,1)-i+1)=getframe(h23);
    aviobj=addframe(aviobj,F(i));
    pause(.2)
end

aviobj=close(aviobj);
%movie(h23, F, 3)
%movie2avi(F, 'testRLE1.avi', 'fps', 5, 'compression', 'RLE', 'quality',
%100) % no MSVC, RLE, 