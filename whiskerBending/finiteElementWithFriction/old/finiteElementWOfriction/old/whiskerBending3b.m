
% Whisker would point to angle base-angle if it weren't deflected by an
% object located at (x_pole, y_pole)
% 
% The deflection is computed iterativley
%
% r - vector from current length element and point of force application
% beta - angle between vertical and r
% force_angle - angle between force and vertical axis

% --------------setting whisker parameters and initialize other params

wh.prop.conical=1;
wh.prop.ym=3.5e9; % young's modulus in pascal
wh.prop.base_radius=1.00e-4; % whisker radius at the base, in meters
wh.prop.L=0.060; % whisker length in meters
wh.sim.n_segments=1000; % # of segments for finite element analysis
wh.sim.segment=wh.prop.L/wh.sim.n_segments; % in meters; total length is n_segments * segment

% Depending on whisker geometry, compute moment of inertia
if ~ wh.prop.conical
    I(1:wh.sim.n_segments)=(pi/4)*wh.prop.base_radius^4; %uniform
else
    I(1:wh.sim.n_segments)=(pi/4)*(wh.prop.base_radius*(1-(1:wh.sim.n_segments)*wh.sim.segment/wh.prop.L)).^4;  %conical
end

wh.prop.I=I;
wh.sim.base_angle = 45; % angle in degrees at the base of the whisker

wh.sim.x_pole=0.030; % x position of the object, in meters
wh.sim.y_pole=0.00; % y position of the object, in meters
wh.sim.delta_force=1/300000; % force increment per iteration in Newtons

% initialize contour in the asbences of forces == straight line
wh.prop.x(1:(wh.sim.n_segments+1))=((1:(wh.sim.n_segments+1))-1)*cos(wh.sim.base_angle *pi/180)*wh.sim.segment; % x-position of segments
wh.prop.y(1:(wh.sim.n_segments+1))=((1:(wh.sim.n_segments+1))-1)*sin(wh.sim.base_angle *pi/180)*wh.sim.segment; % y-position of segmetns
wh.prop.phi(1:wh.sim.n_segments)=wh.sim.base_angle*pi/180; % angles of whisker segments wrt to x-axis

x_old=wh.prop.x;
y_old=wh.prop.y;
phi_old=wh.prop.phi;

figure;
axes('XLim', [0 wh.prop.L], 'YLim', [-0.25*wh.prop.L 0.75*wh.prop.L]) 
line(wh.prop.x, wh.prop.y, 'Color', 'r');
line([wh.sim.x_pole wh.sim.x_pole], [-0.03 0.06], 'Color', 'r')

force=0;
temp=find(wh.prop.x>wh.sim.x_pole);
force_index=temp(1)-1;

% testing stuff
delta=[];

% --------------setting whisker parameters and initialize other params
while wh.prop.y(force_index)>wh.sim.y_pole % iterate over force
    
    temp=find(wh.prop.x>wh.sim.x_pole);
    if isempty(temp)
        'whisker slipped'
        wh.prop.x(1:(wh.sim.n_segments+1))=((1:(wh.sim.n_segments+1))-1)*cos(wh.sim.base_angle *pi/180)*wh.sim.segment; % x-position of segments
        wh.prop.y(1:(wh.sim.n_segments+1))=((1:(wh.sim.n_segments+1))-1)*sin(wh.sim.base_angle *pi/180)*wh.sim.segment; % y-position of segmetns
        wh.prop.phi(1:wh.sim.n_segments)=wh.sim.base_angle*pi/180; % angles of whisker segments wrt to x-axis
        break
    end
    if abs(wh.prop.x(min(temp))-wh.sim.x_pole) > abs(wh.prop.x(min(temp)-1)-wh.sim.x_pole)
        force_index=min(temp)-1;
    else
        force_index=min(temp);
    end
    force_index
    len_whisker=sum( sqrt(diff(wh.prop.x).*diff(wh.prop.x)+diff(wh.prop.y).*diff(wh.prop.y)))
    force=force+wh.sim.delta_force; % force increment in Newtons
    
    for k=2:(force_index-1) % finite element iteration
        r=sqrt(((wh.prop.y(force_index)-wh.prop.y(force_index-k+1))^2+(wh.prop.x(force_index)-wh.prop.x(force_index-k+1))^2));
        gamma=atan((wh.prop.y(force_index)-wh.prop.y(force_index-k+1))/(wh.prop.x(force_index)- wh.prop.x(force_index-k+1)));
        %         %    beta=acot((wh.prop.x(force_index-k+1)-wh.prop.x(force_index))/(wh.prop.y(force_index-k+1)-wh.prop.y(force_index)))+pi/2;
        %beta=abs(wh.prop.phi(force_index)-pi/2-gamma);
        moment=r*force*sin(gamma+pi/2);
        old_phi=wh.prop.phi(force_index-k+1);
        wh.prop.phi(force_index-k+1)=wh.prop.phi(force_index-k)-moment*wh.sim.segment/(wh.prop.ym*wh.prop.I(force_index-k+1));
        wh.prop.phi((force_index-k+2):end)=wh.prop.phi((force_index-k+2):end)-(old_phi-wh.prop.phi(force_index-k+1));
        wh.prop.phi(force_index:end)=wh.prop.phi(force_index-1);
        %         %   wh.prop.phi(force_index-k)=wh.prop.phi(force_index-k+1)-moment*segment/(wh.prop.ym*wh.prop.I(force_index-k+1));
        
        wh.prop.x=[0 cumsum(cos(wh.prop.phi)*wh.sim.segment)];
        wh.prop.y=[0 cumsum(sin(wh.prop.phi)*wh.sim.segment)];
        
        %          wh.prop.x(force_index-k+2)=wh.prop.x(force_index-k+1)+cos(wh.prop.phi(force_index-k+1))*wh.sim.segment;
        %          wh.prop.y(force_index-k+2)=wh.prop.y(force_index-k+1)+sin(wh.prop.phi(force_index-k+1))*wh.sim.segment;
    end
    
    if min(wh.prop.phi) < -pi/2
        'whisker slipped'
        wh.prop.x(1:(wh.sim.n_segments+1))=((1:(wh.sim.n_segments+1))-1)*cos(wh.sim.base_angle *pi/180)*wh.sim.segment; % x-position of segments
        wh.prop.y(1:(wh.sim.n_segments+1))=((1:(wh.sim.n_segments+1))-1)*sin(wh.sim.base_angle *pi/180)*wh.sim.segment; % y-position of segmetns
        wh.prop.phi(1:wh.sim.n_segments)=wh.sim.base_angle*pi/180; % angles of whisker segments wrt to x-axis
        break
    end
    
     
    temp=max(abs((wh.prop.phi-phi_old)))
    if temp > (pi/180)*5 % one degree
        force=force-wh.sim.delta_force;
        wh.sim.delta_force=wh.sim.delta_force/2;
        wh.prop.x=x_old;
        wh.prop.y=y_old;
        wh.prop.phi=phi_old;
    else
        x_old=wh.prop.x;
        y_old=wh.prop.y;
        phi_old=wh.prop.phi;
    end
    if temp < (pi/180)*1
        wh.sim.delta_force=wh.sim.delta_force*2;
    end
  
    line(wh.prop.x, wh.prop.y)
    %pause
    
%   wh.prop.x(force_index:(wh.sim.n_segments+1))=wh.prop.x(force_index-1)+(1:(wh.sim.n_segments-force_index+2))*cos(wh.prop.phi(k))*wh.sim.segment; % MOVED
%    wh.prop.y(force_index:(wh.sim.n_segments+1))=wh.prop.y(force_index-1)+(1:(wh.sim.n_segments-force_index+2))*sin(wh.prop.phi(k))*wh.sim.segment; % MOVED
    delta=[delta (wh.sim.x_pole-wh.prop.x(force_index))/wh.sim.x_pole];
    
end

% force/wh.sim.delta_force % number of force increments used

wh.prop.phi(force_index:(wh.sim.n_segments))=wh.prop.phi(force_index-1);

line(wh.prop.x,wh.prop.y, 'Color', 'k')
line(wh.prop.x(force_index),wh.prop.y(force_index), 'Color', 'r', 'Marker', 'x')
line(wh.sim.x_pole, wh.sim.y_pole, 'Color', 'g', 'Marker', 'o')


%  compute the axial force pushing into the follicle
wh.prop.axial_force=force*cos(wh.prop.phi(1)-wh.prop.phi(force_index-1)-pi/2);

%  compute the moment
gamma=atan((wh.prop.y(force_index)-wh.prop.y(1))/(wh.prop.x(force_index)- wh.prop.x(1)));
r=sqrt(((wh.prop.y(force_index)-wh.prop.y(1))^2+(wh.prop.x(force_index)-wh.prop.x(1))^2));
beta=abs(wh.prop.phi(force_index)-pi/2-gamma);
wh.prop.moment=r*force*sin(beta);

