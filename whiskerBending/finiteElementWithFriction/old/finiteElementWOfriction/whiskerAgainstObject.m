
function wh = whiskerAgainstObject(wh)
%
% Rigid object is on the y-axis; together with follicle
% Calculates the contours of the whisker and the forces acting on it
%
%INPUT
% all units are MKS (meter, Newton, kilogram, etc)
% wh.pole_pos_in_meter - pole-position along y-axis; aligned with follicle; in meters
% wh.angel_at_base_degrees - angle at the base in degrees
% wh.I - 2. moment of inertia; length determines segment length, together with wh.L 
% wh.L - whisker length
% wh.ym - Young's modulus
% wh.friction - friction coefficient
%
% OUTPUT
% wh.x - x coordinate of the whisker (n+1 points)
% wh.y - y coordinate of the whisker (n+1 points)
% wh.theta - theta angle of each contour segment (n segments)
% wh.f_tot - total force acting on whisker
% wh.f_axial - axial force pushing the whisker into the face along theta base
% wh.f_norm - force acting normal to the whisker at the contact site 
% wh.f_friction - frictional force   
% wh.moment - torque acting on the whisker base
% wh.force_index - index where force acts -- i.e. point of contact along the contour
% wh.segment - segment length

r0=[0 wh.pole_pos_in_meter]'; % pole pos in meters (x, y); by definition along y-axis

tolerance=0.000001;% weird paramter - if smaller the fminsearch doesn't work; if larger, results are noisier as expected
initialCondition = [wh.pole_pos_in_meter 1/30000]; % [location of force along the whisker (m), force acting on the whisker (N)], initial condition

theta0=pi*wh.angel_at_base_degrees/180;
rotTheta0=[cos(-theta0) -sin(-theta0); sin(-theta0) cos(-theta0)];
r1=rotTheta0*r0;
rotTheta0=[cos(-theta0) -sin(-theta0); sin(-theta0) cos(-theta0)];
a=r1';

options = optimset('TolFun',tolerance);
[f_fit,fval,exitflag,output] = fminsearch(@(f_fit) whiskerBentByForce(f_fit, wh, a),initialCondition, options);

if fval < tolerance*100 % Kluge; I don't understand why fminsearch 'converges' even though the 'tolerance' is never reached
    [temp, wh]=whiskerBentByForce(f_fit, wh);
    for i=1:size(wh.x, 2)
        temp=[wh.x(i) wh.y(i)]';
        temp=rotTheta0*temp;
        wh.x(i)=temp(1);
        wh.y(i)=temp(2);
    end
    wh.theta=wh.theta+theta0;
    %plot(wh.x*1000, -wh.y*1000, 'r', 'LineWidth', 1.0); % minus is to keep
    %with convention
else
    error('no convergence')
end

