
n_segments=1000;
wh.pole_pos_in_meter=0.005;
wh.angel_at_base_degrees=19;

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

wh = whiskerAgainstObject(wh);

%% plot whisker with force vectors

scale_factor=1/1e-5;

hf=figure;
ha=axes;
axis('equal')
set(ha, 'YDir', 'reverse', 'FontSize', 16, 'LineWidth', 1)%, 'DataAspectRatio', [1 1 1]);
xlabel('Distance (mm)')
ylabel('Distance (mm)')
line(wh.x*1000, wh.y*1000, 'Color', 'k', 'LineWidth', 1)
line(0, 0, 'Color', 'k', 'Marker', 'x', 'MarkerSize', 10)
line(0, wh.pole_pos_in_meter*1000, 'Color', 'k', 'Marker', 'o', 'MarkerSize', 10)
%line(u, v, 'Color', 'r', 'LineWidth', 2)

force_point=[wh.x(wh.force_index)*1000 wh.y(wh.force_index)*1000];
u_tangential=[sin(wh.theta(wh.force_index)) cos(wh.theta(wh.force_index))]; 
u_normal=[cos(wh.theta(wh.force_index)) -sin(wh.theta(wh.force_index))]; 
u_axial=[-sin(wh.theta(1)) -cos(wh.theta(1))];

x_normal=[force_point(1) force_point(1)-u_normal(1)*scale_factor*wh.f_norm];
y_normal=[force_point(2) force_point(2)-u_normal(2)*scale_factor*wh.f_norm];
line(x_normal, y_normal, 'Color', 'r', 'LineStyle', ':', 'LineWidth', 2)
line(x_normal(2), y_normal(2), 'Color', 'r', 'Marker', 'o', 'MarkerSize', 5)

x_tangential=[force_point(1) force_point(1)+u_tangential(1)*scale_factor*wh.f_friction];
y_tangential=[force_point(2) force_point(2)+u_tangential(2)*scale_factor*wh.f_friction];
line(x_tangential, y_tangential, 'Color', 'r', 'LineStyle', ':', 'LineWidth', 2)
line(x_tangential(2), y_tangential(2), 'Color', 'r', 'Marker', 'o', 'MarkerSize', 5)

x_axial=[0 u_axial(1)*scale_factor*wh.f_axial];
y_axial=[0 u_axial(2)*scale_factor*wh.f_axial];
line(x_axial, y_axial, 'Color', 'r', 'LineStyle', ':', 'LineWidth', 2)
line(x_axial(2), y_axial(2), 'Color', 'r', 'Marker', 'o', 'MarkerSize', 5)

%% overlay whiskers without (black) and with (magenta) friction

f_coeff= 0.4;
n_segments=1000;
wh.pole_pos_in_meter=0.008;
wh.angel_at_base_degrees=19;

conical=1; % if 1 then conical whisker with linear taper
wh.ym=5e9; % young's modulus in pascal
base_radius= 33.5e-6; % in meters
wh.L = 16e-3; % whisker lenght in meters


segment=wh.L/n_segments; % in meters; total length is n_segments * segment
if ~ conical
    I(1:n_segments)=(pi/4)*base_radius^4; %uniform
else
    I(1:n_segments)=(pi/4)*(base_radius*(1-(1:n_segments)*segment/wh.L)).^4;  %conical
end
wh.I=I; % 2. moment of intertia

wh.friction = 0;
wh = whiskerAgainstObject(wh);

hf=figure;
ha=axes;
set(ha, 'YDir', 'reverse', 'FontSize', 16, 'LineWidth', 1)%, 'DataAspectRatio', [1 1 1]);
xlabel('Distance (mm)')
ylabel('Distance (mm)')
line(wh.x*1000, wh.y*1000, 'Color', 'k', 'LineWidth', 1)
line(0, wh.pole_pos_in_meter*1000, 'Color', 'k', 'Marker', 'o', 'MarkerSize', 10)

wh.friction = f_coeff;
wh = whiskerAgainstObject(wh);
line(wh.x*1000, wh.y*1000, 'Color', 'm', 'LineWidth', 1)
