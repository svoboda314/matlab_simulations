
pole_pos=5;
force=80;

xy=calcWhiskerShape(pole_pos, force);

h23=figure(23);
plot(xy(1, :), xy(2, :), 'LineWidth', 2) 
axis equal
a_lim=max(xy(1, :));
set(gca, 'YLim', [0, a_lim], 'XLim', [0, a_lim], 'YDir', 'reverse')

hold
xy1=simWhiskerShape(pole_pos, force);
plot(xy1(1, :), xy1(2, :), 'r') 
hold off

h24=figure(24);
fe=2*abs(xy(2, :)-xy1(2, :))./(xy(2, :)+xy1(2, :));
plot(fe)
set(gca, 'YLim', [0, 1])

save('xy_5_2', 'xy1')
