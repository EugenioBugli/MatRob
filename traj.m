syms t omega

p = 3+ [sin(2*t); 1.5*sin(t)]
subplot(5,1,1)
fplot(p(1), p(2), LineWidth=1.5)

title('Position')
xlabel('s')
ylabel('m')

pdot = diff(p)
subplot(5,1,2)
fplot(pdot(1), pdot(2), LineWidth=1.5)

title('Velocity')
xlabel('s')
ylabel('m/s')

pddot = diff(pdot)
subplot(5,1,3)
fplot(pddot(1), pddot(2), LineWidth=1.5)

title('Accelleration')
xlabel('s')
ylabel('m/s^2')


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%% Joint Velocities
pomega = [3+sin(2*0.6124*t) 3+1.5*sin(0.6124*t)]


pomegadot = diff(pomega)
subplot(5,1,4)
fplot(pomegadot, LineWidth=1.5)
title('Velocity')
xlabel('s')
ylabel('m/s')

pomegaddot = diff(pomegadot)

subplot(5,1,5)
fplot(pomegaddot, LineWidth=1.5)
title('Accelleration')
xlabel('s')
ylabel('m/s^2')

%% Bang Coast Bang  = Rest to Rest
syms t

amax = 4;
vmax = 2;

so = -pi/2;
sf = pi/2;

Ts = vmax/amax;
T = ((sf-so)*amax + vmax^2)/(amax*vmax);
intervals = [t >= 0 & t<=Ts; t>=Ts & t<=(T-Ts); t>=(T-Ts) & t<=T];

sigma = so + [0.5*amax*t^2,vmax*t - vmax^2/(2*amax), -0.5*amax*(t-T)^2 + vmax*T - vmax^2/amax]


% position
pos = piecewise(intervals(1),sigma(1),intervals(2),sigma(2),intervals(3),sigma(3));
subplot(3,1,1)
fplot(pos, LineWidth=1.5)
ylim([-max(abs(so),abs(sf)) max(abs(so),abs(sf))]),grid on, xlabel('t'), ylabel('x(t)'), title('Position')


% velocity
sigmadot = diff(sigma)
pos = piecewise(intervals(1),sigmadot(1),intervals(2),sigmadot(2),intervals(3),sigmadot(3));
subplot(3,1,2)
fplot(pos, LineWidth=1.5)
ylim([-max(abs(vo),abs(vf)) max(abs(vo),abs(vf))]), grid on, xlabel('t'), ylabel('dx(t)'), title('Velocity')


% accelleration
sigmaddot = diff(sigmadot)
pos = piecewise(intervals(1),sigmaddot(1),intervals(2),sigmaddot(2),intervals(3),sigmaddot(3));
subplot(3,1,3)
fplot(pos, LineWidth=1.5)
grid on, xlabel('t'), ylabel('ddx(t)'), title('Accelleration')