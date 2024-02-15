syms t A T
int(A*T - A*t,t, [T/2 t - T/2])
%% 
amax = 200;

so = 36.87;
sf = 135;

% Ts = vmax/amax;
% T = ((sf-so)*amax + vmax^2)/(amax*vmax);
T = 1.401
Ts = T/2

intervals = [t >= 0 & t<Ts; t>=Ts & t<=T] % [0, T/2) and [T/2, T]

sigma = so + [0.5*amax*t^2, -0.5*amax*(t-Ts)^2 + amax*Ts*(t-Ts) + 0.125*amax*T^2]

% position
pos = piecewise(intervals(1),sigma(1),intervals(2),sigma(2));
subplot(3,1,1)
fplot(pos, LineWidth=1.5)
% ylim([-max(abs(so),abs(sf)) max(abs(so),abs(sf))])
xlim([0, T])
grid on, xlabel('t'), ylabel('x(t)'), title('Position')


% velocity
vo = 20
sigmadot = diff(sigma)
% sigmadot = vo + sigmadot
pos = piecewise(intervals(1),sigmadot(1),intervals(2),sigmadot(2));
subplot(3,1,2)
fplot(pos, LineWidth=1.5)
xlim([0, T])
grid on, xlabel('t'), ylabel('dx(t)'), title('Velocity')


% accelleration
sigmaddot = diff(sigmadot)
pos = piecewise(intervals(1),sigmaddot(1),intervals(2),sigmaddot(2));
subplot(3,1,3)
fplot(pos, LineWidth=1.5)
xlim([0, T])
grid on, xlabel('t'), ylabel('ddx(t)'), title('Accelleration')