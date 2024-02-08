clc

a1 = 1; 
a2 = 0.5; 
% a3 = 0.2;
DH1(1) = Link([0 0 a1 0]);
DH1(2) = Link([0 0 a2 0]);
% DH(3) = Link([0 0 a3 0]);
th1 = (0:0.01:pi/2) ;
th2 = (-pi/2:0.01:pi/2);
% th3 = (-pi/2:0.05:pi/2) ;
q = {th1,th2};

% full plot (remember to modify the function and delete the color parameter
plotworkspace(DH1, q)

% to use this add the color parameter

% % full mobility of link2 and link1 still at 0
% plotworkspace(DH1, {(0),(-pi/2:0.01:pi/2)}, [0.3010 0.7450 0.9330]) % light-blue
% % full mobility of link1 and link1 still at pi/2
% plotworkspace(DH1, {(0:0.01:pi/2),(pi/2)}, [0.9290 0.6940 0.1250]) % yellow
% % full mobility of link1 and link1 still at 0
% plotworkspace(DH1, {(0:0.01:pi/2),(0)}, [0.6350 0.0780 0.1840]) % wine
% % full mobility of link1 and link1 still at -pi/2
% plotworkspace(DH1, {(0:0.01:pi/2),(-pi/2)}, [0.8500 0.3250 0.0980]) % orange
% % full mobility of link2 and link1 still at pi/2
% plotworkspace(DH1, {(pi/2),(-pi/2:0.01:pi/2)}, [0 0.4470 0.7410]) % blue

% check if this point is inside WS1
x = 1.6
y = -0.2
plot(x,y,'MarkerEdge',[0.4940 0.1840 0.5560], 'LineWidth', 2, 'Marker','.', 'MarkerSize', 25)

% link 1 boundary 
r= a1 ;
teta=-pi:0.01:pi;
x=r*cos(teta);
y=r*sin(teta);
plot(x,y,"mo", 'LineWidth', 0.6)

% outer boundary
r= a1 + a2 ;
teta=-pi:0.01:pi;
x=r*cos(teta);
y=r*sin(teta);
plot(x,y,"k-", 'LineWidth', 0.6)

% inner boundary
r= sqrt(a1^2 + a2^2) ;
teta=-pi:0.01:pi;
x=r*cos(teta);
y=r*sin(teta);
plot(x,y,"k-", 'LineWidth', 0.6)

yline(0)
xline(0)

% WS2
r= abs(a2 - a1);
teta=-pi:0.01:pi;
x=r*cos(teta);
y=r*sin(teta);
plot(x,y,'MarkerFace',[0 0.4470 0.7410], 'LineWidth', 2, 'Marker','.', 'MarkerSize', 25)

hold off

function plotworkspace(DH,q)%, color)
% drawworkspace
%{
This function plots a workspace for a planar n-DOF revolute or prismatic
given DH parameters and the constraints of all variables.

This function uses Robotics Toolbox by Peter Corke which can be
downloaded from :
https://petercorke.com/wordpress/toolboxes/robotics-toolbox
----------------------------------------------
Inputs
DH    DH parameters each row is a link
q     a cell input contains constraints for all variables
        ordered from first link to last link.
---------------------------------------------------------------------

All copyrights go to Mohammad Al-Fetyani
University of Jordan
%}

% L = Link([Th d a alpha])
r = SerialLink(DH);
r.display()
[~,n] = size(DH);


var = sym('q',[n 1]);
assume(var,'real')

% generate a grid of theta1 and theta2,3,4 values
[Q{1:numel(q)}] = ndgrid(q{:}); 
T = simplify(vpa(r.fkine(var),3));
Pos = T.tv;
x(var(:)) = Pos(1);
X = matlabFunction(x);
X = X(Q{:});
y(var(:)) = Pos(2);
Y = matlabFunction(y);
Y = Y(Q{:});
plot(X(:),Y(:),'r.') 
% add color to the parameters of the function to use this 
% % plot(X(:),Y(:),'MarkerEdge',color, 'LineWidth', 1, 'Marker','o')
xlabel('X')
ylabel('Y')
hold on

end