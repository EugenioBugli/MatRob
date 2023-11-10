% Ri = [sqrt(2)/2,0,sqrt(2)/2;
%       0,-1,0;
%       sqrt(2)/2,0,-sqrt(2)/2]
% Rf = x_m(-pi/2)*y_m(pi/3)*z_m(pi/3)
% 
% rot = Ri.'*Rf
% 
% check(rot)
format long
syms q1 q2 q3 d1 a3
syms a1 o1 d1 t1
dh(sym(-pi/2),0,0,0)
y_m(sym(-pi/2))*x_m(sym(-pi/2))
meeee = dh(sym(-pi/2),0,0,0) * direct(sym([0,q2,0]),sym([pi/2,-pi/2,0]),[q1,0,q3],[0,0,0],3)  * [0 1 0 0;0 0 1 0; 1 0 0 0; 0 0 0 1]
%m = direct(sym([0,0,q2,0,0]),sym([-pi/2,pi/2,pi/2,0,0]),[0,q1,0,q3,0],[0,0,0,0,0],5);

function mat = direct(theta,alpha,d,a,num)

    mat = dh(alpha(1),theta(1),d(1),a(1))

    for i = 2:num
        dh(alpha(i),theta(i),d(i),a(i))
        mat = mat * dh(alpha(i),theta(i),d(i),a(i));
        simplify(mat);
    end
end

function mdh = dh(alpha,t,d,a)
    mdh = [cos(t),  -cos(alpha)*sin(t),     sin(alpha)*sin(t),     a*cos(t);
           sin(t),   cos(alpha)*cos(t),    -sin(alpha)*cos(t),     a*sin(t);
           0,               sin(alpha),            cos(alpha),            d;
           0                     0,                 0,            1];
end

function bool = check(m)
    deter = int64(det(m)) == 1;
    nor1 = int64(norm(m(:,1))) == 1;
    nor2 = int64(norm(m(:,2))) == 1;
    nor3 = int64(norm(m(:,3))) == 1;
    orthog = isequal(inv(m),m.');
    bool = deter & nor1 & nor2 & nor3 & orthog;
end

function mat = z_m(a)
mat = [cos(a), -sin(a), 0; sin(a), cos(a), 0; 0, 0, 1];
end
function mat = y_m(a) 
mat = [cos(a), 0, sin(a); 0, 1, 0; -sin(a), 0, cos(a)];
end

function mat = x_m(a) 
mat = [1, 0, 0; 0, cos(a), -sin(a); 0, sin(a), cos(a)];
end

function s = skew(r)
s = [0,-r(3),r(2);
     r(3),0,-r(1);
     -r(2),r(1),0];
end