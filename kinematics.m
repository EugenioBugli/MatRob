Ri = [sqrt(2)/2,0,sqrt(2)/2;
      0,-1,0;
      sqrt(2)/2,0,-sqrt(2)/2]
Rf = x_m(-pi/2)*y_m(pi/3)*z_m(pi/3)

rot = Ri.'*Rf

check(rot)



function mat = direct(theta,alpha,d,num,base)
    mat = base
    for i = 1:num
        mat = mat * dh(alpha[i],theta[i],d[i]);
        simplify(mat);
    end
end

function mdh = dh(a,t,d)
    mdh = [cos(t),  -cos(a)*sin(t),     sin(a)*sin(t),     a*cos(t);
           sin(t),   cos(a)*cos(t),    -sin(a)*cos(t),     a*sin(t);
           0,               sin(a),            cos(a),            d;
           0                     0,                 0,            1]
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