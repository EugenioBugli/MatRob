syms phi
syms theta
syms psi
syms a
%%
mat1 = y_m(deg2rad(45));
mat2 = x_m(deg2rad(-45));
mat3 = y_m(deg2rad(120));

R0i = mat1*mat2*mat3

R0f = [0 sin(pi/3) cos(pi/3); 0 cos(pi/3) -sin(pi/3); -1 0 0]

Rif = R0i.'*R0f

inv_pos(Rif)

w = 1.1*[0.5170; -0.3752; 0.7694]
%%
% vector invariant to rotation --> eigenvector of eigenvalue 1 : [V,D] =
% eig(rot) where D return diagonal matrix with eigenvalues, while V
% represent the eigenvectors of the 3 eigenvalues.

% when you solve the direct problem remember to express the vector respect
% to 0RF so you need to compute this : Ri * direction
R2 = 0.5*[sqrt(2) 1 -1;
      0 -sqrt(2) -sqrt(2);
      -sqrt(2) 1 -1]

Rvia = [sqrt(6)/4 sqrt(2)/4 -sqrt(2)/2;
        -sqrt(6)/4 -sqrt(2)/4 -sqrt(2)/2;
        -0.5 sqrt(3)/2 0]
% Ri = [sqrt(2)/2,0,sqrt(2)/2;
%       0,-1,0;
%       sqrt(2)/2,0,-sqrt(2)/2]
% Rf = x_m(-pi/2)*y_m(pi/3)*z_m(pi/3)
% 
% rot = Ri.'*Rf
% 
% dir = inv_pos(rot)
% 
% Ri * dir
% Remember that when you return the direction vector you need to express it
% in frame 0 coordinates : Ri * vector obtained from inv_pos

% r = (1/sqrt(2)) * [1;-1;0]
% Rx = x_m(deg2rad(-90))
% r_n = Rx * r
% Rot = dir_pos(deg2rad(60),r_n)
% dir = inv_pos(Rot)
function rot = dir_pos(angle,dir)
    display(dir)
    display(angle)
    display(rad2deg(angle))
    rot = dir * dir.' + (eye(3)-(dir*dir.'))*cos(angle) + skew(dir)*sin(angle);
end

function direction = inv_pos(r_mat)
    if int64(det(r_mat)) == 1
        y = sqrt((r_mat(1,2) - r_mat(2,1))^2 + (r_mat(1,3) - r_mat(3,1))^2 + (r_mat(2,3) - r_mat(3,2))^2);
        x = (trace(r_mat)-1);
        t = atan2(y,x);
        if y == 0.0
            display("Singularity")
            dir_1 = [+sqrt((r_mat(1,1)+1)/2);+sqrt((r_mat(2,2)+1)/2);+sqrt((r_mat(3,3)+1)/2)];
            dir_2 = - dir_1;
            direction = [dir_1,dir_2];
            display("Remember check 2rxry = R12, 2rxrz = R13, 2ryrz = R23")
            angles = atan2(y,x);
            if angles == 0.0
                display("No direction solution")
                direction = NaN
            end
        else
            angles = [atan2(y,x),atan2(-y,x)];
            dir = [r_mat(3,2)-r_mat(2,3);
                   r_mat(1,3)-r_mat(3,1);
                   r_mat(2,1)-r_mat(1,2)];
            dir_1 = dir / (2*sin(atan2(y,x)));
            k = -y;
            dir_2 = dir / (2*sin(atan2(k,x)));
            direction = [dir_1, dir_2];
        end
        display("Radiant")
        display(angles)
        display("Degree")
        display(rad2deg(angles))
    else
        display("Det is not equal to 1")
    end
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