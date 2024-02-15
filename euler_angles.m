clc 
syms phi
syms theta
syms psi
syms alpha beta gamma

mat1 = y_m(deg2rad(45));
mat2 = x_m(deg2rad(-45));
mat3 = y_m(deg2rad(120));

R0i = mat1*mat2*mat3

R0f = [0 sin(pi/3) cos(pi/3); 0 cos(pi/3) -sin(pi/3); -1 0 0]

Rif = R0i.'*R0f

display("Remember to use inverse order if the axes are fixed (i.e RPY angles)")
display("You may use wrapToPi to transform your angles in [-pi,pi] (note that we use (-pi,pi])")

Rin = [0 1 0; 0.5 0 sqrt(3)/2; sqrt(3)/2 0 -0.5]
% Rfin = [sqrt(2)/2, -sqrt(2)/2, 0; -0.5, -0.5, -sqrt(2)/2; 0.5, 0.5, -sqrt(2)/2]
% 
% R = Rin.' * Rfin

rot_eul(Rin, 'YXZ')
[a, b] = rotm2eul(Rin, 'YXZ')
rad2deg(a)
rad2deg(b)

% R = [0 -sqrt(2)/2 sqrt(2)/2;
%     1 0 0;
%     0 sqrt(2)/2 sqrt(2)/2]
% R2 = 0.5*[sqrt(2) 1 -1;
%       0 -sqrt(2) -sqrt(2);
%       -sqrt(2) 1 -1]
% 
% Rvia = [sqrt(6)/4 sqrt(2)/4 -sqrt(2)/2;
%         -sqrt(6)/4 -sqrt(2)/4 -sqrt(2)/2;
%         -0.5 sqrt(3)/2 0]
% 
% angles = rot_eul(Rvia, "ZYX");
% [a,b] = rotm2eul(Rvia, "ZYX")
% display("mtlb function with degrees")
% rad2deg(a)
% rad2deg(b)
% 
% R = x_m(sym(pi/4))*y_m(sym(-pi/4))*z_m(sym(-pi/2))
% [V,D] = eig(R)

% Ri = [0 0.5 -sqrt(3)/2; -1 0 0; 0 sqrt(3)/2 0.5]
% R = Ri.' * eye(3)
% an = rot_eul(R, "ZYZ");
% [a,b] = rotm2eul(R, "ZYZ")
% wrapToPi(a)

% Rp = z_m(alpha)*y_m(beta)
% rp = z_m(deg2rad(45))*y_m(-deg2rad(60))
% rot = x_m(phi)*y_m(theta)*z_m(psi)
% 
% an = rot_eul(rp, "XYZ");
% [a,b] = rotm2eul(rp, "XYZ")
% 
% Rfin1 = x_m(0.8861)*y_m(-0.6591)*z_m(1.1071)
% Rfin2 = x_m(-2.2555)*y_m(-2.4825)*z_m(-2.0344)

function angles = rot_eul(R,seq)
    %----------------------------------
    
    % X FIRST

    %----------------------------------
    if seq == "XYX"
        xtheta = R(1,1);
        ytheta = sqrt(R(2,1)^2 + R(3,1)^2);
        if ytheta ~= 0.0
            ytheta1 = ytheta;
            ytheta2 = -ytheta;
            
            theta1 = atan2(ytheta1,xtheta);
            theta2 = atan2(ytheta2, xtheta);
            
            xphi1 = -R(3,1) / ytheta1;
            xphi2 = -R(3,1) / ytheta2;
            yphi1 = R(2,1) / ytheta1;
            yphi2 = R(2,1) / ytheta2;

            phi1 = atan2(yphi1, xphi1);
            phi2 = atan2(yphi2, xphi2);

            xpsi1 = R(1,3) / ytheta1;
            xpsi2 = R(1,3) / ytheta2;
            ypsi1 = R(1,2) / ytheta1;
            ypsi2 = R(1,2) / ytheta2;

            psi1 = atan2(ypsi1, xpsi1);
            psi2 = atan2(ypsi2, xpsi2);
            
            A = [phi1; theta1; psi1]
            B = [phi2; theta2; psi2]

            angles = 0;
        end
        theta = atan2(ytheta, xtheta)
        if theta == 0.0
            
        end
    %----------------------------------
    elseif seq == "XZX"
        xtheta = R(1,1);
        ytheta = sqrt(R(2,1)^2 + R(3,1)^2);
        if ytheta ~= 0.0
            ytheta1 = ytheta;
            ytheta2 = -ytheta;
            
            theta1 = atan2(ytheta1,xtheta);
            theta2 = atan2(ytheta2, xtheta);
            
            xphi1 = R(2,1) / ytheta1;
            xphi2 = R(2,1) / ytheta2;
            yphi1 = R(3,1) / ytheta1;
            yphi2 = R(3,1) / ytheta2;

            phi1 = atan2(yphi1, xphi1);
            phi2 = atan2(yphi2, xphi2);

            xpsi1 = -R(1,2) / ytheta1;
            xpsi2 = -R(1,2) / ytheta2;
            ypsi1 = R(1,3) / ytheta1;
            ypsi2 = R(1,3) / ytheta2;

            psi1 = atan2(ypsi1, xpsi1);
            psi2 = atan2(ypsi2, xpsi2);
            
            A = [phi1; theta1; psi1]
            B = [phi2; theta2; psi2]

            angles = 0;
        end
        theta = atan2(ytheta, xtheta)
        if theta == 0.0
            
        end
    %----------------------------------
    elseif seq == "XYZ"

        ytheta = R(1,3);
        xtheta = sqrt(R(2,3)^2 + R(3,3)^2);
        if xtheta ~= 0.0
            xtheta1 = xtheta;
            xtheta2 = -xtheta;
            
            theta1 = atan2(ytheta,xtheta1);
            theta2 = atan2(ytheta, xtheta2);
            
            xphi1 = R(3,3) / xtheta1;
            xphi2 = R(3,3) / xtheta2;
            yphi1 = -R(2,3) / xtheta1;
            yphi2 = -R(2,3) / xtheta2;

            phi1 = atan2(yphi1, xphi1);
            phi2 = atan2(yphi2, xphi2);

            xpsi1 = R(1,1) / xtheta1;
            xpsi2 = R(1,1) / xtheta2;
            ypsi1 = -R(1,2) / xtheta1;
            ypsi2 = -R(1,2) / xtheta2;

            psi1 = atan2(ypsi1, xpsi1);
            psi2 = atan2(ypsi2, xpsi2);
            
            A = [phi1; theta1; psi1]
            B = [phi2; theta2; psi2]

            angles = 0;
        end
        theta = atan2(ytheta, xtheta)
        if theta == 0.0
            
        end
    %----------------------------------
       elseif seq == "XZY"

        ytheta = -R(1,2);
        xtheta = sqrt(R(2,2)^2 + R(3,2)^2);
        if xtheta ~= 0.0
            xtheta1 = xtheta;
            xtheta2 = -xtheta;
            
            theta1 = atan2(ytheta,xtheta1);
            theta2 = atan2(ytheta, xtheta2);
            
            xphi1 = R(2,2) / xtheta1;
            xphi2 = R(2,2) / xtheta2;
            yphi1 = R(3,2) / xtheta1;
            yphi2 = R(3,2) / xtheta2;

            phi1 = atan2(yphi1, xphi1);
            phi2 = atan2(yphi2, xphi2);

            xpsi1 = R(1,1) / xtheta1;
            xpsi2 = R(1,1) / xtheta2;
            ypsi1 = R(1,3) / xtheta1;
            ypsi2 = R(1,3) / xtheta2;

            psi1 = atan2(ypsi1, xpsi1);
            psi2 = atan2(ypsi2, xpsi2);
            
            A = [phi1; theta1; psi1]
            B = [phi2; theta2; psi2]

            angles = 0;
        end
        theta = atan2(ytheta, xtheta)
        if theta == 0.0
            
        end
    %----------------------------------
    
    % Y FIRST

    %----------------------------------
    elseif seq == "YXY"
        xtheta = R(2,2);
        ytheta = sqrt(R(2,1)^2 + R(2,3)^2);
        if ytheta ~= 0.0
            ytheta1 = ytheta;
            ytheta2 = -ytheta;
            
            theta1 = atan2(ytheta1,xtheta);
            theta2 = atan2(ytheta2, xtheta);
            
            xphi1 = R(3,2) / ytheta1;
            xphi2 = R(3,2) / ytheta2;
            yphi1 = R(1,2) / ytheta1;
            yphi2 = R(1,2) / ytheta2;

            phi1 = atan2(yphi1, xphi1);
            phi2 = atan2(yphi2, xphi2);

            xpsi1 = -R(2,3) / ytheta1;
            xpsi2 = -R(2,3) / ytheta2;
            ypsi1 = R(2,1) / ytheta1;
            ypsi2 = R(2,1) / ytheta2;

            psi1 = atan2(ypsi1, xpsi1);
            psi2 = atan2(ypsi2, xpsi2);
            
            A = [phi1; theta1; psi1]
            B = [phi2; theta2; psi2]

            angles = 0;
        end
        theta = atan2(ytheta, xtheta)
        if theta == 0.0
            
        end
    %----------------------------------    
    elseif seq == "YZY"
        xtheta = R(2,2);
        ytheta = sqrt(R(2,1)^2 + R(2,3)^2);
        if ytheta ~= 0.0
            ytheta1 = ytheta;
            ytheta2 = -ytheta;
            
            theta1 = atan2(ytheta1,xtheta);
            theta2 = atan2(ytheta2, xtheta);
            
            xphi1 = -R(1,2) / ytheta1;
            xphi2 = -R(1,2) / ytheta2;
            yphi1 = R(3,2) / ytheta1;
            yphi2 = R(3,2) / ytheta2;

            phi1 = atan2(yphi1, xphi1);
            phi2 = atan2(yphi2, xphi2);

            xpsi1 = R(2,1) / ytheta1;
            xpsi2 = R(2,1) / ytheta2;
            ypsi1 = R(2,3) / ytheta1;
            ypsi2 = R(2,3) / ytheta2;

            psi1 = atan2(ypsi1, xpsi1);
            psi2 = atan2(ypsi2, xpsi2);
            
            A = [phi1; theta1; psi1]
            B = [phi2; theta2; psi2]

            angles = 0;
        end
        theta = atan2(ytheta, xtheta)
        if theta == 0.0
            
        end
    %----------------------------------
    elseif seq == "YXZ"

        ytheta = -R(2,3);
        xtheta = sqrt(R(1,3)^2 + R(3,3)^2);
        if xtheta ~= 0.0
            xtheta1 = xtheta;
            xtheta2 = -xtheta;
            
            theta1 = atan2(ytheta,xtheta1);
            theta2 = atan2(ytheta, xtheta2);
            
            xphi1 = R(3,3) / xtheta1;
            xphi2 = R(3,3) / xtheta2;
            yphi1 = R(1,3) / xtheta1;
            yphi2 = R(1,3) / xtheta2;

            phi1 = atan2(yphi1, xphi1);
            phi2 = atan2(yphi2, xphi2);

            xpsi1 = R(2,2) / xtheta1;
            xpsi2 = R(2,2) / xtheta2;
            ypsi1 = R(2,1) / xtheta1;
            ypsi2 = R(2,1) / xtheta2;

            psi1 = atan2(ypsi1, xpsi1);
            psi2 = atan2(ypsi2, xpsi2);
            
            A = [phi1; theta1; psi1]
            B = [phi2; theta2; psi2]

            angles = 0;
        end
        theta = atan2(ytheta, xtheta)
        if theta == 0.0
            
        end
    %---------------------------------- 
    elseif seq == "YZX"

        ytheta = R(2,1);
        xtheta = sqrt(R(1,1)^2 + R(3,1)^2);
        if xtheta ~= 0.0
            xtheta1 = xtheta;
            xtheta2 = -xtheta;
            
            theta1 = atan2(ytheta,xtheta1);
            theta2 = atan2(ytheta, xtheta2);
            
            xphi1 = R(1,1) / xtheta1;
            xphi2 = R(1,1) / xtheta2;
            yphi1 = -R(3,1) / xtheta1;
            yphi2 = -R(3,1) / xtheta2;

            phi1 = atan2(yphi1, xphi1);
            phi2 = atan2(yphi2, xphi2);

            xpsi1 = R(2,2) / xtheta1;
            xpsi2 = R(2,2) / xtheta2;
            ypsi1 = -R(2,3) / xtheta1;
            ypsi2 = -R(2,3) / xtheta2;

            psi1 = atan2(ypsi1, xpsi1);
            psi2 = atan2(ypsi2, xpsi2);
            
            A = [phi1; theta1; psi1]
            B = [phi2; theta2; psi2]

            angles = 0;
        end
        theta = atan2(ytheta, xtheta)
        if theta == 0.0
            
        end
    %----------------------------------
    
    % Z FIRST

    %----------------------------------
    elseif seq == "ZXZ"
        xtheta = R(3,3);
        ytheta = sqrt(R(1,3)^2 + R(2,3)^2);
        if ytheta ~= 0.0
            ytheta1 = ytheta;
            ytheta2 = -ytheta;
            
            theta1 = atan2(ytheta1,xtheta);
            theta2 = atan2(ytheta2, xtheta);
            
            xphi1 = -R(2,3) / ytheta1;
            xphi2 = -R(2,3) / ytheta2;
            yphi1 = R(1,3) / ytheta1;
            yphi2 = R(1,3) / ytheta2;

            phi1 = atan2(yphi1, xphi1);
            phi2 = atan2(yphi2, xphi2);

            xpsi1 = R(3,2) / ytheta1;
            xpsi2 = R(3,2) / ytheta2;
            ypsi1 = R(3,1) / ytheta1;
            ypsi2 = R(3,1) / ytheta2;

            psi1 = atan2(ypsi1, xpsi1);
            psi2 = atan2(ypsi2, xpsi2);
            
            A = [phi1; theta1; psi1]
            B = [phi2; theta2; psi2]

            angles = 0;
        end
        theta = atan2(ytheta, xtheta)
        if theta == 0.0
            
        end
    %----------------------------------
    elseif seq == "ZYZ"
        xtheta = R(3,3);
        ytheta = sqrt(R(1,3)^2 + R(2,3)^2);
        if ytheta ~= 0.0
            ytheta1 = ytheta;
            ytheta2 = -ytheta;
            
            theta1 = atan2(ytheta1,xtheta);
            theta2 = atan2(ytheta2, xtheta);
            
            xphi1 = R(1,3) / ytheta1;
            xphi2 = R(1,3) / ytheta2;
            yphi1 = R(2,3) / ytheta1;
            yphi2 = R(2,3) / ytheta2;

            phi1 = atan2(yphi1, xphi1);
            phi2 = atan2(yphi2, xphi2);

            xpsi1 = -R(3,1) / ytheta1;
            xpsi2 = -R(3,1) / ytheta2;
            ypsi1 = R(3,2) / ytheta1;
            ypsi2 = R(3,2) / ytheta2;

            psi1 = atan2(ypsi1, xpsi1);
            psi2 = atan2(ypsi2, xpsi2);
            
            A = [phi1; theta1; psi1]
            B = [phi2; theta2; psi2]

            angles = 0;
        end
        theta = atan2(ytheta, xtheta)
        if theta == 0.0
            
        end
        %----------------------------------
        elseif seq == "ZXY"

        ytheta = R(3,2);
        xtheta = sqrt(R(1,2)^2 + R(2,2)^2);
        if xtheta ~= 0.0
            xtheta1 = xtheta;
            xtheta2 = -xtheta;
            
            theta1 = atan2(ytheta,xtheta1);
            theta2 = atan2(ytheta, xtheta2);
            
            xphi1 = R(2,2) / xtheta1;
            xphi2 = R(2,2) / xtheta2;
            yphi1 = -R(1,2) / xtheta1;
            yphi2 = -R(1,2) / xtheta2;

            phi1 = atan2(yphi1, xphi1);
            phi2 = atan2(yphi2, xphi2);

            xpsi1 = R(3,3) / xtheta1;
            xpsi2 = R(3,3) / xtheta2;
            ypsi1 = -R(3,1) / xtheta1;
            ypsi2 = -R(3,1) / xtheta2;

            psi1 = atan2(ypsi1, xpsi1);
            psi2 = atan2(ypsi2, xpsi2);
            
            A = [phi1; theta1; psi1]
            B = [phi2; theta2; psi2]

            angles = 0;
        end
        theta = atan2(ytheta, xtheta)
        if theta == 0.0
            
        end
        %----------------------------------
        elseif seq == "ZYX"

        ytheta = -R(3,1);
        xtheta = sqrt(R(1,1)^2 + R(2,1)^2);
        if xtheta ~= 0.0
            xtheta1 = xtheta;
            xtheta2 = -xtheta;
            
            theta1 = atan2(ytheta,xtheta1);
            theta2 = atan2(ytheta, xtheta2);
            
            xphi1 = R(1,1) / xtheta1;
            xphi2 = R(1,1) / xtheta2;
            yphi1 = R(2,1) / xtheta1;
            yphi2 = R(2,1) / xtheta2;

            phi1 = atan2(yphi1, xphi1);
            phi2 = atan2(yphi2, xphi2);

            xpsi1 = R(3,3) / xtheta1;
            xpsi2 = R(3,3) / xtheta2;
            ypsi1 = R(3,2) / xtheta1;
            ypsi2 = R(3,2) / xtheta2;

            psi1 = atan2(ypsi1, xpsi1);
            psi2 = atan2(ypsi2, xpsi2);
            
            A = [phi1; theta1; psi1]
            B = [phi2; theta2; psi2]

            angles = 0;
        end
        theta = atan2(ytheta, xtheta)
        if theta == 0.0
            
        end
    %----------------------------------
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




