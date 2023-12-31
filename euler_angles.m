syms phi
syms theta
syms psi

mat1 = x_m(phi);
mat2 = z_m(theta);
mat3 = y_m(psi);

display("Remember to use inverse order if the axes are fixed (i.e RPY angles)")

mat1 * mat2 * mat3

R = [0.7292 -0.6124 -0.2803;
    0.5732 0.3536 0.7392;
    0.3536 -0.7071 0.6124]

angles = rot_eul(R, "XZY");
[a,b] = rotm2eul(R, "XZY")

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




