syms phi
syms theta
syms psi

mat1 = z_m(phi);
mat2 = y_m(theta);
mat3 = z_m(psi);

display("Remember to use inverse order if the axes are fixed (i.e RPY angles)")

mat1 * mat2 * mat3

function mat = z_m(a)
mat = [cos(a), -sin(a), 0; sin(a), cos(a), 0; 0, 0, 1];
end
function mat = y_m(a) 
mat = [cos(a), 0, sin(a); 0, 1, 0; -sin(a), 0, cos(a)];
end

function mat = x_m(a) 
mat = [1, 0, 0; 0, cos(a), -sin(a); 0, sin(a), cos(a)];
end




