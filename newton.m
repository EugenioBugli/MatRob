clc
format long
syms q1 q2 q3 q4 q5 q6
syms a1 a2 a3 a4
syms d0 d1 d2 d3 d4 de l1 l2 l3 l4 N L M N d A B C D K dtcp h p L1 L2
syms alpha beta gamma

%% Newton
q0 = pi/4*[-1;1;1];
r = [L*cos(q1) + N*cos(q1+q2)*cos(q3);
     L*sin(q1) + N*sin(q1+q2)*cos(q3);
     M + N*sin(q3)]

% display('as far as the det is not 0 you are ok !')
% jacobian(r, [q1,q2,q3])
% simplify(det(jacobian(r, [q1,q2,q3])))
r = subs(r, [L,M,N], [0.5,0.5,0.5])
pos = [0.3; -0.3; 0.7]

error = 1e-3
NewtonMethod(r, q0, pos, error, [q1,q2,q3])

%% Gradient

q0 = [0;0;1];

r = [q3*cos(q2)*cos(q1);
     q3*cos(q2)*sin(q1);
     d1 + q3*sin(q2)]

% display('as far as the det is not 0 you are ok !')
% jacobian(r, [q1,q2,q3])
% simplify(det(jacobian(r, [q1,q2,q3])))

r = subs(r, [d1], [0.5])
pos = [1;1;1]

error = 1e-5
minjoint = 1e-6
limitcount = 15
GradientMethod(r, q0, pos, 0.5, error, minjoint, limitcount, [q1,q2,q3])
GradientMethod(r, q0, pos, 1, error, minjoint, limitcount, [q1,q2,q3])
GradientMethod(r, q0, pos, 0.7, error, minjoint, limitcount, [q1,q2,q3])







function NewtonMethod(r, q0, pos, error, param)
    display('START Newton Method')
    e = [1;1;1];
    counter = 0;
    qi = q0;
    final = [q0];
    errornorms = [];
    while norm(e) >= error
        display(counter)
        j = jacobian(r, param);
        jac = subs(j, param, qi.');
        err = pos - subs(r, param, qi.');
        n = size(jac);

        if ~isequal(n(1), n(2))
            display('Redundant')
            invjac = pinv(jac);
        else
            invjac = inv(jac);
        end

        qf = double(qi + invjac*err)
        counter = counter + 1;
        qi = qf;
        e = double(err);
        errornorms = [errornorms norm(e)];
        final = [final qf];
    end
    display("each column will represent q from iter 0 to "+ counter)
    final
    errornorms
    y_plot = [0:counter-1];
    plot(y_plot, errornorms, LineWidth=1.5)
    title('Convergence Newton Method')
    grid on,xlabel('iterations'), ylabel('error norm')
    double(qf)
end

function GradientMethod(r, q0, pos, a, error, minjointincrem, limitcount, param)
    display('START Gradient Method')
    e = [1;1;1];
    counter = 0;
    qf = q0 + pi;
    qi = q0;
    final = [q0];
    errors = [];
    errornorms = [];
    jointerr = 1
    % remember to change something here in case of having orientation in r
    % possibly do cartesian error separate for position and orientation
    while norm(e) >= error && jointerr >= minjointincrem && counter < limitcount
        display(counter)
        j = jacobian(r, param);
        jac = subs(j, param, qi.');
        err = pos - subs(r, param, qi.');
        jacT = jac.';

        qf = double(qi + a*jacT*err)
        
        jointerr = norm(qf-qi)
        counter = counter + 1;
        qi = qf;
        e = double(err);
        errors = [errors err];
        errornorms = [errornorms norm(e)];
        final = [final qf];
    end
    display("each column will represent q from iter 0 to "+ counter)
    final
    errornorms
    final(1:3:end-3)
    y_plot = [0:counter-1];
    subplot(4,1,1)
    plot(y_plot, errornorms, LineWidth=1.5)
    title('Convergence Gradient Method')
    grid on,xlabel('iterations'), ylabel('error norm')
    hold on

    subplot(4,1,2)
    plot(y_plot, final(1:3:end-3), LineWidth=1.5)
    title('Joint q1')
    grid on,xlabel('iterations'), ylabel('q1')
    hold on
    
    subplot(4,1,3)
    plot(y_plot, final(2:3:end-3), LineWidth=1.5)
    title('Joint q2')
    grid on,xlabel('iterations'), ylabel('q2')
    hold on

    subplot(4,1,4)
    plot(y_plot, final(3:3:end-3), LineWidth=1.5)
    title('Joint q3')
    grid on,xlabel('iterations'), ylabel('q3')
    hold on

    % plot more joint variables by just adding this:
    % you need to change that for all the subplots
    % % subplot(#total num of joints,1,i+1)
    % % plot(y_plot, final(i:3:end-3), LineWidth=1.5)
    % % title('Joint qi')
    % % grid on,xlabel('iterations'), ylabel('qi')
    % % hold on
    % with i as number of that joint variable

end