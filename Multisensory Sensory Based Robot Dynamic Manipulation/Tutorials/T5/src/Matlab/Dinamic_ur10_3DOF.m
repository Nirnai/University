function qdd = Dinamic_ur10_3DOF(u)
    
    q1 = u(1);
    q2 = u(2);
    q3 = u(3);
    q = [q1; q2; q3];
    q1p = u(4);
    q2p = u(5);
    q3p = u(6);
    qp=[q1p; q2p; q3p];

    %% External forces
    % Damping
    b1 = 1;
    b2 = 1;
    b3 = 1;
    
    Beta = zeros(3);
    Beta(1,1) = b1;
    Beta(2,2) = b2;
    Beta(3,3) = b3;

    % Euler Langrang Matrices
    [M,C,G] = EulerLagrangeMatrices(q,qp);

    % Dynamics Equation
    qdd=(M)\(-C*qp-G-Beta*qp);
end

