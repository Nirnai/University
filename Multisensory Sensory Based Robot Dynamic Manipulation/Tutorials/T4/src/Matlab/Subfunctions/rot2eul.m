function [ phi_d,theta_d,psi_d ] = rot2eul( R )

d = sqrt(1-R(3,1)^2);

if(abs(R(3,1)) ~= 1)
    theta_d= atan2(-R(3,1), d);
    phi_d = atan2(R(2,1),R(1,1));
    psi_d = atan2(R(3,2),R(3,3));
else
    % Gimbal Lock case occurs when theta_d = +- pi/2
    % Gimbal Lock case has infinite solutions. It is convinient
    % to set phi_d to 0 to get one solution
    phi_d = 0;
    if(R(3,1) == -1)
        theta_d = pi/2;
        psi_d = phi_d + atan2(R(1,2),R(1,3));
    else
        theta_d = -pi/2;
        psi_d = -phi_d + atan2(-R(1,2),-R(1,3));
    end
end

end

