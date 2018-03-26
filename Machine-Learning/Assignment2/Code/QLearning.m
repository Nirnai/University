function [ initstate, Q, delta, rew ] = QLearning(delta, rew, Q, eps, initstate, gamma, alpha, steps )
states = 16;
ctr = 1;
converged= 0;
while(ctr < steps)
    rdm = rand;
    if rdm > eps
        [~, a] = max(Q(states,:));
    else
        a = randi(4);
    end

    [state_next,r] = SimulateRobot(states,a,delta,rew);
    
%   save old Q to check for convergence    
    Q_old = Q;
    
    Q(states,a) = Q(states,a) + alpha*(r+gamma*max(Q(state_next,:))-Q(states,a));
    converged = converged + sum(sum(abs(Q-Q_old)));
    % update state
    states = state_next;
    ctr = ctr+1;  
end

end

