function [ state_sequence ] = test_greedy( initstate, Q, delta, rew)

ctr = 1;
state = initstate;
state_sequence= zeros(16,1);
while(ctr < 17)
    [~, a] = max(Q(state,:));
    [s_next,~] = SimulateRobot(state,a,delta,rew);
    state_sequence(ctr) = state;
    state= s_next;
    ctr = ctr+1;
end

end

