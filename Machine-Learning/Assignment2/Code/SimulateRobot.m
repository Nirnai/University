function [ newstate, reward ] = SimulateRobot( state, action, delta, rew )
newstate = delta(state,action);
reward = rew(state,action);
end

