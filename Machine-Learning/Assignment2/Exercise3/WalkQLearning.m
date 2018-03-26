function [ ] = WalkQLearning(s)
%%
rew = zeros(16,4);

rew(1,:) = [0 -1 0 -1];
rew(2,:) = [0 0 -1 -1];
rew(3,:) = [0 0 -1 -1];
rew(4,:) = [-1 -1 0 -1];

rew(5,:) = [-1 -1 -1 0];
% rew(6,:) = [-1 1 1 -1];
% rew(7,:) = [1 -1 1 -1];
rew(8,:) = [0 1 0 0];

rew(9,:) = [-1 -1 0 -1];
% rew(10,:) = [-1 1 -1 -1];
% rew(11,:) = [1 -1 -1 -1];
rew(12,:) = [0 1 0 -1];

rew(13,:) = [0 -1 0 -1];
rew(14,:) = [-1 0 0 1];
rew(15,:) = [-1 -1 0 1];
rew(16,:) = [0 -1 0 -1];


    delta = [2 4 5 13; ...
             1 3 6 14; ...
             4 2 7 15; ...
             3 1 8 16; ...
             6 8 1 9;  ...
             5 7 2 10; ...
             8 6 3 11; ...
             7 5 4 12; ...
             10 12 13 5; ...
             9 11 14 6; ...
             12 10 15 7; ...
             11 9 16 8; ...
             14 16 9 1; ...
             13 15 10 2; ...
             16 14 11 3; ...
             15 13 12 4];


%% Q-Learning
Q = zeros(16,4); 
eps = 0.01;
initstate = s;
gamma = 0.99;
alpha = 0.1;

[ initstate, Q, delta, rew ] = QLearning( delta, rew, Q, eps, initstate, gamma, alpha, 100000 );

%% testing with greedy
[ state_sequence ] = test_greedy( initstate, Q, delta, rew);
walkshow(state_sequence');
end