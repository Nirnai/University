function [ ] = WalkPolicyIteration(s)
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
    
    %% Policy Iteration
    policy = ceil(rand(16,1)*4);
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



    %%
    gamma = 0.99;
    converged = 1;
    while(converged ~= 0)
        temp = policy;
        % a)
        r = zeros(16,1);
        A = zeros(16,16);
        for state = 1:16
            r(state,1) = rew(state, policy(state,1));
            A(state,delta(state, policy(state,1))) = 1;
        end

        % LSG
        value = -inv((gamma*A-eye(16))) * r;

        % b)
        Vp = zeros(16,4);
        for state = 1:16
            for action = 1:4
                Vp(state,action) = rew(state,action) + gamma * value(delta(state,action));
            end
        end
        [~,policy] = max(Vp,[],2);
        converged = sum(abs(temp - policy));
        %disp(converged);
    end


    %% Test Policy
    state_sequence = zeros(1,16);
    state_sequence(1) = s;
    for i = 2:16
        state_sequence(i) = delta(state_sequence(i-1), policy(state_sequence(i-1)));
    end

    walkshow(state_sequence);
end