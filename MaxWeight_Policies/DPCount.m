% Parameter settings
N = 2;              
A_max = 15;           % Maximum value of AoI
L_i = [6, 2];        
p_i = [1, 1];     
    
% Initialize state list and mapping
stateID = 0;
stateMap = containers.Map();
stateList = {};

% Initialize initial state
initialState = zeros(1, 3*N);
for i = 1:N
    initialState(3*(i-1)+1) = 1; % AoI
    initialState(3*(i-1)+2) = 0; % c_i
    initialState(3*(i-1)+3) = 0; % d_i
end

stateKeyStr = stateKey(initialState);
stateID = stateID + 1;
stateMap(stateKeyStr) = stateID;
stateList{stateID} = initialState;

% Action space
actions = 1:N; % Actions are the indices of sources to schedule
numActions = length(actions);

% Precompute transition probabilities and costs
[P, C, stateList, stateMap] = computeTransitionAndCost(N, L_i, A_max, p_i, stateList, stateMap, actions);

% Update number of states
numStates = length(stateList);

% Initialize policy
policy = ones(numStates, 1);

% Policy iteration algorithm
maxIter = 1000;
epsilon = 1e-6;
avgCost = 0;

for iter = 1:maxIter
    % Policy evaluation
    [J, avgCost_new] = policyEvaluation(policy, P, C, numStates);

    % Policy improvement
    policy_stable = true;
    for sID = 1:numStates
        old_action = policy(sID);
        action_costs = zeros(numActions, 1);
        for a = actions
            % Immediate cost
            cost = C(sID, a);
            % Expected value of next state
            nextStates = P{sID, a}.nextStates;
            probs = P{sID, a}.probs;
            expected_J = sum(probs .* J(nextStates));
            action_costs(a) = cost - avgCost_new + expected_J;
        end
        % Choose action that minimizes the relative value function
        min_cost = min(action_costs);
        min_actions = find(abs(action_costs - min_cost) < 1e-6);
        % If multiple actions have the same cost, choose the source with the maximum AoI
        s = stateList{sID};
        aois = s(1:3:end); % Extract AoI values
        if length(min_actions) > 1
            [~, idx] = max(aois(min_actions));
            min_action = min_actions(idx);
        else
            min_action = min_actions(1);
        end
        policy(sID) = min_action;
        if old_action ~= min_action
            policy_stable = false;
        end
    end
    if policy_stable
        fprintf('Policy is stable after %d iterations\n', iter);
        break;
    end
    % Update average cost
    avgCost = avgCost_new;
end

% Output optimal policy (optional)
% fprintf('Optimal scheduling policy:\n');
% for sID = 1:numStates
%     s = stateList{sID};
%     fprintf('State [%s]: Schedule source %d\n', num2str(s), policy(sID));
% end

% Perform multiple simulation experiments to verify if packet switching exists
numSimulations = 100; % Number of simulations
simulationLength = 10000; % Number of time steps per simulation

packetSwitchingCount = 0; % Record the number of times packet switching occurs
totalSwitchStatesMap = containers.Map(); % Map to record states and switch counts

for sim = 1:numSimulations
    [switchOccurred, switchStatesMap] = simulateSystem(N, L_i, p_i, A_max, policy, stateList, stateMap, initialState, simulationLength);
    if switchOccurred
        packetSwitchingCount = packetSwitchingCount + 1;
    end
    % Merge switchStatesMap into totalSwitchStatesMap
    keys = switchStatesMap.keys;
    for k = 1:length(keys)
        key = keys{k};
        count = switchStatesMap(key);
        if isKey(totalSwitchStatesMap, key)
            totalSwitchStatesMap(key) = totalSwitchStatesMap(key) + count;
        else
            totalSwitchStatesMap(key) = count;
        end
    end
end

% Output results
fprintf('In %d simulations, packet switching occurred %d times.\n', numSimulations, packetSwitchingCount);
if packetSwitchingCount > 0
    fprintf('Packet switching occurs when executing the policy.\n');
    % Output the states where switching occurred and counts
    disp('States where switching occurred and the number of times it switched in that state:');
    switchStatesKeys = totalSwitchStatesMap.keys;
    for i = 1:length(switchStatesKeys)
        stateKeyStr = switchStatesKeys{i};
        count = totalSwitchStatesMap(stateKeyStr);
        stateVec = str2num(stateKeyStr);
        fprintf('State [%s]: Switched %d times\n', num2str(stateVec), count);
    end
else
    fprintf('No packet switching occurs when executing the policy.\n');
end

% Auxiliary function: simulate system operation
function [switchOccurred, switchStatesMap] = simulateSystem(N, L_i, p_i, A_max, policy, stateList, stateMap, initialState, simulationLength)
    currentState = initialState;
    sKey = stateKey(currentState);
    sID = stateMap(sKey);
    switchOccurred = false;
    ongoingSource = 0; % Currently transmitting source (cumulative successful transmissions not reaching L_i)
    switchStatesMap = containers.Map();

    for t = 1:simulationLength
        % Get action for the current state
        action = policy(sID);

        % Check if packet switching occurs
        if ongoingSource == 0
            % Start a new transmission
            ongoingSource = action;
        elseif action ~= ongoingSource
            % If the source is switched before transmission is completed, record packet switching
            switchOccurred = true;
            % Record the current state
            sKey = stateKey(currentState);
            if isKey(switchStatesMap, sKey)
                switchStatesMap(sKey) = switchStatesMap(sKey) + 1;
            else
                switchStatesMap(sKey) = 1;
            end
            % Update ongoingSource to the new action
            ongoingSource = action;
            % Continue simulation to collect more data
        end

        % Determine transmission result based on success probability
        p_succ = p_i(action);
        if rand() < p_succ
            outcome = 'success';
        else
            outcome = 'failure';
        end

        % Update state
        [nextState] = updateState(currentState, action, outcome, N, L_i, A_max);
        % Update ongoingSource
        idx_c = 3*(ongoingSource-1)+2;
        c_i = nextState(idx_c);
        if c_i == 0
            % Transmission of the current source is completed
            ongoingSource = 0;
        end

        % Update state ID
        sKey = stateKey(nextState);
        if isKey(stateMap, sKey)
            sID = stateMap(sKey);
        else
            % If the new state is not in stateMap, add it
            sID = length(stateList) + 1;
            stateList{sID} = nextState;
            stateMap(sKey) = sID;
        end

        currentState = nextState;
    end
end

% Auxiliary function: update state
function [s_next] = updateState(s, a, outcome, N, L_i, A_max)
    s_next = s;
    % Update the state of the scheduled source
    i = a;
    idx_a = 3*(i-1)+1;
    idx_c = idx_a + 1;
    idx_d = idx_a + 2;
    a_i = s(idx_a);
    c_i = s(idx_c);
    d_i = s(idx_d);

    if strcmp(outcome, 'success')
        c_i_new = c_i + 1;
        if c_i == 0 % First successful transmission, start timing
            d_i_new = 1;
        else
            d_i_new = min(d_i + 1, A_max); % Increment d_i, limit to maximum value
        end
        if c_i_new == L_i(i)
            % Threshold reached, update AoI
            a_i_new = min(d_i_new, A_max); % Update AoI, limit to maximum value
            c_i_new = 0;
            d_i_new = 0;
        else
            a_i_new = min(a_i + 1, A_max);
        end
    else % Transmission failure
        a_i_new = min(a_i + 1, A_max);
        c_i_new = c_i;
        if c_i > 0
            d_i_new = min(d_i + 1, A_max); % Increment d_i, limit to maximum value
        else
            d_i_new = 0;
        end
    end
    s_next(idx_a) = a_i_new;
    s_next(idx_c) = c_i_new;
    s_next(idx_d) = d_i_new;

    % Update AoI for other sources
    for j = 1:N
        if j ~= i
            idx_a_j = 3*(j-1)+1;
            idx_c_j = idx_a_j + 1;
            idx_d_j = idx_a_j + 2;
            a_j = s(idx_a_j);
            c_j = s(idx_c_j);
            d_j = s(idx_d_j);
            a_j_new = min(a_j + 1, A_max);
            c_j_new = c_j;
            if c_j > 0
                d_j_new = min(d_j + 1, A_max); % Increment d_j, limit to maximum value
            else
                d_j_new = 0;
            end
            s_next(idx_a_j) = a_j_new;
            s_next(idx_c_j) = c_j_new;
            s_next(idx_d_j) = d_j_new;
        end
    end
end

% Auxiliary function: compute transition probabilities and costs
function [P, C, stateList, stateMap] = computeTransitionAndCost(N, L_i, A_max, p_i, stateList, stateMap, actions)
    numStates = length(stateList);
    P = cell(numStates, N);
    C = zeros(numStates, N);
    queue = 1:numStates; % Queue for breadth-first search

    processed = false(1, numStates);

    while ~isempty(queue)
        sID = queue(1);
        queue(1) = [];
        if processed(sID)
            continue;
        end
        processed(sID) = true;

        s = stateList{sID};
        for a = actions % Action is the source to schedule
            % Immediate cost is the sum of AoI of all sources
            aois = s(1:3:end);
            C(sID, a) = sum(aois);
            nextStates = [];
            probs = [];
            % Probability of transmission success and failure
            p_succ = p_i(a);
            p_fail = 1 - p_succ;
            for outcome = {'success', 'failure'}
                if strcmp(outcome{1}, 'success')
                    prob = p_succ;
                else
                    prob = p_fail;
                end
                s_next = s;
                % Update the state of the scheduled source
                i = a;
                idx_a = 3*(i-1)+1;
                idx_c = idx_a + 1;
                idx_d = idx_a + 2;
                a_i = s(idx_a);
                c_i = s(idx_c);
                d_i = s(idx_d);
                if strcmp(outcome{1}, 'success')
                    c_i_new = c_i + 1;
                    if c_i == 0 % First successful transmission, start timing
                        d_i_new = 1;
                    else
                        d_i_new = min(d_i + 1, A_max); % Increment d_i, limit to maximum value
                    end
                    if c_i_new == L_i(i)
                        % Threshold reached, update AoI
                        a_i_new = min(d_i_new, A_max); % Update AoI, limit to maximum value
                        c_i_new = 0;
                        d_i_new = 0;
                    else
                        a_i_new = min(a_i + 1, A_max);
                    end
                else % Transmission failure
                    a_i_new = min(a_i + 1, A_max);
                    c_i_new = c_i;
                    if c_i > 0
                        d_i_new = min(d_i + 1, A_max); % Increment d_i, limit to maximum value
                    else
                        d_i_new = 0;
                    end
                end
                s_next(idx_a) = a_i_new;
                s_next(idx_c) = c_i_new;
                s_next(idx_d) = d_i_new;
                % Update AoI for other sources
                for j = 1:N
                    if j ~= i
                        idx_a_j = 3*(j-1)+1;
                        idx_c_j = idx_a_j + 1;
                        idx_d_j = idx_a_j + 2;
                        a_j = s(idx_a_j);
                        c_j = s(idx_c_j);
                        d_j = s(idx_d_j);
                        a_j_new = min(a_j + 1, A_max);
                        c_j_new = c_j;
                        if c_j > 0
                            d_j_new = min(d_j + 1, A_max); % Increment d_j, limit to maximum value
                        else
                            d_j_new = 0;
                        end
                        s_next(idx_a_j) = a_j_new;
                        s_next(idx_c_j) = c_j_new;
                        s_next(idx_d_j) = d_j_new;
                    end
                end
                % Ensure s_next is a row vector
                s_next = s_next(:)';
                % Get the ID of the next state
                key_next = stateKey(s_next);
                if isKey(stateMap, key_next)
                    s_next_ID = stateMap(key_next);
                else
                    % If state not recorded, add to state list and mapping
                    s_next_ID = length(stateList) + 1;
                    stateList{s_next_ID} = s_next;
                    stateMap(key_next) = s_next_ID;
                    % Add the new state to the queue
                    queue(end+1) = s_next_ID;
                    processed(s_next_ID) = false;
                    % Expand C and P matrices
                    C(s_next_ID, N) = 0;
                    P{s_next_ID, N} = [];
                end
                nextStates = [nextStates; s_next_ID];
                probs = [probs; prob];
            end
            P{sID, a}.nextStates = nextStates;
            P{sID, a}.probs = probs;
        end
    end
end

% Auxiliary function: policy evaluation
function [J, avgCost] = policyEvaluation(policy, P, C, numStates)
    % Initialize linear equation system
    A = zeros(numStates, numStates);
    b = zeros(numStates, 1);
    for sID = 1:numStates
        a = policy(sID);
        cost = C(sID, a);
        nextStates = P{sID, a}.nextStates;
        probs = P{sID, a}.probs;
        A(sID, sID) = 1;
        for idx = 1:length(nextStates)
            s_next_ID = nextStates(idx);
            prob = probs(idx);
            A(sID, s_next_ID) = A(sID, s_next_ID) - prob;
        end
        b(sID) = cost;
    end
    % Add average cost variable avgCost and construct augmented matrix
    A_aug = [A, -ones(numStates, 1); ones(1, numStates), 0];
    b_aug = [b; 0];
    % Solve the linear equation system
    x = A_aug \ b_aug;
    J = x(1:numStates);
    avgCost = x(end);
end

% Auxiliary function: generate state key
function key = stateKey(state)
    % Ensure state is a row vector
    state = state(:)';
    key = mat2str(state);
end