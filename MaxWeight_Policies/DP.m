% 参数设置
N = 2;        % 源的数量
L = 5;        % 成功传输次数的阈值
A_max = 10;    % AoI 的最大值
p_succ = 0.5; % 成功传输的概率

% 生成所有可能的状态
a_values = 1:A_max;
c_values = 0:L-1;

% 初始化状态列表和映射
stateID = 0;
stateMap = containers.Map(); % 从状态向量到状态ID的映射
stateList = {}; % 从状态ID到状态向量的映射

% 生成所有状态
[stateList, stateMap, numStates] = generateStates(N, a_values, c_values);

% 动作空间
actions = 1:N; % 动作为要调度的源的索引
numActions = length(actions);

% 预计算转移概率和成本
[P, C] = computeTransitionAndCost(N, L, A_max, p_succ, stateList, stateMap, actions);

% 初始化策略
policy = ones(numStates, 1);

% 策略迭代算法
maxIter = 1000;
epsilon = 1e-6;
avgCost = 0;

for iter = 1:maxIter
    % 策略评估
    [J, avgCost_new] = policyEvaluation(policy, P, C, numStates);

    % 策略改进
    policy_stable = true;
    for sID = 1:numStates
        old_action = policy(sID);
        action_costs = zeros(numActions, 1);
        for a = actions
            % 即时成本
            cost = C(sID, a);
            % 预期下一个状态的值
            nextStates = P{sID, a}.nextStates;
            probs = P{sID, a}.probs;
            expected_J = sum(probs .* J(nextStates));
            action_costs(a) = cost - avgCost_new + expected_J;
        end
        % 选择最小化相对值函数的动作
        min_cost = min(action_costs);
        min_actions = find(abs(action_costs - min_cost) < 1e-6);
        % 如果有多个动作代价相同，选择 AoI 最大的源
        s = stateList{sID};
        aois = s(1:2:end); % 提取 AoI 值
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
        fprintf('策略在第 %d 次迭代后稳定\n', iter);
        break;
    end
    % 更新平均成本
    avgCost = avgCost_new;
end

% 输出最优策略
fprintf('最优调度策略：\n');
for sID = 1:numStates
    s = stateList{sID};
    fprintf('状态 [%s]: 调度源 %d\n', num2str(s), policy(sID));
end

% 可视化策略
visualizePolicy(stateList, policy, N, A_max, L);

% 辅助函数：生成所有状态
function [stateList, stateMap, numStates] = generateStates(N, a_values, c_values)
    % 初始化
    stateList = {};
    stateMap = containers.Map();
    % 为每个源创建状态空间
    args = cell(1, 2*N);
    for i = 1:N
        args{2*i - 1} = a_values;
        args{2*i} = c_values;
    end
    % 使用 allcomb 生成所有可能的状态组合
    stateVectors = allcomb(args{:});
    numStates = size(stateVectors, 1);
    for idx = 1:numStates
        state = stateVectors(idx, :);
        % 确保状态为行向量
        state = state(:)';
        key = stateKey(state);
        if ~isKey(stateMap, key)
            stateID = length(stateList) + 1;
            stateList{stateID} = state;
            stateMap(key) = stateID;
        end
    end
end

% 辅助函数：计算转移概率和成本
function [P, C] = computeTransitionAndCost(N, L, A_max, p_succ, stateList, stateMap, actions)
    numStates = length(stateList);
    P = cell(numStates, length(actions));
    C = zeros(numStates, length(actions));
    for sID = 1:numStates
        s = stateList{sID};
        for a = actions % 动作为要调度的源
            % 初始化
            nextStates = [];
            probs = [];
            % 即时成本为所有源的 AoI 之和
            aois = s(1:2:end);
            C(sID, a) = sum(aois);
            % 成功和失败两种结果
            for outcome = {'success', 'failure'}
                if strcmp(outcome{1}, 'success')
                    prob = p_succ;
                else
                    prob = 1 - p_succ;
                end
                s_next = s;
                % 更新调度源的状态
                i = a;
                a_i = s(2*(i-1)+1);
                c_i = s(2*(i-1)+2);
                if strcmp(outcome{1}, 'success')
                    c_i_new = c_i + 1;
                    if c_i_new == L
                        c_i_new = 0;
                        a_i_new = 1; % AoI 重置为 1
                    else
                        a_i_new = min(a_i + 1, A_max); % AoI 增加 1
                    end
                else
                    c_i_new = c_i; % 失败时，计数器不变
                    a_i_new = min(a_i + 1, A_max); % AoI 增加 1
                end
                s_next(2*(i-1)+1) = a_i_new;
                s_next(2*(i-1)+2) = c_i_new;
                % 更新其他源的 AoI
                for j = 1:N
                    if j ~= i
                        a_j = s(2*(j-1)+1);
                        c_j = s(2*(j-1)+2);
                        a_j_new = min(a_j + 1, A_max); % AoI 增加 1
                        c_j_new = c_j; % 计数器不变
                        s_next(2*(j-1)+1) = a_j_new;
                        s_next(2*(j-1)+2) = c_j_new;
                    end
                end
                % 确保 s_next 为行向量
                s_next = s_next(:)';
                % 获取下一个状态的 ID
                key_next = stateKey(s_next);
                if isKey(stateMap, key_next)
                    s_next_ID = stateMap(key_next);
                else
                    % 如果状态未记录，添加到状态列表和映射中
                    s_next_ID = length(stateList) + 1;
                    stateList{s_next_ID} = s_next;
                    stateMap(key_next) = s_next_ID;
                end
                nextStates = [nextStates; s_next_ID];
                probs = [probs; prob];
            end
            P{sID, a}.nextStates = nextStates;
            P{sID, a}.probs = probs;
        end
    end
end

% 辅助函数：策略评估
function [J, avgCost] = policyEvaluation(policy, P, C, numStates)
    % 初始化线性方程组
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
    % 增加平均成本变量 avgCost，构建增广矩阵
    A_aug = [A, -ones(numStates, 1); ones(1, numStates), 0];
    b_aug = [b; 0];
    % 解线性方程组
    x = A_aug \ b_aug;
    J = x(1:numStates);
    avgCost = x(end);
end

% 辅助函数：生成状态的键值
function key = stateKey(state)
    % 确保状态为行向量
    state = state(:)';
    key = mat2str(state);
end

% 辅助函数：生成所有组合
function combs = allcomb(varargin)
    % 检查输入参数
    narginchk(1, Inf);
    % 确保所有输入都是数值数组
    args = varargin;
    n = nargin;
    for i = 1:n
        if ~isnumeric(args{i})
            error('所有输入参数都应为数值数组。');
        end
    end
    % 使用 ndgrid 生成网格
    [F{1:n}] = ndgrid(args{:});
    % 将结果转换为组合
    for i = n:-1:1
        G(:, i) = F{i}(:);
    end
    combs = G;
end

% 辅助函数：可视化策略
function visualizePolicy(stateList, policy, N, A_max, L)
    % 提取所有状态
    numStates = length(stateList);
    aois = zeros(numStates, N);
    cs = zeros(numStates, N);
    for sID = 1:numStates
        s = stateList{sID};
        aois(sID, :) = s(1:2:end);
        cs(sID, :) = s(2:2:end);
    end

    % 可视化1：固定 AoI，绘制不同 c 值下的策略
    % 假设固定 AoI 为 [a1, a2]
    a_fixed = [3, 6]; % 您可以调整固定的 AoI 值
    idx_a_fixed = (aois(:,1) == a_fixed(1)) & (aois(:,2) == a_fixed(2));
    cs_fixed = cs(idx_a_fixed, :);
    policy_fixed_a = policy(idx_a_fixed);

    figure;
    % 创建网格
    [C1, C2] = meshgrid(0:L-1, 0:L-1);
    Z = zeros(L, L);
    for i = 1:size(cs_fixed,1)
        c1 = cs_fixed(i,1);
        c2 = cs_fixed(i,2);
        action = policy_fixed_a(i);
        Z(c1+1, c2+1) = action;
    end
    imagesc(0:L-1, 0:L-1, Z);
    set(gca,'YDir','normal');
    xlabel('c_1');
    ylabel('c_2');
    title(sprintf('固定 AoI = [%d, %d] 下的策略', a_fixed(1), a_fixed(2)));
    colorbar('Ticks', [1,2], 'TickLabels', {'调度源1','调度源2'});
    colormap([0 0 1; 1 0 0]); % 蓝色和红色

    % 可视化2：固定 c，绘制不同 AoI 值下的策略
    % 假设固定 c 为 [c1, c2]
    c_fixed = [0, 0]; % 您可以调整固定的 c 值
    idx_c_fixed = (cs(:,1) == c_fixed(1)) & (cs(:,2) == c_fixed(2));
    aois_fixed_c = aois(idx_c_fixed, :);
    policy_fixed_c = policy(idx_c_fixed);

    figure;
    % 创建网格
    [A1, A2] = meshgrid(1:A_max, 1:A_max);
    Z = zeros(A_max, A_max);
    for i = 1:size(aois_fixed_c,1)
        a1 = aois_fixed_c(i,1);
        a2 = aois_fixed_c(i,2);
        action = policy_fixed_c(i);
        Z(a1, a2) = action;
    end
    imagesc(1:A_max, 1:A_max, Z);
    set(gca,'YDir','normal');
    xlabel('AoI_1');
    ylabel('AoI_2');
    title(sprintf('固定 c = [%d, %d] 下的策略', c_fixed(1), c_fixed(2)));
    colorbar('Ticks', [1,2], 'TickLabels', {'调度源1','调度源2'});
    colormap([0 0 1; 1 0 0]); % 蓝色和红色
end