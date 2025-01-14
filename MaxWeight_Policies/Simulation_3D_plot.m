% function Simulation_3D_plot()
%     clear
%     clc
% 
%     global probability cost alpha arrival optQ optMU optMU_buff optMU_NO_buff L
% 
%     tic;
% 
%     %% 固定参数设置
%     M = 2;                       % 源的个数
%     p_fixed = 0.5;               % 信道可靠度 p=0.5
%     arrival = [1; 1];            % 到达率(或其他场景相关参数)，保持原值
%     num_iterations = 10;         % 蒙特卡洛仿真迭代次数(可根据需要调整)
%     T = 1;                       % 时隙长度(保留原先的设定)
%     h1 = 1;                      % 初始 age
%     z1 = 0;                      % 初始队列(或其他状态)，需 < h1
% 
%     % x轴: L_1 ∈ [1,30]
%     L1_values = 2:2:30;     
%     % y轴: alpha_1 ∈ [1,10]
%     alpha1_values = 1:10;
% 
%     % 分别用于保存四个策略下的 AoI 结果
%     % 这里我们称为 Z_XXX(i, j) = 在 L1_values(i), alpha1_values(j) 组合时的 平均AoI
%     Z_rand_switch         = zeros(length(L1_values), length(alpha1_values));
%     Z_rand_no_switch      = zeros(length(L1_values), length(alpha1_values));
%     Z_rand_no_switch_sim  = zeros(length(L1_values), length(alpha1_values));
%     Z_lower_bound         = zeros(length(L1_values), length(alpha1_values));
% 
%     %% 主循环
%     for i = 1:length(L1_values)
%         for j = 1:length(alpha1_values)
% 
%             % 赋值 L, alpha
%             L = [L1_values(i); 2];          % L_1, L_2=5
%             alpha = [alpha1_values(j); 10];  % alpha_1, alpha_2=5
%             probability = p_fixed * [1; 1]; % 信道可靠度 p=0.5
% 
%             cost = zeros(M,1);
% 
%             % ========== 1. 下界 (Lower Bound) ==========
%             Lower_Bound_find_q(M);   
%             % 存储下界
%             Z_lower_bound(i,j) = sum(alpha.*(1+L./optQ)) / (2*M);
% 
%             % ========== 2. Switching Randomized Policy ==========
%             % 先计算 E[Y^2] (在原代码中用于随机策略相关分析)
%             E_Y2 = zeros(M,1);
%             for idx = 1:M
%                 E_Y2(idx) = 1 + 2*(L(idx)-1) + ...
%                             ((L(idx)-1)*(L(idx)-probability(idx))) / probability(idx);
%             end
% 
%             % 更新 optMU_buff (与原代码一致)
%             Random_find_mu_buffered_BerBer1(M);
% 
%             % 原代码给出解析公式:
%             % (sum(sqrt(alpha./p.*(3*L.^2 - L)./(2.*L))))^2 / M + sum(alpha.*3/2) / M
%             sum_term = sum( sqrt( alpha ./ probability .* ...
%                                   (3*L.^2 - L) ./ (2.*L) ) );
%             Z_rand_switch(i,j) = (sum_term)^2 / M + sum( alpha .* 3/2 ) / M;
% 
%             % ========== 3. No-Switching Random Policy ==========
%             % （1）解析解
%             mu_opt_no_switch = optimize_mu_for_random_policy(M, probability, L, alpha, E_Y2);
%             Z_rand_no_switch(i,j) = ...
%                 weighted_AoI(mu_opt_no_switch, probability, L, E_Y2, alpha) - alpha1_values(j)/2-0.5 ;
% 
%             % （2）仿真
%             K = 5000* M;  % 仿真时隙数(可根据需求调整)
%             [~, cost_no_switch_sim] = ...
%                 No_Switch_Random_simulation_stochastic(K, T, M, h1, z1, ...
%                                                        num_iterations, mu_opt_no_switch);
%             Z_rand_no_switch_sim(i,j) = cost_no_switch_sim;
%         end
%     end
% 
%     %% 绘图：四个曲面叠加在一张 3D 图中
%     % x轴 = L_1, y轴 = alpha_1, z轴 = 对应策略的 AoI
%     [X, Y] = meshgrid(L1_values, alpha1_values);
%     % 注意: Z(i,j) 对应 X(i), Y(j)。若 Z_rand_switch(i,j) 是 i->L1, j->alpha1，
%     % 那么要么对 Z 矩阵做转置，要么在 surf() 时使用转置矩阵。
% 
%     figure('Name','3D Comparison of Policies','NumberTitle','off');
%     hold on;
% 
%     % 让曲面带点透明度，便于相互对比
%     s1 = surf(X, Y, Z_rand_switch', ...
%         'FaceAlpha',0.6, 'EdgeColor','none', 'FaceColor',[0.1 0.7 0.9]); 
%     s2 = surf(X, Y, Z_rand_no_switch', ...
%         'FaceAlpha',0.6, 'EdgeColor','none', 'FaceColor',[0.9 0.1 0.5]);
%     s3 = surf(X, Y, Z_rand_no_switch_sim', ...
%         'FaceAlpha',0.6, 'EdgeColor','none', 'FaceColor',[0.1 0.9 0.1]);
%     s4 = surf(X, Y, Z_lower_bound', ...
%         'FaceAlpha',0.6, 'EdgeColor','none', 'FaceColor',[0.9 0.7 0.1]);
% 
%     xlabel('L_1');
%     ylabel('\alpha_1');
%     zlabel('Expected Weighted Sum AoI');
%     legend([s1,s2,s3,s4], ...
%         {'Switching Random','No-Switch (Analytic)', ...
%          'No-Switch (Simulation)','Lower Bound'}, ...
%         'Location','northeast');
%     title('3D Comparison: Four Policies');
%     view(45, 45);  % 调整视角 (可自行调整角度)
% 
%     hold off;
%     toc;
% end


function Simulation_3D_plot()
    clear
    clc

    global probability cost alpha arrival optQ optMU optMU_buff optMU_NO_buff L

    tic;

    %% 固定参数设置
    M = 2;                       % 源的个数
    p_fixed = 0.5;               % 信道可靠度 p=0.5
    arrival = [1; 1];            % 到达率(或其他场景相关参数)，保持原值
    num_iterations = 10;         % 蒙特卡洛仿真迭代次数(可根据需要调整)
    T = 1;                       % 时隙长度(保留原先的设定)
    h1 = 1;                      % 初始 age
    z1 = 0;                      % 初始队列(或其他状态)，需 < h1
    
    % x轴: L_1 ∈ [2,4,6,...,30]
    L1_values = 2:2:30;     
    % y轴: alpha_1 ∈ [1,2,3,...,10]
    alpha1_values = 1:10;
    
    % 分别用于保存三个策略下的 AoI 结果
    % 这里我们称为 Z_XXX(i, j) = 在 L1_values(i), alpha1_values(j) 组合时的 平均AoI
    Z_rand_switch    = zeros(length(L1_values), length(alpha1_values));  % Switching Random
    Z_rand_no_switch = zeros(length(L1_values), length(alpha1_values));  % No-Switching Random(解析)
    Z_lower_bound    = zeros(length(L1_values), length(alpha1_values));  % Lower Bound
    
    %% 主循环
    for i = 1:length(L1_values)
        for j = 1:length(alpha1_values)
            
            % 赋值 L, alpha
            L = [L1_values(i); 2];           % L_1, L_2=2 (原文你写了, 改成你需要的值即可)
            alpha = [alpha1_values(j); 10];  % alpha_1, alpha_2=10 (原文写了10)
            probability = p_fixed * [1; 1];  % 信道可靠度 p=0.5
            
            cost = zeros(M,1);
            
            % ========== 1. 下界 (Lower Bound) ==========
            Lower_Bound_find_q(M);   
            Z_lower_bound(i,j) = sum(alpha .* (1 + L./optQ)) / (2*M);
            
            % ========== 2. Switching Randomized Policy ==========
            % 先计算 E[Y^2] (在原代码中用于随机策略相关分析)
            E_Y2 = zeros(M,1);
            for idx = 1:M
                E_Y2(idx) = 1 + 2*(L(idx)-1) + ...
                            ((L(idx)-1)*(L(idx)-probability(idx))) / probability(idx);
            end
            
            % 更新 optMU_buff (与原代码一致)
            Random_find_mu_buffered_BerBer1(M);
            
            % 原代码的解析公式:
            % (sum(sqrt(alpha./p.*(3*L.^2 - L)./(2.*L))))^2 / M + sum(alpha.*3/2) / M
            sum_term = sum( sqrt( alpha ./ probability .* ...
                                  (3*L.^2 - L) ./ (2.*L) ) );
            Z_rand_switch(i,j) = (sum_term)^2 / M + sum(alpha .* 3/2) / M;
            
            % ========== 3. No-Switching Random Policy (解析) ==========
            mu_opt_no_switch = optimize_mu_for_random_policy(M, probability, L, alpha, E_Y2);
            % 此处做了一点修正，原来代码中: 
            %   Z_rand_no_switch(i,j) = weighted_AoI(...) - alpha1_values(j)/2 - 0.5
            % 如果你不想做这个额外减法，可以去掉
            Z_rand_no_switch(i,j) = weighted_AoI(mu_opt_no_switch, probability, L, E_Y2, alpha);
            
            % --------------- 已删除「No-Switching Random Policy 的仿真」部分 ---------------
        end
    end

    %% 绘图：三个曲面叠加在一张 3D 图中
    % x轴 = L_1, y轴 = alpha_1, z轴 = 对应策略的 AoI
    [X, Y] = meshgrid(L1_values, alpha1_values);
    
    figure('Name','3D Comparison of Policies','NumberTitle','off');
    hold on;
    
    % 让曲面带点透明度，便于相互对比
    % 1) Switching Random
    s1 = surf(X, Y, Z_rand_switch', ...
        'FaceAlpha',0.6, 'EdgeColor','none', 'FaceColor',[0.1 0.7 0.9]); 
    % 2) No-Switching Random (解析)
    s2 = surf(X, Y, Z_rand_no_switch', ...
        'FaceAlpha',0.6, 'EdgeColor','none', 'FaceColor',[0.9 0.1 0.5]);
    % 3) Lower Bound
    s3 = surf(X, Y, Z_lower_bound', ...
        'FaceAlpha',0.6, 'EdgeColor','none', 'FaceColor',[0.9 0.7 0.1]);
    
    xlabel('L_1');
    ylabel('\alpha_1');
    zlabel('Expected Weighted Sum AoI');
    legend([s1, s2, s3], ...
        {'Switching Random','No-Switching Random (Analytic)','Lower Bound'}, ...
        'Location','northeast');
    title('3D Comparison: Three Policies');
    view(45, 45);  % 调整视角 (可自行调整角度)
    
    hold off;
    toc;
end