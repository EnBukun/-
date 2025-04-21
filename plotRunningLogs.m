f1 = figure(1);
f2 = figure(2);
clf (f1, "reset")  % 重置图窗1
clf (f2, "reset")   % 重置图窗2
clc                 % 清空命令窗口
clear               % 清空工作区变量

% --------------- 加载五次运行的数据文件 ---------------
first_run_distances = load('workingDirectory/first_run_distance.txt');
first_run_thetas = load('workingDirectory/first_run_theta.txt');
first_run_sideline_distances = load('workingDirectory/first_run_side.txt');
first_run_crossline_distances = load('workingDirectory/first_run_cross.txt');
first_run_current_velocity = load('workingDirectory/first_run_current_velocity.txt');
first_run_target_velocity = load('workingDirectory/first_run_target_velocity.txt');

second_run_distances = load('workingDirectory/second_run_distance.txt');
second_run_thetas = load('workingDirectory/second_run_theta.txt');
second_run_sideline_distances = load('workingDirectory/second_run_side.txt');
second_run_crossline_distances = load('workingDirectory/second_run_cross.txt');
second_run_current_velocity = load('workingDirectory/second_run_current_velocity.txt');
second_run_target_velocity = load('workingDirectory/second_run_target_velocity.txt');

third_run_distances = load('workingDirectory/third_run_distance.txt');
third_run_thetas = load('workingDirectory/third_run_theta.txt');
third_run_sideline_distances = load('workingDirectory/third_run_side.txt');
third_run_crossline_distances = load('workingDirectory/third_run_cross.txt');
third_run_current_velocity = load('workingDirectory/third_run_current_velocity.txt');
third_run_target_velocity = load('workingDirectory/third_run_target_velocity.txt');

fourth_run_distances = load('workingDirectory/fourth_run_distance.txt');
fourth_run_thetas = load('workingDirectory/fourth_run_theta.txt');
fourth_run_sideline_distances = load('workingDirectory/fourth_run_side.txt');
fourth_run_crossline_distances = load('workingDirectory/fourth_run_cross.txt');
fourth_run_current_velocity = load('workingDirectory/fourth_run_current_velocity.txt');
fourth_run_target_velocity = load('workingDirectory/fourth_run_target_velocity.txt');

fifth_run_distances = load('workingDirectory/fifth_run_distance.txt');
fifth_run_thetas = load('workingDirectory/fifth_run_theta.txt');
fifth_run_sideline_distances = load('workingDirectory/fifth_run_side.txt');
fifth_run_crossline_distances = load('workingDirectory/fifth_run_cross.txt');
fifth_run_current_velocity = load('workingDirectory/fifth_run_current_velocity.txt');
fifth_run_target_velocity = load('workingDirectory/fifth_run_target_velocity.txt');

% --------------- 数据预处理 ---------------
% 去除零值并调整角度数据
first_run_distances = nonzeros(first_run_distances);        % 毫米单位
first_run_thetas = first_run_thetas(1 : size(first_run_distances)); % 弧度单位
first_run_sideline_distances = nonzeros(first_run_sideline_distances);
first_run_crossline_distances = nonzeros(first_run_crossline_distances);
first_run_current_velocity = nonzeros(first_run_current_velocity); % 米/秒单位
first_run_target_velocity = first_run_target_velocity(1 : size(first_run_current_velocity));
first_run_thetas = first_run_thetas * 1.015;  % 角度校准系数

second_run_distances = nonzeros(second_run_distances);
second_run_thetas = second_run_thetas(1 : size(second_run_distances));
second_run_sideline_distances = nonzeros(second_run_sideline_distances);
second_run_crossline_distances = nonzeros(second_run_crossline_distances);
second_run_current_velocity = nonzeros(second_run_current_velocity);
second_run_target_velocity = second_run_target_velocity(1 : size(second_run_current_velocity));
second_run_thetas = second_run_thetas * 1.015;

third_run_distances = nonzeros(third_run_distances);
third_run_thetas = third_run_thetas(1 : size(third_run_distances));
third_run_sideline_distances = nonzeros(third_run_sideline_distances);
third_run_crossline_distances = nonzeros(third_run_crossline_distances);
third_run_current_velocity = nonzeros(third_run_current_velocity);
third_run_target_velocity = third_run_target_velocity(1 : size(third_run_current_velocity));
third_run_thetas = third_run_thetas * 1.015;

fourth_run_distances = nonzeros(fourth_run_distances);
fourth_run_thetas = fourth_run_thetas(1 : size(fourth_run_distances));
fourth_run_sideline_distances = nonzeros(fourth_run_sideline_distances);
fourth_run_crossline_distances = nonzeros(fourth_run_crossline_distances);
fourth_run_current_velocity = nonzeros(fourth_run_current_velocity);
fourth_run_target_velocity = fourth_run_target_velocity(1 : size(fourth_run_current_velocity));
fourth_run_thetas = fourth_run_thetas * 1.015;

fifth_run_distances = nonzeros(fifth_run_distances);
fifth_run_thetas = fifth_run_thetas(1 : size(fifth_run_distances));
fifth_run_sideline_distances = nonzeros(fifth_run_sideline_distances);
fifth_run_crossline_distances = nonzeros(fifth_run_crossline_distances);
fifth_run_current_velocity = nonzeros(fifth_run_current_velocity);
fifth_run_target_velocity = fifth_run_target_velocity(1 : size(fifth_run_current_velocity));
fifth_run_thetas = fifth_run_thetas * 1.015;

% --------------- 可视化模块 ---------------
% 图1：首航次路径信息
figure(1);
plotCourseInformation(first_run_distances, first_run_thetas,...
    first_run_sideline_distances, first_run_crossline_distances);
title('1走目');

% 图2：后续航次路径对比
figure(2)
subplot(2, 2, 1);
plotCourseInformation(second_run_distances, second_run_thetas,...
    second_run_sideline_distances, second_run_crossline_distances);
title('2走目')

subplot(2, 2, 2);
plotCourseInformation(third_run_distances, third_run_thetas,...
    third_run_sideline_distances, third_run_crossline_distances);
title('3走目')

subplot(2, 2, 3);
plotCourseInformation(fourth_run_distances, fourth_run_thetas,...
    fourth_run_sideline_distances, fourth_run_crossline_distances);
title('4走目')

subplot(2, 2, 4);
plotCourseInformation(fifth_run_distances, fifth_run_thetas,...
    fifth_run_sideline_distances, fifth_run_crossline_distances);
title('5走目')

% 图3：首航次运动学分析
figure(3);
subplot(2, 1, 1);
plotRadius(first_run_distances, first_run_thetas);  % 转弯半径分析
title('转弯半径')

subplot(2, 1, 2);
plotdTheta(first_run_distances, first_run_thetas);  % 角速度分析
title('角速度变化')

% 图4：速度规划策略对比
figure(4)
subplot(2, 2, 1);
plotVelocityTable(first_run_distances, first_run_thetas, 6.0, 2.5, 1000, 'linear');
title('线性速度规划')

subplot(2, 2, 2);
plotVelocityTable(first_run_distances, first_run_thetas, 6.0, 2.5, 1000, 'sigmoid');
title('S型速度规划')

subplot(2, 2, 3);
plotVelocityTable(first_run_distances, first_run_thetas, 6.0, 2.5, 1000, 'sigmoid');
title('优化S型速度规划')

subplot(2, 2, 4);
plotVelocityTable(first_run_distances, first_run_thetas, 6.5, 2.5, 1000, 'sigmoid');
title('调参后速度规划')

% 图8-9：速度跟踪性能
figure(8)
plot(1:length(first_run_current_velocity), first_run_current_velocity)
hold on
plot(1:length(first_run_target_velocity), first_run_target_velocity)
hold off
legend('实际速度', '目标速度')
title('首航次速度跟踪')

figure(9)
subplot(2, 2, 1);
plot(1:length(second_run_current_velocity), second_run_current_velocity)
hold on
plot(1:length(second_run_target_velocity), second_run_target_velocity)
hold off
legend('实际速度', '目标速度')
title('第二航次速度跟踪')

subplot(2, 2, 2);
plot(1:length(third_run_current_velocity), third_run_current_velocity)
hold on
plot(1:length(third_run_target_velocity), third_run_target_velocity)
hold off
legend('实际速度', '目标速度')
title('第三航次速度跟踪')

subplot(2, 2, 3);
plot(1:length(fourth_run_current_velocity), fourth_run_current_velocity)
hold on
plot(1:length(fourth_run_target_velocity), fourth_run_target_velocity)
hold off
legend('实际速度', '目标速度')
title('第四航次速度跟踪')

subplot(2, 2, 4);
plot(1:length(fifth_run_current_velocity), fifth_run_current_velocity)
hold on
plot(1:length(fifth_run_target_velocity), fifth_run_target_velocity)
hold off
legend('实际速度', '目标速度')
title('第五航次速度跟踪')

% --------------- 自定义函数定义 ---------------
function plotCourseInformation(distances, thetas, sidelines, crosslines)
% 绘制完整的路径信息图，包含：
% - 机器人中心轨迹（蓝色散点）
% - 传感器轨迹（青色散点）
% - 边线标记（品红色星号）
% - 跨线标记（红色叉号）
    
    % 初始化位置参数
    x = 0;  % 初始x坐标(mm)
    y = 0;  % 初始y坐标(mm)
    th = 0; % 初始航向角(rad)
    total_distance = 0;  % 累计行驶距离
    
    % 预分配内存
    robot_positions = zeros(length(distances), 2);  % 机器人中心轨迹
    sensor_positions = zeros(length(distances), 2);  % 传感器轨迹
    sideline_positions = zeros(length(sidelines), 2); % 边线标记位置
    crossline_positions = zeros(length(crosslines), 2); % 跨线标记位置
    
    pivot_length = 110;  % 轴距(mm)
    sideline_idx = 1;    % 边线标记计数器
    crossline_idx = 1;   % 跨线标记计数器

    % 轨迹计算循环
    for i = 1:length(distances)
        % 使用圆弧运动模型更新位置
        delta_theta = thetas(i);
        chord_length = distances(i);
        
        % 更新中心坐标
        x = x + chord_length * cos(th + delta_theta/2);
        y = y + sin(th + delta_theta/2) * chord_length;
        th = th + delta_theta;
        
        % 存储轨迹点
        robot_positions(i,:) = [x, y];
        sensor_positions(i,:) = [x + pivot_length*cos(th), y + pivot_length*sin(th)];
        
        % 边线标记检测
        if ~isempty(sidelines)
            if sideline_idx <= length(sidelines) &&...
               abs(total_distance - sidelines(sideline_idx)) < 10
                sideline_positions(sideline_idx,:) = [x, y + 100*cos(th)];
                sideline_idx = sideline_idx + 1;
            end
        end
        
        % 跨线标记检测
        if ~isempty(crosslines)
            if crossline_idx <= length(crosslines) &&...
               abs(total_distance - crosslines(crossline_idx)) < 10
                crossline_positions(crossline_idx,:) = [x, y];
                crossline_idx = crossline_idx + 1;
            end
        end
        
        total_distance = total_distance + distances(i);
    end

    % 可视化绘制
    hold on
    scatter(robot_positions(:,1), robot_positions(:,2), 6, 'filled',...
        'MarkerFaceColor', [0 0.4470 0.7410]); % 蓝色轨迹点
    scatter(sensor_positions(:,1), sensor_positions(:,2), 6, 'filled',...
        'MarkerFaceColor', [0.3010 0.7450 0.9330]); % 青色传感器轨迹
    
    % 绘制边线标记（品红色星号）
    if ~isempty(sidelines)
        scatter(sideline_positions(:,1), sideline_positions(:,2), 100,...
            'm*', 'LineWidth', 1.5);
        text(sideline_positions(:,1)+50, sideline_positions(:,2)+50,...
            string(1:length(sidelines)), 'Color', 'm');
    end
    
    % 绘制跨线标记（红色叉号）
    if ~isempty(crosslines)
        scatter(crossline_positions(:,1), crossline_positions(:,2), 100,...
            'rx', 'LineWidth', 1.5);
        text(crossline_positions(:,1)+50, crossline_positions(:,2)+50,...
            string(1:length(crosslines)), 'Color', 'r');
    end
    
    hold off
    axis equal  % 保持坐标轴比例一致
    grid on     % 显示网格
    xlabel('X位置 (mm)')
    ylabel('Y位置 (mm)')
end

function plotRadius(distances, thetas)
% 绘制转弯半径变化曲线
% 输入：平移距离数组，旋转角度数组
% 输出：半径变化曲线图
    
    % 防止除零错误
    thetas(thetas == 0) = 1e-5;  
    
    % 半径计算（使用弦长公式）
    radius = abs(distances ./ thetas);  
    radius(radius > 2000) = 2000;  % 设置显示上限
    
    % 可视化设置
    plot(radius, 'LineWidth', 1.5)
    ylim([0 2100])
    grid on
    xlabel('时间步')
    ylabel('转弯半径 (mm)')
end

function plotdTheta(distances, thetas)
% 绘制归一化角速度变化曲线
% 输入：平移距离数组，旋转角度数组
% 输出：角速度变化曲线图
    
    % 防止除零错误
    thetas(thetas == 0) = 1e-5;  
    
    % 计算角速度（弧度/mm）
    angular_velocity = abs(thetas ./ distances);  
    
    % 可视化
    plot(angular_velocity, 'LineWidth', 1.5)
    grid on
    xlabel('时间步')
    ylabel('角速度 (rad/mm)')
end

function plotVelocityTable(distances, thetas, Vmax, Vmin, R_straight, func_type)
% 速度规划表生成器
% 输入参数：
%   distances  - 各段平移距离
%   thetas     - 各段旋转角度
%   Vmax       - 直线段最大速度 (m/s)
%   Vmin       - 弯道最小速度 (m/s)
%   R_straight - 直线判定阈值 (mm)
%   func_type  - 速度插值函数类型 ('linear','quadratic','sigmoid')
    
    % 计算瞬时转弯半径
    thetas(thetas == 0) = 1e-5;
    radius = abs(distances ./ thetas);
    radius(radius > R_straight) = R_straight;
    
    % 根据函数类型生成速度规划
    switch func_type
        case 'linear'
            velocity = (Vmax - Vmin)/R_straight * radius + Vmin;
        case 'quadratic'
            velocity = 1e-3*(Vmax - Vmin)/R_straight * radius.^2 + Vmin;
        case 'sigmoid'
            k = 12/R_straight;  % 斜率调节因子
            velocity = (Vmax - Vmin)./(1 + exp(-k*(radius - R_straight/2))) + Vmin;
        otherwise
            error('不支持的函数类型');
    end
    
    % 可视化
    plot(velocity, 'LineWidth', 1.5)
    ylim([Vmin-0.5 Vmax+0.5])
    grid on
    xlabel('时间步')
    ylabel('规划速度 (m/s)')
end