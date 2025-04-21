clear;  % 清空工作区变量
clc;    % 清空命令窗口

% ================= 参数配置 =================
straight = 1000;        % 直线基准距离 (单位: mm)
adjust_x = 25;          % Sigmoid水平缩放因子 (控制曲线起始位置)
adjust_y = 0;           % 垂直偏移量 (调整曲线起始点)
gain = 0.45;            % 曲线陡峭度系数 (值越大曲线越陡)

% ================= 计算生成 =================
factor = straight/adjust_x;             % 实际缩放因子
radius = 0:0.1:straight;                % 生成半径采样点 (0-1000mm,步长0.1mm)
v_max = 6.5;                            % 最大速度 (m/s)
v_min = 2.5;                            % 最小速度 (m/s)

% Sigmoid速度曲线计算
velo_sigmoid = (1 ./ (1 + exp(-(gain/factor)*radius + (adjust_x/2)*gain)))... 
               * (v_max - v_min) + v_min + adjust_y;

% 对比曲线生成
linear_velo = linspace(v_min, v_max, length(radius));           % 线性速度曲线
quadratic_velo = 1e-3 .* radius.^2 * ((v_max - v_min)/straight) + v_min; % 二次函数速度曲线

% ================= 可视化 =================
figure('Name','速度规划策略对比','NumberTitle','off')  % 创建命名图形窗口

% 绘制主要曲线
plot(radius, velo_sigmoid, 'LineWidth', 2)        % Sigmoid曲线
hold on
plot(radius, linear_velo, '--', 'LineWidth', 1.5) % 线性曲线
plot(radius, quadratic_velo, ':', 'LineWidth', 2)  % 二次曲线
hold off

% 图形美化设置
grid on
xlabel('转弯半径 (mm)', 'FontSize', 12)
ylabel('规划速度 (m/s)', 'FontSize', 12)
title('\bf 不同速度规划策略对比', 'FontSize', 14)  % 使用加粗字体
legend('Sigmoid函数', '线性函数', '二次函数',...
       'Location', 'southeast')                  % 图例放置在东南方
set(gca, 'FontSize', 11)                         % 设置坐标轴字体大小
xlim([0 straight])                               % X轴范围固定

% ================= 关键点标记 =================
hold on
% 标记Sigmoid曲线特征点
plot(0, v_min, 'ro', 'MarkerSize', 8, 'LineWidth', 2)       % 起始点
plot(straight, v_max, 'rh', 'MarkerSize', 8, 'LineWidth', 2) % 终点
text(50, v_min+0.2, '起始速度', 'Color', 'r', 'FontSize', 10)
text(straight-150, v_max-0.3, '极限速度', 'Color', 'r', 'FontSize', 10)
hold off