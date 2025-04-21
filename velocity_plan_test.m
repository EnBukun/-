clear 

% 定义基础参数
r_max = 1000;            % 最大转弯半径(mm)
radius = 0:0.1:r_max;    % 生成半径采样点(0-1000mm,步长0.1)
v_max = 7;               % 直线段最大速度(m/s)
v_min = 2.5;             % 最小速度(m/s)

% 分段函数参数
r1 = 400;                % 第一临界半径
r2 = 800;                % 第二临界半径
ratio1 = 0.4;            % 第一区间调整系数
ratio2 = 0.9;            % 第二区间调整系数

% 计算基础线性函数参数
a = (v_max - v_min) / r_max; % 整体线性斜率

% 第一区间参数计算
v1 = linearFunc(a, r1, 0, v_min);    % 基础线性函数在r1处的值
vv1 = (v1-v_min)*ratio1 + v_min;     % 调整后的第一区间终点速度
a1 = (vv1 - v_min) / r1;             % 第一区间新斜率

% 第二区间参数计算
v2 = linearFunc(a, r2, 0, v_min);    % 基础线性函数在r2处的值 
vv2 = (v2-v_min)*ratio2 + v_min;     % 调整后的第二区间终点速度
a2 = (vv2 - vv1)/(r2 - r1);          % 第二区间新斜率

% 第三区间参数计算
v3 = v_max;              % 最终速度设为最大值
vv3 = v_max;             
a3 = (vv3 - vv2)/(r_max - r2); % 第三区间斜率

% 初始化速度曲线存储数组
velo = zeros(1, length(radius));      % 条件分段线性函数
velo1 = zeros(1, length(radius));     % 基础线性函数
velo2 = zeros(1, length(radius));     % 二次函数

% 速度曲线计算
for i = 1 : length(radius)
    % 条件分段线性函数
    if radius(i) < r1
        velo(i) = linearFunc(a1, radius(i), 0, v_min); 
    elseif radius(i) < r2
        velo(i) = linearFunc(a2, radius(i), r1, vv1);
    else
        velo(i) = linearFunc(a3, radius(i), r2, vv2);
    end
    
    % 基础线性函数
    velo1(i) = linearFunc(a, radius(i), 0, v_min);
    
    % 二次函数
    velo2(i) = 1e-3 * radius(i)^2 * ((v_max - v_min)/r_max) + v_min;
end

% Sigmoid函数参数调整
adjust_x = 25;       % 水平方向位置调整系数
adjust_y = 0;        % 垂直方向偏移量
gain = 0.30;         % 曲线陡峭度控制参数
factor = r_max/adjust_x;

% Sigmoid函数速度曲线
velo_sigmoid = (1 ./ (1 + exp(-(gain/factor)*radius + (adjust_x/2)*gain)))...
              * (v_max - v_min) + v_min + adjust_y;

% 多项式拟合
n = 4;  % 多项式阶数
x0 = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]; % 采样点X
y0 = [2.5, 2.5, 2.8, 3.0, 3.3, 4.3, 4.5, 5.0, 5.5, 6.3, 7];   % 采样点Y
p = polyfit(x0, y0, n);                % 多项式系数拟合
velo_poly = polyval(p, radius);         % 多项式曲线生成

% 可视化绘图
figure
plot(radius, velo, 'LineWidth', 2)      % 条件分段线性函数
hold on
plot(radius, velo1, '--', 'LineWidth', 1.5) % 基础线性函数
plot(radius, velo2, ':', 'LineWidth', 2)    % 二次函数
plot(radius, velo_sigmoid, '-.', 'LineWidth', 2) % Sigmoid函数
plot(radius, velo_poly, 'LineWidth', 2)      % 多项式曲线
scatter(x0, y0, 100, 'filled')          % 原始采样点显示

% 图例设置
legend('条件分段线性', '基础线性', '二次函数', 'Sigmoid调整', '四次多项式',...
       '采样点', 'Location','northwest')
xlabel('转弯半径 (mm)')
ylabel('规划速度 (m/s)')
grid on
xlim([0 1000])
set(gca, 'FontSize', 12)

% 线性计算函数
function v = linearFunc(a, r, r_shift, b)
    % 参数说明:
    % a - 斜率
    % r - 当前半径
    % r_shift - 截距偏移量
    % b - 基础截距
    v = a * (r - r_shift) + b;
end