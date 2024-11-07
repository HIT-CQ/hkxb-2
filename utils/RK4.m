function [t, state] = RK4(odefun, tspan, x0, h)
    t = tspan(1):h:tspan(2); % 时间向量
    n = length(t); % 步数
    state = zeros(length(x0), n); % 初始化状态矩阵
    state(:, 1) = x0; % 设置初始条件

    for i = 1:n-1
        k1 = h * odefun(t(i), state(:, i));
        k2 = h * odefun(t(i) + h/2, state(:, i) + k1/2);
        k3 = h * odefun(t(i) + h/2, state(:, i) + k2/2);
        k4 = h * odefun(t(i) + h, state(:, i) + k3);
        
        state(:, i+1) = state(:, i) + (k1 + 2*k2 + 2*k3 + k4) / 6; % 更新状态
    end
    t = t.';
    state = state.';
end