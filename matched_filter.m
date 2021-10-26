% 可调参数：过采样倍数oversamp、高斯滤波器BT、相关长度L

clear all;
close all;

oversamp = 8; % 过采样倍数，可调
bit_rate = 16e6;  % 符号速率（这个值不重要）
Tb = 1/bit_rate;  % 符号时间
BbTb = 0.3;       % BT参数，可调
Bb = BbTb / Tb;
fs = bit_rate * oversamp; % 采样率
dt = 1/fs;        % 采样间隔
L = 5;            % 相关长度L，可调

% 生成高斯脉冲响应波形g(t)
t = -L/2 * Tb : 1/fs: L/2 * Tb - 1/fs;
t = t + 1/fs/2;
g = erfc(2 * pi * Bb * (t - Tb / 2) / sqrt(log(2)) / sqrt(2)) / 2 - erfc(2 * pi * Bb * (t + Tb / 2) / sqrt(log(2)) / sqrt(2)) / 2;

% 做积分
for i = 1:length(g)

    if i == 1
        q_g(i) = g(i) * dt;
    else
        q_g(i) = q_g(i - 1) + g(i) * dt;
    end
end
q_g = q_g / 2 / Tb;

% figure;
% plot(q_g);

s0(1:length(q_g)) = sin(pi*q_g);
s0(length(q_g) + 1 : 2 * length(q_g)) = cos(pi*q_g); % 共2*L*Tb范围

% 计算c0和c1
t0 = 1 : (L + 1) * oversamp; % 0~(L+1)*Tb范围内
c0 = s0(t0);
for i = 1 : L - 1
    c0 = c0 .* s0(t0 + i * oversamp);
end

t1 = 1 : (L - 1) * oversamp; % 0~(L-1)*Tb范围内
c1 = s0(t1);
for i = 1 : L - 1
    if i == 1
        c1 = c1 .* s0(t1 + (i + L) * oversamp);
    else
        c1 = c1 .* s0(t1 + i * oversamp);
    end
end

figure;
plot(c0);
hold on;
plot(c1);
