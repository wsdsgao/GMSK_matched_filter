% �ɵ�����������������oversamp����˹�˲���BT����س���L

clear all;
close all;

oversamp = 8; % �������������ɵ�
bit_rate = 16e6;  % �������ʣ����ֵ����Ҫ��
Tb = 1/bit_rate;  % ����ʱ��
BbTb = 0.3;       % BT�������ɵ�
Bb = BbTb / Tb;
fs = bit_rate * oversamp; % ������
dt = 1/fs;        % �������
L = 5;            % ��س���L���ɵ�

% ���ɸ�˹������Ӧ����g(t)
t = -L/2 * Tb : 1/fs: L/2 * Tb - 1/fs;
t = t + 1/fs/2;
g = erfc(2 * pi * Bb * (t - Tb / 2) / sqrt(log(2)) / sqrt(2)) / 2 - erfc(2 * pi * Bb * (t + Tb / 2) / sqrt(log(2)) / sqrt(2)) / 2;

% ������
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
s0(length(q_g) + 1 : 2 * length(q_g)) = cos(pi*q_g); % ��2*L*Tb��Χ

% ����c0��c1
t0 = 1 : (L + 1) * oversamp; % 0~(L+1)*Tb��Χ��
c0 = s0(t0);
for i = 1 : L - 1
    c0 = c0 .* s0(t0 + i * oversamp);
end

t1 = 1 : (L - 1) * oversamp; % 0~(L-1)*Tb��Χ��
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
