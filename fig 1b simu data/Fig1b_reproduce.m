load('one_bit_simu');
F_sample = 1000e6;
x=(0:1/F_sample:1e-6);
figure;plot(x*1e6,f)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18)
title('Raw signal')
xlabel('t (ms)');
ylabel('Amplitude');

figure;plot(x*1e6,ff)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18)
title('Noisy signal')
xlabel('t (ms)');
ylabel('Amplitude');

ffb = sign(ff);
figure;plot(x*1e6,ffb)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18)
title('1-bit signal')
xlabel('t (ms)');
%% temporal low-pass filters
cutoffFreq = 12e6; % initial 12 MHz
% Design Butterworth low-pass filter
filterOrder = 4; % Filter order, can be adjusted as needed
[b, a] = butter(filterOrder, cutoffFreq / (F_sample / 2), 'low'); % Low-pass filter
ffb = filtfilt(b, a, ffb);
%% temporal high-pass filters
cutoffFreq = 3e6; % initial 3 MHz
% Design Butterworth high-pass filter
filterOrder = 4; % Filter order, can be adjusted as needed
[b, a] = butter(filterOrder, cutoffFreq / (F_sample / 2), 'high'); % Low-pass filter
ffb = filtfilt(b, a, ffb);
figure;plot(x*1e6,ffb)
hold on;plot(x*1e6,f)
set(gca, 'FontSize', 18); 
set(gca, 'FontName', 'Times New Roman')
xlabel('t (ms)');
ylabel('Amplitude');
legend('1-bit restored','Raw')