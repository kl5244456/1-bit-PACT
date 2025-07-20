%% for weak signal
% signal produce
F_sample = 1000e6;
x=(0:1/F_sample:1e-6);
f = zeros(size(x));
for ii=5:1:10
    f = f+(-1)^ii*sin(2*pi*ii*1e6*x);
end
f = f/max(abs(f));
% figure;plot(x,f)

v = VideoWriter('compare.mp4', 'MPEG-4');
v.FrameRate = 10;
open(v);

blueColor = [0, 114, 189]/255;
orangeColor = [217, 83, 25]/255;

g = zeros(size(f));
h = zeros(size(f));
figure('Color', 'w', 'Position', [100, 100, 1200, 600]); 

tiled = tiledlayout(1, 2);
tiled.TileSpacing = 'compact';
tiled.Padding = 'compact';
for ii=1:100
% noise; 2 for strong noise simulation
mag = 2;
noise = mag*randn(size(f));
ff = f+noise;
% figure;plot(x,ff)
ffb = sign(ff);
ffb1 = ffb;
% figure;plot(x,ffb)

% temporal low-pass filters
cutoffFreq = 12e6; % Cutoff frequency 4 MHz
% Design Butterworth low-pass filter
filterOrder = 4; % Filter order, can be adjusted as needed
[b, a] = butter(filterOrder, cutoffFreq / (F_sample / 2), 'low'); % Low-pass filter
ffb = filtfilt(b, a, ffb);

% temporal high-pass filters
cutoffFreq = 3e6; % Cutoff frequency 4 MHz
% Design Butterworth high-pass filter
filterOrder = 4; % Filter order, can be adjusted as needed
[b, a] = butter(filterOrder, cutoffFreq / (F_sample / 2), 'high'); % Low-pass filter
ffb = filtfilt(b, a, ffb);
g = g + ffb;
g1 = g/max(abs(g));
h = h + ff;
h1 = h/max(abs(h));
fontSize = 18;
 
    nexttile(1);
plot(x*1e6,h1, 'Color', blueColor, 'LineWidth', 2);
hold on;plot(x*1e6,f, 'Color', orangeColor, 'LineWidth', 2);
title(sprintf('Direct averaging: #%d', ii), ...
        'FontName', 'Times New Roman', 'FontSize', 18);
hold off;
xlim([0 1]);
ylim([-1 1]);
set(gca, 'FontSize', 18); 
set(gca, 'FontName', 'Times New Roman')
xlabel('t (ms)');
ylabel('Amplitude');
legend('Ensemble average','Raw', 'Location', 'northwest')
legend('boxoff')

    nexttile(2);
plot(x*1e6,g1, 'Color', blueColor, 'LineWidth', 2);
hold on;plot(x*1e6,f, 'Color', orangeColor, 'LineWidth', 2);
title(sprintf('1-bit: #%d', ii), ...
        'FontName', 'Times New Roman', 'FontSize', 18);
hold off;
xlim([0 1]);
ylim([-1 1]);
set(gca, 'FontSize', 18); 
set(gca, 'FontName', 'Times New Roman')
xlabel('t (ms)');
ylabel('Amplitude');
legend('Ensemble average','Raw', 'Location', 'northwest')
legend('boxoff')
drawnow;
% write
    frame = getframe(gcf);
    writeVideo(v, frame);
    pause(0.1);
end
close(v);
disp('video saved');