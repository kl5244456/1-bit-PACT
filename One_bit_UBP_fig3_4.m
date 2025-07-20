% clear
% clc
tic
load('Arrangement1k20241018_after_calibration.mat');
frame_num = 100;
steps = frame_num;
theta = (0:steps - 1)'  / steps * pi / 4;
% --- after calibration ----
N_transducer = sum(trans_select);
transducer_x_r = zeros(steps * N_transducer, 1);
transducer_y_r = zeros(steps * N_transducer, 1);
transducer_z_r = zeros(steps * N_transducer, 1);
for k = 1:steps
    [xr, yr, zr] = rotationAboutAxis20230830(...
        transducer_x_ms, transducer_y_ms, transducer_z_ms,...
        PointOnAx, Ax, theta(k));
    transducer_x_r((k - 1) * N_transducer + 1:k * N_transducer) = xr;
    transducer_y_r((k - 1) * N_transducer + 1:k * N_transducer) = yr;
    transducer_z_r((k - 1) * N_transducer + 1:k * N_transducer) = zr;
end

transducer.x = transducer_x_r;
transducer.y = transducer_y_r;
transducer.z = transducer_z_r;

R = sqrt(transducer_x_r.^2 + transducer_y_r.^2 + transducer_z_r.^2);
transducer.x_norm = - transducer_x_r ./ R;
transducer.y_norm = - transducer_y_r ./ R;
transducer.z_norm = - transducer_z_r ./ R;
ds = sqrt(R.^2 - transducer.z.^2);
%%
Path = '/scratch/haitang/Binary paper data/new hand data/';
filename= 'hand_100_30dB_h';
load([Path, filename]);
%%
Sinogram = VOLTAGE(id_rearrange, :, 1:frame_num);
Sinogram = double(permute(Sinogram(trans_select, :, :), [2, 1, 3]));
%%
Sinogram = reshape(Sinogram, size(Sinogram, 1), []);
Sinogram = Sinogram.*ds';
%% add noise to the sinogram (0 for weak features' reconstruction; )
Sinogram1 = Sinogram;
noise_level = 0.1*max(Sinogram1(:)); % Noise standard deviation
noise = noise_level * randn(size(Sinogram1)); % Artificial gaussian noise
Sinogram = Sinogram1 + noise; % Noisy signal
%% binary
minV = min(Sinogram(:));
maxV = max(Sinogram(:));
Sinogram(Sinogram>0) = maxV * sign(Sinogram(Sinogram>0));
Sinogram(Sinogram<0) = abs(minV) * sign(Sinogram(Sinogram<0));
%%
csys.px = 2e-4; csys.x0 =-0.08; csys.x1 =  0.08-csys.px; 
csys.py = 2e-4; csys.y0 =-0.08; csys.y1 =  0.08-csys.py; 
csys.pz = 2e-4; csys.z0 =-0.00; csys.z1 =  0.025-csys.pz; 

csys = formCsys20221108(csys);

F_sample = 20e6;
acq_delay = -1/F_sample;  
Device      = 1;
%% temporal low-pass filters
cutoffFreq = 5e6; % Cutoff frequency 4 MHz
% Design Butterworth low-pass filter
filterOrder = 4; % Filter order, can be adjusted as needed
[b, a] = butter(filterOrder, cutoffFreq / (F_sample / 2), 'low'); % Low-pass filter
Sinogram = filtfilt(b, a, Sinogram);
%% temporal high-pass filters
cutoffFreq = 0.5e6; % Cutoff frequency 4 MHz
% Design Butterworth high-pass filter
filterOrder = 4; % Filter order, can be adjusted as needed
[b, a] = butter(filterOrder, cutoffFreq / (F_sample / 2), 'high'); % Low-pass filter
Sinogram = filtfilt(b, a, Sinogram);
%%
file_name = 'Code/CudaNote.txt';
V_sound   = 1530;
Batch = 3072 * 4 * 4 * 4 *4*4;
%%
message = 'UBP';
Mode = 2;

image = cudaTmDomnRecon3D20240618(32, 32, 32, 32,...
                                         Sinogram, Batch,...
                                         transducer, csys,...
                                         V_sound, F_sample, acq_delay,...
                                         Mode, file_name, Device, message);

toc
%% imaging
temp = double(image);   
temp(temp<0) = 0;
temp = abs(temp);
temp = temp/max(temp(:));
figure;
all_slices = cat3d(temp);
imshow(all_slices, []); colormap hot;colorbar;caxis([0 1])