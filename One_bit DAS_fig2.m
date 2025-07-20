clearvars;

% =========================================================================
% SIMULATION
% =========================================================================

% assign the grid size and create the computational grid
PML_size = 20;              % size of the PML in grid points
Nx = 512 - 2 * PML_size;    % number of grid points in the x direction
Ny = 512 - 2 * PML_size;    % number of grid points in the y direction
x = 10e-2;                  % total grid size [m]
y = 10e-2;                  % total grid size [m]
dx = x / Nx;                % grid point spacing in the x direction [m]
dy = y / Ny;                % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

disc_magnitude = 50; % [Pa]
disc_x_pos = 80*2;    % [grid points]
disc_y_pos = 110*2;    % [grid points]
disc_radius = 4;    % [grid points]
disc_3 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

disc_magnitude = 1; % [Pa]
disc_x_pos = 110*2;    % [grid points]
% disc_y_pos = 119*2;    % [grid points]
disc_y_pos = 115*2;    % [grid points]
disc_radius = 4;    % [grid points]
disc_2 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

source.p0 = disc_2 + disc_3;
p0 = source.p0;

% define the properties of the propagation medium
medium.sound_speed = 1560;  % [m/s]
medium.density = 1000;                    % [kg/m^3]
% medium.alpha_coeff = 2.4;      % [dB/(MHz^y cm)]
% medium.alpha_power = 2.9;
% medium.BonA = 6;

% define a centered Cartesian circular sensor
sensor_radius = 4.5e-2;     % [m]
sensor_angle = pi;      % [rad]
sensor_pos = [0, 0];        % [m]
num_sensor_points = 283*2;
cart_sensor_mask = makeCartCircle(sensor_radius, num_sensor_points, sensor_pos, sensor_angle);

% assign to sensor structure
sensor.mask = cart_sensor_mask;

% create the time array
kgrid.makeTime(medium.sound_speed,0.1);
t_end = 54e-6;
Nt = round(t_end/kgrid.dt);
kgrid.setTime(Nt,kgrid.dt);

% set the input options
input_args = {'Smooth', false, 'PMLInside', false, 'PlotPML', false};

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% add noise to the recorded sensor data
signal_to_noise_ratio = 40;	% [dB]
sensor_data = addNoise(sensor_data, signal_to_noise_ratio, 'peak');

% filter the sensor data using a Gaussian filter
Fs = 1/kgrid.dt;        % [Hz]
center_freq = 2e6;      % [Hz]
bandwidth = 100;        % [%]
sensor_data = gaussianFilter(sensor_data, Fs, center_freq, bandwidth);
%% VISUALISATION
% plot the initial pressure and sensor distribution
figure;
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, p0 + cart2grid(kgrid, cart_sensor_mask), [-1, 1]);
colormap(getColorMap);
ylabel('x/mm');
xlabel('y/mm');
axis image;
xlim([-50 50]);ylim([-50 50])
xticks(-50:25:50)
yticks(-50:25:50)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
title('Sensor array & Sources')

%% 1-bit DAS
kgrid.dt = 1.3581e-08;
sensor_data1 = sign(sensor_data);
ss=zeros(size(p0));
for ii=1:472
    for jj=1:472
        for kk=1:size(sensor_data,1)%number of sensors
            if ceil(sqrt(((ii-236)*dx-cart_sensor_mask(1,kk)).^2+((jj-236)*dx-cart_sensor_mask(2,kk)).^2)/medium.sound_speed/kgrid.dt) == 0
                ss(ii,jj)=ss(ii,jj) + 0;
            elseif ceil(sqrt(((ii-236)*dx-cart_sensor_mask(1,kk)).^2+((jj-236)*dx-cart_sensor_mask(2,kk)).^2)/medium.sound_speed/kgrid.dt) > size(sensor_data1,2)
                ss(ii,jj)=ss(ii,jj) + 0;
            else
        ss(ii,jj)=ss(ii,jj)+sensor_data1(kk, ceil(sqrt(((ii-236)*dx-cart_sensor_mask(1,kk)).^2+((jj-236)*dx-cart_sensor_mask(2,kk)).^2)/medium.sound_speed/kgrid.dt));
            end
            end
    end
end

% Nx = 512 - 2 * PML_size;    % number of grid points in the x direction
% Ny = 512 - 2 * PML_size;    % number of grid points in the y direction
% x = 10e-2;                  % total grid size [m]
% y = 10e-2;                  % total grid size [m]
% dx = x / Nx;                % grid point spacing in the x direction [m]
% dy = y / Ny;                % grid point spacing in the y direction [m]
% kgrid = kWaveGrid(Nx, dx, Ny, dy);
ss(ss<0) = 1;
figure;
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, abs(ss/max(abs(ss(:)))));
colormap(getColorMap);
ylabel('x/mm');
xlabel('y/mm');
axis image;
colorbar
xlim([-50 50]);ylim([-50 50])
xticks(-50:25:50)
yticks(-50:25:50)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18)
title('1-bit Recon (linear)')
caxis([0 1])

figure;
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, 10*log10(abs(ss/max(abs(ss(:))))));
colormap(getColorMap);
ylabel('x/mm');
xlabel('y/mm');
axis image;
xlim([-50 50]);ylim([-50 50])
xticks(-50:25:50)
yticks(-50:25:50)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18)
colorbar
title('1-bit Recon (log)')
caxis([-12 0])
%% conventional DAS
kgrid.dt = 1.3581e-08;
ss=zeros(size(p0));
for ii=1:472
    for jj=1:472
        for kk=1:size(sensor_data,1)%number of sensors
            if ceil(sqrt(((ii-236)*dx-cart_sensor_mask(1,kk)).^2+((jj-236)*dx-cart_sensor_mask(2,kk)).^2)/medium.sound_speed/kgrid.dt) == 0
                ss(ii,jj)=ss(ii,jj) + 0;
            elseif ceil(sqrt(((ii-236)*dx-cart_sensor_mask(1,kk)).^2+((jj-236)*dx-cart_sensor_mask(2,kk)).^2)/medium.sound_speed/kgrid.dt) > size(sensor_data1,2)
                ss(ii,jj)=ss(ii,jj) + 0;
            else
        ss(ii,jj)=ss(ii,jj)+sensor_data(kk, ceil(sqrt(((ii-236)*dx-cart_sensor_mask(1,kk)).^2+((jj-236)*dx-cart_sensor_mask(2,kk)).^2)/medium.sound_speed/kgrid.dt));
            end
            end
    end
end

% Nx = 512 - 2 * PML_size;    % number of grid points in the x direction
% Ny = 512 - 2 * PML_size;    % number of grid points in the y direction
% x = 10e-2;                  % total grid size [m]
% y = 10e-2;                  % total grid size [m]
% dx = x / Nx;                % grid point spacing in the x direction [m]
% dy = y / Ny;                % grid point spacing in the y direction [m]
% kgrid = kWaveGrid(Nx, dx, Ny, dy);
ss(ss<0) = 1;
figure;
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, abs(ss/max(abs(ss(:)))));
colormap(getColorMap);
ylabel('x/mm');
xlabel('y/mm');
axis image;
colorbar
xlim([-50 50]);ylim([-50 50])
xticks(-50:25:50)
yticks(-50:25:50)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18)
title('Conventional Recon (linear)')
caxis([0 1])

figure;
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, 10*log10(abs(ss/max(abs(ss(:))))));
colormap(getColorMap);
ylabel('x/mm');
xlabel('y/mm');
axis image;
xlim([-50 50]);ylim([-50 50])
xticks(-50:25:50)
yticks(-50:25:50)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18)
colorbar
title('Conventional Recon (log)')
caxis([-30 0])