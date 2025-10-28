%% ========================================================================
%  Figure 2 & 3: Calculating exceptional point surfaces
%  ------------------------------------------------------------------------
%  This script computes and visualizes the exceptional point surface
%  forming the elliptic umbilic neighborhood topology of a pair of linear
%  exceptional points in nonlinear parameter space.
%  ------------------------------------------------------------------------
%  The second part of the script uses the data to create a video of the
%  umbilic "growing" out of the two exceptional points lying in the plane
%  forming the linear parameterspace
%
%  Author: [Jan Wingenbach]
% ========================================================================
clear; close all; clc;

%% ------------------------ USER PARAMETERS -------------------------------
get_video = 0;                  %create video of umbilic
x0 = 10; y0 = 10;
width = 1200; height = 1000;
az = 140; el = 40; lw = 3;
color = [0, 0.7, 0.7];

linecolors = [
    0.0011    0.0001    0.0116
    0.4703    0.1098    0.4286
    0.9307    0.4135    0.1437
    0.9913    1.0000    0.6517
];

% download "inferno.mat" from here:
% https://github.com/JanWI1/nonlinear-exceptional-points-topology/code
% otherwise set cm to 'jet'
load('inferno.mat','inferno');
cm = inferno;

%% ------------------ UMBILIC PARAMETER DEFINITIONS -----------------------
alpha_step = 5e-2;
alpha_min = 0;
alpha_max_index = 201;
num_teta = 200;

% the teta-offset ensures avoidance of divergence
offset = 1e-12;

% Define four distinct linspace intervals for teta
teta_intervals = {
    @(alpha_gamma) linspace(3 * pi / 2 + offset, atan(-1 / alpha_gamma) + 2 * pi - offset, num_teta),
    @(alpha_gamma) linspace(atan(-1 / alpha_gamma) + pi + offset , 3 * pi / 2 - offset, num_teta),
    @(alpha_gamma) linspace(atan(-1 / alpha_gamma) + offset,  pi/2 - offset, num_teta),
    @(alpha_gamma) linspace( pi/ 2 + offset, atan(-1 / alpha_gamma) + pi - offset, num_teta)
    };

% Pre-allocate arrays
alpha_gamma_ar = zeros(1, alpha_max_index);
surface_data = cell(length(teta_intervals), 2); % Store data for each teta interval and its surfaces
beta_gamma_ar = zeros(alpha_max_index, num_teta);       % Normalized
delta_m_gamma_ar = zeros(alpha_max_index, num_teta);    % Normalized

%% -------------------- COMPUTE UMBILIC CROSS-SECTION ---------------------
for teta_idx = 1:length(teta_intervals)
    beta_gamma_unnorm = zeros(alpha_max_index, num_teta);
    delta_m_gamma_unnorm = zeros(alpha_max_index, num_teta);
    delta_p_gamma_unnorm = zeros(alpha_max_index, num_teta);

    for i = 1:alpha_max_index
        if teta_idx == 2
            alpha_gamma = -i * alpha_step + alpha_min;
            teta = teta_intervals{teta_idx}(alpha_gamma);
            beta_gamma = sqrt(1 ./ (sin(teta).^3 .* (sin(teta) - alpha_gamma .* cos(teta))));
        elseif teta_idx == 1
            alpha_gamma = i * alpha_step + alpha_min;
            teta = teta_intervals{teta_idx}(alpha_gamma);
            beta_gamma = sqrt(1 ./ (sin(teta).^3 .* (sin(teta) - alpha_gamma .* cos(teta))));
            s = 2 * (sqrt(1 + alpha_gamma^2) - 1);
            beta_gamma_ar(i, :) = beta_gamma/ s;
        elseif teta_idx == 3
            alpha_gamma = -i * alpha_step + alpha_min;
            teta = teta_intervals{teta_idx}(alpha_gamma);
            beta_gamma = -sqrt(1 ./ (sin(teta).^3 .* (sin(teta) - alpha_gamma .* cos(teta))));
        elseif teta_idx == 4
            alpha_gamma = i * alpha_step + alpha_min;
            teta = teta_intervals{teta_idx}(alpha_gamma);
            beta_gamma = -sqrt(1 ./ (sin(teta).^3 .* (sin(teta) - alpha_gamma .* cos(teta))));
        end

        alpha_gamma_ar(i) = alpha_gamma;
        p = beta_gamma .* sin(teta);

        % Calculate delta_m_gamma
        delta_m_gamma = (1 ./ p - p ./ (1 - sqrt(1 - p.^2))) .* (alpha_gamma .* p + beta_gamma .* cos(teta));
        s = 2 * (sqrt(1 + alpha_gamma^2) - 1);
        delta_m_gamma_ar(i, :) = delta_m_gamma / s;

        % Store unnormalized values
        beta_gamma_unnorm(i, :) = beta_gamma;
        delta_m_gamma_unnorm(i, :) = delta_m_gamma;
        delta_p_gamma_unnorm(i, :) = -delta_m_gamma;
    end

    % Store surface data for plotting
    surface_data{teta_idx, 1} = beta_gamma_unnorm;
    surface_data{teta_idx, 2} = delta_m_gamma_unnorm;
end
delta_p_gamma_ar = -delta_m_gamma_ar;

%% ---------------------- PLOTTING SECTION -------------------------------

figure(1); clf; hold on;
colormap(cm);

%% Umbilic surface in δ–β-α space
for teta_idx = 1:length(teta_intervals)
    beta_gamma_unnorm = surface_data{teta_idx, 1};
    delta_m_gamma_unnorm = surface_data{teta_idx, 2};

    delta_p_gamma_unnorm = -delta_m_gamma_unnorm;
    if teta_idx == 2 || teta_idx ==3
        surf(beta_gamma_unnorm, -alpha_gamma_ar, delta_m_gamma_unnorm, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
        surf(beta_gamma_unnorm, -alpha_gamma_ar, delta_p_gamma_unnorm, 'FaceAlpha', 0.8, 'EdgeColor', 'none');

        edgelines(beta_gamma_unnorm,delta_p_gamma_unnorm,delta_m_gamma_unnorm,alpha_max_index,lw,'minus',-alpha_gamma_ar)
    elseif teta_idx == 1 || teta_idx ==4
        surf(beta_gamma_unnorm, alpha_gamma_ar, delta_m_gamma_unnorm, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
        surf(beta_gamma_unnorm, alpha_gamma_ar, delta_p_gamma_unnorm, 'FaceAlpha', 0.8, 'EdgeColor', 'none');

        edgelines(beta_gamma_unnorm,delta_p_gamma_unnorm,delta_m_gamma_unnorm,alpha_max_index,lw,'plus',alpha_gamma_ar);
    end
end

hold off;box on;set(gca, 'linewidth', lw - 1);shading interp;
axis tight;grid on;view(az, el);set(gca, 'fontsize', 18);

% Define axis limits
x_min = -10.5; x_max = 10.5;
y_min = -10;   y_max = 10;
z_min = -8.5;  z_max = 8.5;

x_mid = (x_min + x_max) / 2;y_mid = (y_min + y_max) / 2;z_mid = (z_min + z_max) / 2;
xlim([x_min x_max]);ylim([y_min y_max]);zlim([z_min z_max]);

% Set ticks at ends and center
xticks([x_min x_mid x_max]);
yticks([y_min y_mid y_max]);
zticks([z_min z_mid z_max]);
set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);set(gca,'ZTickLabel',[]);
set(gcf, 'position', [x0, y0, width, height]);

%% Umbilic cross-sections line in δ–β space
w_values = [0.1,1,10];
figure(2);clf;hold on;
line([-10 10], [0 0], 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', lw); % Horizontal line
line([0 0], [-10 10], 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', lw); % Vertical line
line([0.5 10], [0 0], 'Color', color, 'LineStyle', '-', 'LineWidth', 3*lw); % Vertical line
for j = 1:length(w_values)
    plot(beta_gamma_ar(alpha_gamma_ar==w_values(j),:)-beta_gamma_ar(alpha_gamma_ar==w_values(j),1),delta_m_gamma_ar(alpha_gamma_ar==w_values(j),:),'color',linecolors(j,:), 'LineWidth', lw*2);
    plot(beta_gamma_ar(alpha_gamma_ar==w_values(j),:)-beta_gamma_ar(alpha_gamma_ar==w_values(j),1),delta_p_gamma_ar(alpha_gamma_ar==w_values(j),:),'color',linecolors(j,:), 'LineWidth', lw*2);
end
hold off;ylim([-0.5 0.5]);xlim([-0.2 0.6]);box on;set(gca,'linewidth',lw-1);
set(gcf,'position',[x0,y0,width,height]);xlabel("beta'/(\gamma s)");ylabel('delta/(\gamma s)');set(gca,'fontsize',18)

%% Video visualizing "growth" of the umbilic
if get_video == 1
    profile on;
    % Video settings
    video_filename = 'umbilic_growth.mp4';
    v = VideoWriter(video_filename, 'MPEG-4');
    v.Quality = 100;  % Set quality from 0 (lowest) to 100 (highest)
    v.FrameRate = 15; % Adjust as needed
    open(v);

    % Get alpha index limits (assuming all alpha vectors are same length)
    n_alpha = length(alpha_gamma_ar);

    figure(3); clf; hold on;
    set(gcf, 'Color', 'w');set(gca, 'Color', 'w'); 
    fill3([-3 3 3 -3], [0 0 0 0], [-3 -3 3 3], [0.2 0.6 0.8], 'FaceAlpha', 0.3); % RGB color with 30% opacity
    colormap(cm);view(az, el);axis off;
    set(gca, 'linewidth', lw - 1);set(gca, 'fontsize', 18);
    xlim([x_min x_max]); ylim([y_min y_max]); zlim([z_min z_max]);
    xticks([x_min x_mid x_max]);yticks([y_min y_mid y_max]);zticks([z_min z_mid z_max]);
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);set(gca,'ZTickLabel',[]);
    plot3(1, 0, 0, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');   % Dot at (0,1,0)
    plot3(-1, 0, 0, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');  % Dot at (0,-1,0)
    set(gcf, 'position', [x0, y0, width, height]);
    for alpha_idx = 2:n_alpha
        for teta_idx = 1:length(teta_intervals)
            beta_gamma_unnorm = surface_data{teta_idx, 1};
            delta_m_gamma_unnorm = surface_data{teta_idx, 2};
            delta_p_gamma_unnorm = -delta_m_gamma_unnorm;

            % Extract 1-row slices and replicate to make them 2D
            z_m = delta_m_gamma_unnorm(1:alpha_idx, :);
            z_p = delta_p_gamma_unnorm(1:alpha_idx, :);
            beta = beta_gamma_unnorm(1:alpha_idx, :);
            % Get alpha value and make Y matrix
            alpha_val = alpha_gamma_ar(1:alpha_idx);

            step = max(10, round(alpha_idx / 30)); % adaptive stride
            index_plot = 1:step:alpha_idx;
            if teta_idx == 2 || teta_idx == 3
                A = -alpha_val;
                surf(beta, A, z_m, 'FaceAlpha', 0.8, 'EdgeColor', 'none');surf(beta, A, z_p, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
                X = [];Y = []; Zm = []; Zp = [];
                w_values = alpha_gamma_ar(index_plot);
                for j = 1:length(index_plot)
                    idx = index_plot(j);
                    X = [X; beta(idx, :) NaN];
                    Y = [Y; -w_values(j) * ones(1, size(beta, 2)) NaN];
                    Zm = [Zm; z_m(idx, :) NaN];
                    Zp = [Zp; z_p(idx, :) NaN];
                end
                % Single batch plot
                [row,~] = size(Y);
                for i = 1:row
                    plot3(X(i,:), Y(i,:), Zm(i,:), 'k', 'LineWidth', lw-2);
                    plot3(X(i,:), Y(i,:), Zp(i,:), 'k', 'LineWidth', lw-2);
                end

            else
                A = alpha_val;
                surf(beta, A, z_m, 'FaceAlpha', 0.8, 'EdgeColor', 'none');surf(beta, A, z_p, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
                X = []; Y = []; Zm = []; Zp = [];
                w_values = alpha_gamma_ar(index_plot);
                for j = 1:length(index_plot)
                    idx = index_plot(j);
                    X = [X; beta(idx, :) NaN];
                    Y = [Y; w_values(j) * ones(1, size(beta, 2)) NaN];
                    Zm = [Zm; z_m(idx, :) NaN];
                    Zp = [Zp; z_p(idx, :) NaN];
                end
                % Single batch plot
                [row,~] = size(Y);
                for i = 1:row
                    plot3(X(i,:), Y(i,:), Zm(i,:), 'k', 'LineWidth', lw-2);
                    plot3(X(i,:), Y(i,:), Zp(i,:), 'k', 'LineWidth', lw-2);
                end
            end

        end
        frame = getframe(gca);
        writeVideo(v, frame);
    end
    close(v);
end

%% FUNCTION: highlight edges of the umbilic
function F = edgelines(beta_gamma_unnorm,delta_p_gamma_unnorm,delta_m_gamma_unnorm,alpha_max_index,lw,direction,alpha_gamma_ar)
w_values = 1:1:10;

% Highlight the segments in the u-v plane for w = -1 and -2
for j = 1:length(w_values)
    index_plot = round(alpha_max_index*j/max(w_values));
    beta = beta_gamma_unnorm(index_plot,:);
    if strcmp(direction, 'minus')
        plot3(beta,-w_values(j) * ones(size(beta)),delta_m_gamma_unnorm(index_plot,:)','k', 'LineWidth', lw-1);
        plot3(beta,-w_values(j) * ones(size(beta)),delta_p_gamma_unnorm(index_plot,:)','k', 'LineWidth', lw-1);
    elseif strcmp(direction, 'plus')
        plot3(beta,w_values(j) * ones(size(beta)),delta_m_gamma_unnorm(index_plot,:)','k', 'LineWidth', lw-1);
        plot3(beta,w_values(j) * ones(size(beta)),delta_p_gamma_unnorm(index_plot,:)','k', 'LineWidth', lw-1);
    end
end

if strcmp(direction, 'plus')
    % Find the minimum value and its index
    [~, minIndex] = min(min(delta_m_gamma_unnorm));
    % Find the maximum value and its index
    [~, maxIndex] = max(max(delta_p_gamma_unnorm));

elseif strcmp(direction, 'minus')
    % Find the minimum value and its index
    [~, minIndex] = min(min(delta_p_gamma_unnorm));
    % Find the maximum value and its index
    [~, maxIndex] = max(max(delta_m_gamma_unnorm));
end

color = [0, 0.3, 1];
% Highlight the edges of the surface with black lines
if strcmp(direction, 'plus')
    plot3(beta_gamma_unnorm(:, end), alpha_gamma_ar(:), delta_m_gamma_unnorm(:,end), 'color',color, 'LineWidth', lw);
elseif strcmp(direction, 'minus')
    plot3(beta_gamma_unnorm(:, 1), alpha_gamma_ar(:), delta_p_gamma_unnorm(:,1),'color',color, 'LineWidth', lw);
end
plot3(beta_gamma_unnorm(:, maxIndex), alpha_gamma_ar(:), delta_p_gamma_unnorm(:,maxIndex), 'color',color, 'LineWidth', lw);
plot3(beta_gamma_unnorm(:, minIndex), alpha_gamma_ar(:), delta_m_gamma_unnorm(:,minIndex), 'color',color, 'LineWidth', lw);

end
