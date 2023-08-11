% Plot ground motion A, V, D and Sa curve for certain purposes, such as
% a time point to check the non-convergence at ANSYS
% 
% Input
% accel = acceleration vector;
% dt = sampling time interval;
% select_point = selected point, e.g. for Non-convergence Point (s);
% gm_fig = switch of figure of the gournd motion: 1 = on;
% Sa_fig = switch of figure of Sa: 1 = on;
% 
% e.g. 
% [accel,~,~,~] = getAmpDtPEER(recordFolder,recordName);
% [~,dt,~,~] = getAmpDtPEER(recordFolder,recordName);
% selectPoint = 2.5;
% gmFig = 1;
% SaFig = 1;
% [time,accel,vel,dsp,Sa] = plot_GM_Sa(accel,dt,select_point,...
%     gm_fig,Sa_fig);

function [time,accel,vel,dsp,Sa] = plotGMandSa(accel,dt,selectPoint,...
    gmFig,SaFig)
    
    %%% GM time series and A,V,S
%     [accel, dt, ~, ~] = getAmpDtPEER(recordFolder,recordName);   % Acceleration

    timemax = size(accel,1) * dt;
    time = (0: dt: timemax - dt)';

    vel = cumtrapz(time,accel);
    dsp = cumtrapz(time,vel);   
    
    if gmFig == 1
        %%% Figure
        accel_max = max(abs(accel));
        vel_max = max(abs(vel));
        dsp_max = max(abs(dsp));
        waveFig = figure;
%         set(waveFig,'Position',[100,380,500,600])

        % Acceleration
        subplot(3,1,1);
        plot(time,accel);
        hold on
        for i = 1:length(selectPoint)
            plot([selectPoint(i),selectPoint(i)],[-1.1*accel_max,1.1*accel_max]);
            hold on
        end
        set(gca,'xlim',[0 timemax]);     % x limit
        set(gca,'ylim',[-1.1*accel_max 1.1*accel_max]);     % y limit
        %set(gca,'xticklabel',[]);        % not display x tick label
        set(gca,'box','on');             % box
        legend('','Peak X','Peak Z');          % legend
        ylabel('Acceleration');          % label
        grid on;

        % Velocity
        subplot(3,1,2);
        plot(time,vel);
        hold on
        for i = 1:length(selectPoint)
            plot([selectPoint(i),selectPoint(i)],[-1.1*vel_max,1.1*vel_max]);
            hold on
        end
        set(gca,'xlim',[0 timemax]);     % x limit
        set(gca,'ylim',[-1.1*vel_max 1.1*vel_max]);     % y limit
        set(gca,'box','on');             % box
        legend('','Peak X','Peak Z');          % legend
        ylabel('Velocity');          % label
        grid on;

        % Displacement
        subplot(3,1,3);
        plot(time,dsp);
        hold on
        for i = 1:length(selectPoint)
            plot([selectPoint(i),selectPoint(i)],[-1.1*dsp_max,1.1*dsp_max]);
            hold on
        end
        set(gca,'xlim',[0 timemax]);     % x limit
        set(gca,'ylim',[-1.1*dsp_max 1.1*dsp_max]);     % y limit
        set(gca,'box','on');             % box
        legend('','Peak X','Peak Z');          % legend
        ylabel('Displacement');          % label
        grid on;

        % X axis label
        xlabel('Time');
    end
    
    if SaFig == 1
        %%% Spectrum Response
        dT = 0.01;
        kesi = 0.05;
        [~, Sa, ~] = responseSpectrum(accel,dt,kesi,1,dT,1,'A',0); 
    else
        Sa = NaN;
    end

end