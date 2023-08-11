% Plot a Response Spectrum (red = pseudo-£¬blue = absolution-,default for acceleration with normalization,)
%
% Input:
% 1) wave = Ground Motion (g);
% 2) kesi = Damping ratio
% 3) dt = Sampling period (s)
% 4) abs_psd = Output abs(1) or psd(0), Default: 1,abs
% 5) dT = Natural period Interval (s)
% 6) fig = 0 or 1, Default: 1 = plot the figure
% 7) variable = 'A','V','D' for Acceleration(default),Velocity,Displacement
% 8) normalize = 1,0 for yes,no, Default: 1 for 'A', 0 for 'V' and 'D'
%
% Output:
% 1) Figure S-T;
% 2) peak_psd = pseudo-;
% 3) peak_abs = absolution-


function [T, peak_abs, peak_psd] = responseSpectrum(wave,dt,kesi,abs_psd,T,fig,variable,normalize)
    % Default  parameter
    if nargin < 4
        abs_psd = 1; % Output abs (1), psd (0) or both (other values)
        load("periodPEER.mat")
        T = periodPEER'; % Natural period Interval
        fig = 1; % Default: plot the figure
        variable = 'A'; % Default: acceleration response spectrum
        normalize = 1; % Default: normalization ON 
    elseif nargin < 5
        load("periodPEER.mat")
        T = periodPEER'; % Natural period Interval
        fig = 1; % Default: plot the figure
        variable = 'A'; % Default: acceleration response spectrum
        normalize = 1; % Default: normalization ON 
    elseif nargin < 6
        fig = 1; % Default: plot the figure
        variable = 'A'; % Default: acceleration response spectrum
        normalize = 1; % Default: normalization ON 
    elseif nargin < 7
        variable = 'A'; % Default: acceleration response spectrum
        normalize = 1; % Default: normalization ON
    elseif nargin < 8
        normalize = 0; % Default: normalization OFF
    end

    gg = 9.80; % gravity acceleration
    xg=wave*gg; % transfer from unit g to m/s^2
    xg_peak=max(xg); % PGA
%     T=1e-6: dT: 10; % natural period

    % zero matrix
    var_peak_abs = zeros(1,length(T));
    var_peak_psd = zeros(1,length(T));
        
    for k = 1:1:size(T,2)  % for each natural period
        omiga = 2*pi/T(k); % structural circle frequency

        % SISO state-space model: SDOF
        % state-space model: state variable x = [displacement velocity]'
        A =[0,1;-omiga^2,-2*kesi*omiga]; % state matrix
        B = [0;-1]; % input matrix
        C = eye(2); % output matrix
        D = zeros(2,1); % direct transfer matrix
        sys = ss(A,B,C,D);

        % Time
        n = length(xg); % time step
        t = 0:dt:(n-1)*dt;  % time sequence

        % Response
        y = lsim(sys,xg,t);  %  response of the dynamic system model: 1)PGD; 2) PGV
        dis_peak0 = max(y(:,1)); % peak displacement

        if variable == 'D'
            var_peak_psd = []; % null pseudo for displacement
            var_peak_abs(k) = max(y(:,1));  % peak displacement vector
        elseif variable == 'V'
            var_peak_psd(k) = omiga * dis_peak0;  % pseudo-velocity
            var_peak_abs(k) = max(y(:,2)); % spectral absolute velocity
        elseif variable == 'A'
            var_peak_psd(k) = omiga^2 * dis_peak0;  % pseudo-acceleration
            var_peak_abs(k) = max(-2 * kesi * omiga * y(:,2)...
                - omiga^2 * y(:,1)); % spectral absolute acceleartion
        end
    end

    % Normalization
    if normalize == 1
        peak_psd = (var_peak_psd/xg_peak)'; % normalization
        peak_abs = (var_peak_abs/xg_peak)'; % normalization
        yL = 'alpha';
    else
        peak_psd = var_peak_psd' ./ gg;
        peak_abs = var_peak_abs' ./ gg;
        yL = 'S';
    end
    
    % Figure

    %
    
    % Figure
    if fig == 1
        if size(var_peak_psd,1) == 0 % Displacement
            figure
            plot(T,peak_abs,'-r') 
            title(['Absolute Spectrum' '-' variable])
            xlabel('T(s)')
            ylabel([yL '-' variable])
            legend('Abs')
        else % Velocity or Acceleration
            if abs_psd == 1
                disp('Plot spectra in the abs');
                figure
                plot(T,peak_abs,'-blue') 
                title(['Absolute Spectrum' '-' variable])
                xlabel('T(s)')
                ylabel([yL '-' variable])
                legend('Abs')
            elseif abs_psd == 0
                disp('Plot spectra in the psd');
                figure
                plot(T,peak_psd,'-r') 
                title(['Pseudo Spectrum' '-' variable])
                xlabel('T(s)')
                ylabel([yL '-' variable])
                legend('Pseudo')
            else
                disp('Plot spectra in both the abs and psd');
                figure
                plot(T,peak_psd,'-r',T,peak_abs,'-blue') 
                title(['Pseudo & Absolute Spectrum' '-' variable])
                xlabel('T(s)')
                ylabel([yL '-' variable])
                legend('Pseudo','Abs')
            end
        end
    end
end