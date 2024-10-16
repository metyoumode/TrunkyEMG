clear all;
close all;
clc;

% setting movements parameters
opengl('software');

warning('off')

% Movemnts type
movements = {'FrontFlex','LatFLex'};

% Assistance and Resistance levels
conditions = {'ImpLowAss', 'ImpMediumAss', 'ImpHighAss', 'ResConst', 'ResInteraction'};

% Muscles engaged in EMG analysis
nameMuscles = {'L.Thoracic Es', 'R.Rect.Abdom.Up.', 'L.Rect.Abdom.Up.', 'R.Ext.Oblique', ...
               'L.Ext.Oblique','R.Scm', 'L.Scm','L.Lat.Dorsi', 'R.Lat.Dorsi','R.Thoracic Es'};
% nameMuscles = {'L.Lat.Dorsi', 'R.Lat.Dorsi', 'R.Thoracic Es', 'L.Thoracic Es', ...
%                'R.Rect.Abdom.Up.', 'L.Rect.Abdom.Up.', 'R.Ext.Oblique', 'L.Ext.Oblique', ...
%                'R.Scm', 'L.Scm'};
MuscleName={'latDorsi','estensori','abs','RightOblique','LeftOblique','sternocleido'};
numMuscles = length(nameMuscles);

%% crd reading

%data directory
datadir = fullfile('C:\Users\Mattia\OneDrive - Fondazione Istituto Italiano Tecnologia\Desktop\poli\tesi\Aquisition\Dave_emg', 'c3d');

% Structure initialization for storing all C3D data
allData = struct();

% Cycle for all the assistance/resistance conditions
for c = 1:length(conditions)
    condition = conditions{c};
    
    % Cycle for all the movements 
    for a = 1:length(movements)
        movement = movements{a};
        
        % Creation of the file name
        fname = sprintf('%s%s_dave.c3d', condition, movement);
        
        % Read the .c3d file
        data = readc3d(fullfile(datadir, fname));
        
        % Parameters extraction
        fs = data.Parameter.ANALOG.RATE.data; %prendo solo la frequenza di camponamento
        %il resto stica
        % nbits = data.Parameter.ANALOG.BITS.data;
        % gain = data.Parameter.ANALOG.GAIN.data;
        % unit = data.Parameter.ANALOG.UNITS.data'; 
        
        % Store the data in the structure allData
        % The structure is indexed for condition and movements
        allData.(condition).(movement) = struct('data', data,'fs', fs);
    end
end
% Now, you can access to data for a certain condition and angle using:
% allData.<condition>.Angle_<angle>.data, ad esempio: allData.AssBassa.Angle_30.data

%% Reading .txt data

%data directory
datadir = fullfile('C:\Users\Mattia\OneDrive - Fondazione Istituto Italiano Tecnologia\Desktop\poli\tesi\Aquisition\Dave_emg', 'txt');

% Structure initialization for storing all TXT data
allTxtData = struct();

% Cycle for all conditions
for c = 1:length(conditions)
    condition = conditions{c};
    
    % Cycle for all angles
    for a = 1:length(movements)
        movement = movements{a};
        
        % Creation of the file .txt name
        fname = sprintf('%s%s_dave.txt', condition, movement);
        filepath = fullfile(datadir, fname);

        % Apri il file
        fileID = fopen(filepath, 'r');     
        
        % Salta le prime quattro righe che sono i muscoli
        for i = 1:4
            fgetl(fileID);
        end

        % Leggi il resto del file in una tabella
        data = readtable(filepath);

        % Chiudi il file
        fclose(fileID);

        % Remove the 1st column that is the time
        data = data(:, 2:11);
        %remove extra raws
        data=data(1:73535,:);
       

        % Imposta i nomi delle variabili con i nomi dei muscoli
        data.Properties.VariableNames = nameMuscles;

        % Numero di campioni e tempo
        nsamples =height(data); % Usa height invece di length per le tabelle
        time = (0:(nsamples-1)) / fs;

        % Memorizzazione dei dati nella struttura allTxtData
        allTxtData.(condition).(movement) = struct('data', data, ...
            'nsamples', nsamples, ...
            'time', time, ...
            'fs', fs);
    end
end
% Now, you can access to data for a certain condition and angle using:
% allTxtData.<condition>.Angle_<angle>.data, ad esempio: allTxtData.AssBassa.Angle_30.data

%% Reading MVC value

%qui non viene fatto ma probabilmente non sarebbe male filtrarlo prima

% Initialize a row vector for storing MVC values
%mvc value must be the same number or the muscles
mvc = zeros(1, numMuscles);
% Load MVC data for each muscle
for i = 1:length(MuscleName) %I will read less mvc than muscles, since we 
    %simmetrical muscles with the same data for the mvc
    fname = sprintf('mvc_%s.txt', MuscleName{i});
    filepath = fullfile(datadir, fname);

    % Apri il file
    fileID = fopen(filepath, 'r');

    % Salta le prime quattro righe
    for j = 1:4
        fgetl(fileID);
    end

    % Leggi il resto del file in una tabella
    data = readtable(filepath);

    % Chiudi il file
    fclose(fileID);

    % Remove the 1st column that is the time
    data = data(:, 2:11);
 

  
    % Numero di campioni e tempo
  nsamples = 73535; % Usa height invece di length per le tabelle
    time = (0:(nsamples-1)) / fs; % fc deve essere definito, viene dai file .c3d

    if strcmp(MuscleName{i},'latDorsi')
        muscleData1 = data{:, 8};
        muscleData2 = data{:, 9};
        mvc(1) = max(muscleData1);
        mvc(2) = max(muscleData2);

    elseif strcmp(MuscleName{i},'abs')
        muscleData1 = data{:, 2};
        muscleData2 = data{:, 3};
        mvc(5) = max(muscleData1);
        mvc(6) = max(muscleData2);
    elseif strcmp(MuscleName{i},'estensori')
        muscleData1 = data{:, 1};
        muscleData2 = data{:, 10};
        mvc(3) = max(muscleData1);
        mvc(4) = max(muscleData2);
    elseif strcmp(MuscleName{i},'sternocleido')
        muscleData1 = data{:, 6};
        muscleData2 = data{:, 7};
        mvc(9) = max(muscleData1);
        mvc(10) = max(muscleData2);
    elseif strcmp(MuscleName{i},'RightOblique')
        muscleData1 = data{:, 4};
        mvc(7) = max(muscleData1);
    elseif strcmp(MuscleName{i},'LeftOblique')
        muscleData1 = data{:, 5};
        mvc(8) = max(muscleData1);
    end



end

%% Plot Raw Data

clear fileID filepath fname angle datadir 
figure;

% Iterate through each condition
for c = 1:length(conditions)
    condition = conditions{c};
    
    % Iterate through each angle
    for a = 1:length(movements)
        movement = movements{a};
        
        % Access the raw EMG data for this condition and angle
        emgData = allTxtData.(condition).(movement).data;
        time = allTxtData.(condition).(movement).time;
        
        % Create a new figure for each condition and angle
        figure('Name', sprintf('%s - %s', condition, movement), 'NumberTitle', 'off');
        
        % Plot each muscle signal in a subplot
        for m = 1:numMuscles
            subplot(4, 4, m);
            plot(time, emgData{:, m}, "Color", "#7E2F8E");
            title(nameMuscles{m}, 'Interpreter', 'none');
            xlabel('Time [s]');
            ylabel('Signal [\muV]');
            grid on;
        end
        
        % Add a main title for the entire figure
        sgtitle(sprintf('Condition: %s, Movement: %s', condition, movement));
    end
end

%% FFT
fs = 2000;
f_FFT = fs*(0:(nsamples/2))/nsamples;


% Iterate through each condition
for c = 1:length(conditions)
    condition = conditions{c};    
    % Iterate through each angle
    for a = 1:length(movements)
        movement = movements{a};
        
        % Access the raw EMG data for this condition and angle
        emgData = allTxtData.(condition).(movement).data;
        time = allTxtData.(condition).(movement).time;
        
        % Create a new figure for each condition and angle
        figure('Name', sprintf('%s - %s', condition, movement), 'NumberTitle', 'off');
        
        % Plot each muscle signal in a subplot
        for m = 1:numMuscles
            % Compute the 2-sided spectrum [-Fmax:+Fmax]
            p1 = fft(emgData{:,m});

            p1 = abs(p1/nsamples);
            
            % Compute the single-sided spectrum by taking the positive part
            % of the 2-sided spectrum and multiply by 2
            p1 = p1(1:nsamples/2+1);
            p1(2:end-1) = 2*p1(2:end-1);
        
            subplot(4,4,m)
            plot(f_FFT,p1);
            set(gcf, 'Position', [100, 100, 1600, 900]);  % Aumenta la dimensione della figura

            title(nameMuscles{m}, 'Interpreter', 'none');        
            xlabel('Frequency [Hz]');
            ylabel('Intensity');
            hold on;
            grid on;
        end
        
        % Add a main title for the entire figure
        sgtitle(sprintf('Condition: %s, Movement: %s', condition, movement));
        
    end
end

%% Parametri del filtro
fs = 2000;  % Frequenza di campionamento
low_cutoff = 10;  % Frequenza di taglio inferiore 
high_cutoff = 500;  % Frequenza di taglio superiore
notch_freq = 50;  % Frequenza di disturbo (es. 50 Hz per l'interferenza di rete)
notch_bandwidth = 5;  % Larghezza di banda del filtro notch (5 Hz)

% Progettazione del filtro passa banda Butterworth
[b_BP, a_BP] = butter(6, [low_cutoff, high_cutoff] / (fs / 2), 'bandpass');
% [b_LP, a_LP] = butter(6, high_cutoff / (fs / 2), 'low');

% Progettazione del filtro notch
[d_Notch, c_Notch] = iirnotch(notch_freq / (fs / 2), notch_bandwidth / (fs / 2));

% Filtraggio e plot dei dati
figure;

% Iterate through each condition
for c = 1:length(conditions)
    condition = conditions{c};
    
    % Iterate through each angle
    for a = 1:length(movements)
        movement = movements{a};
        
        % Access the raw EMG data for this condition and angle
        emgData = allTxtData.(condition).(movement).data;
        rec_signal = zeros(size(emgData,1), numMuscles);

        time = allTxtData.(condition).(movement).time;        
        
        % Create a new figure for each condition and angle
        figure('Name', sprintf('%s - %s', condition, movement), 'NumberTitle', 'off');
        
        % Plot each muscle signal in a subplot
        for m = 1:numMuscles

                % Segnale EMG grezzo
                raw_signal = emgData{:, m};
                
                 if m==2 || m==3
                % Progettazione del filtro passa banda Butterworth
                [b_BP, a_BP] = butter(6, [35, 200] / (fs / 2), 'bandpass');

                % Progettazione del filtro notch
                [d_Notch, c_Notch] = iirnotch(notch_freq / (fs / 2), 5 / (fs / 2));
                [num,den] = iirnotch(1.4 / (fs / 2), 20 / (fs / 2));

                filtered_signal = filtfilt(b_BP, a_BP, raw_signal);

                % Applicazione del filtro notch
                filtered_signal = filtfilt(d_Notch, c_Notch, filtered_signal);
                filtered_signal=filtfilt(num,den,filtered_signal);
            else
                % Filtraggio del segnale con passa banda Butterworth
                %                 filtered_signal = filtfilt(b_LP, a_LP, raw_signal);
                filtered_signal = filtfilt(b_BP, a_BP, raw_signal);

                % Applicazione del filtro notch
                filtered_signal = filtfilt(d_Notch, c_Notch, filtered_signal);

            end
                
                % Full wave rectification
                rec_signal(:,m) = abs(filtered_signal);

                subplot(4, 4, m)

                % Plot RMS envelope normalized to MVC
                plot(time, raw_signal, 'Color', "#7E2F8E", 'LineWidth', 1);  % Increase line width
                hold on

                plot(time, filtered_signal, 'Color', "#77AC30", 'LineWidth', 1);  % Increase line width

                title(nameMuscles{m}, 'Interpreter', 'none');
                xlabel('Time [s]');
                ylabel('Signal [\muV]');

                % Create the legend and reduce its font size
                %legend({'Raw', 'Filtered'}, 'FontSize', 8);  % Decrease font size

                grid on;
        end

        % Add a main title for the entire figure
        sgtitle(sprintf('Condition: %s, Movement: %s', condition, movement));
    end
end

%% RMS evelope and mvc normalization
%FIR filter 

% Define parameters
window_size_ms = 300;  % Window size in milliseconds

% Convert window size to samples
window_size_samples = round(window_size_ms * fs / 1000);

% Filtraggio e plot dei dati
figure;
%signal_norm=zeros(104032,10);
% Iterate through each condition
for c = 1:length(conditions)
    condition = conditions{c};
    
    % Iterate through each angle
    for a = 1:length(movements)
        movement = movements{a};
        
        % Access the raw EMG data for this condition and angle
        emgData = allTxtData.(condition).(movement).data;
        rec_signal = zeros(size(emgData,1), numMuscles);

        time = allTxtData.(condition).(movement).time;        
        
        % Create a new figure for each condition and angle
        figure('Name', sprintf('%s - %s', condition, movement), 'NumberTitle', 'off');
        
        % Plot each muscle signal in a subplot
        for m = 1:numMuscles
             % Segnale EMG grezzo
                raw_signal = emgData{:, m};
                
                   if m==2 || m==3
                % Progettazione del filtro passa banda Butterworth
                [b_BP, a_BP] = butter(6, [35, 200] / (fs / 2), 'bandpass');

                % Progettazione del filtro notch
                [d_Notch, c_Notch] = iirnotch(notch_freq / (fs / 2), 5 / (fs / 2));
                [num,den] = iirnotch(1.4 / (fs / 2), 20 / (fs / 2));

                filtered_signal = filtfilt(b_BP, a_BP, raw_signal);

                % Applicazione del filtro notch
                filtered_signal = filtfilt(d_Notch, c_Notch, filtered_signal);
                filtered_signal=filtfilt(num,den,filtered_signal);
            else
                % Filtraggio del segnale con passa banda Butterworth
                %                 filtered_signal = filtfilt(b_LP, a_LP, raw_signal);
                filtered_signal = filtfilt(b_BP, a_BP, raw_signal);

                % Applicazione del filtro notch
                filtered_signal = filtfilt(d_Notch, c_Notch, filtered_signal);

            end
                
                % Full wave rectification
                rec_signal(:,m) = abs(filtered_signal);
                
              % RMS envelope
             % singleRec_signal=rec_signal(:,m);
                % Step 1: Square the EMG signal
                squared_signal =rec_signal .^2;
                
                % Step 2: Apply the moving average on the squared signal
                rms_moving_avg = movmean(squared_signal, window_size_samples);
                
                % Step 3: Take the square root to get the RMS envelope
                rms_envelope = sqrt(rms_moving_avg);
        
                % Normalizatio to MVC
                signal_norm = rms_envelope/mvc(m);
                        
                subplot(4,4,m)

                % Plot RMS envelope normlized to MVC
                plot(time, signal_norm(:,m), 'Color', "#A2142F");                
                
                title(nameMuscles{m}, 'Interpreter', 'none');
                xlabel('Time [s]');
                ylabel('Signal []');
                legend({'Normalized to MVC'});
                grid on;
        end
        
        % Add a main title for the entire figure
        sgtitle(sprintf('Condition: %s, movement: %s', condition, movement));
    end
end

% %% Plot dei confronti delle varie modalità
% 
% % Definisci le condizioni specifiche
% specificConditions = {'LowAss','MediumAss','HighAss','ConstRes','InteractionRes'};
% 
% % Iterate through each angle
% for a = 1:length(movements)
%     movement = movements{a};
%     
%     % Create a new figure for each angle
%     figure('Name', sprintf('Comparison - %s', movement), 'NumberTitle', 'off');
%     
%     % Iterate through each muscle to create a subplot
%     for m = 1:numMuscles
%         subplot(4, 4, m);
%         hold on; % Hold the plot to overlay signals from different conditions
%         
%         % Iterate through each specific condition
%         for c = 1:length(conditions)
%             condition = conditions{c};
%             
%             % Access the raw EMG data for this condition and angle
%             emgData = allTxtData.(condition).(movement).data;
%             time = allTxtData.(condition).(movement).time;
%             
%             % Segnale EMG grezzo
%             raw_signal = emgData{:, m};
%             
%             % Filtraggio del segnale con passa banda Butterworth
%             filtered_signal = filtfilt(b_BP, a_BP, raw_signal);
%     
%             % Applicazione del filtro notch
%             filtered_signal = filtfilt(d_Notch, c_Notch, filtered_signal);
%             
%             % Full wave rectification
%             rec_signal = abs(filtered_signal);
%             
%             % RMS envelope
%             squared_signal = raw_signal.^2;
%             rms_moving_avg = movmean(squared_signal, window_size_samples);
%             rms_envelope = sqrt(rms_moving_avg);
%     
%             % Normalization to MVC
%             signal_norm = rms_envelope / mvc(m);
%                     
%             % Plot RMS envelope normalized to MVC
%             plot(time, signal_norm, 'DisplayName', condition);
%         end
%         
%         title(nameMuscles{m}, 'Interpreter', 'none');
%         xlabel('Time [s]');
%         ylabel('Signal []');
%         grid on;
%         hold off; % Release the plot for the next subplot
%     end
%     
%     % Add a main title for the entire figure
%     sgtitle(sprintf('Comparison of Conditions - Movement: %s', movement));
% end
%% Plot dei confronti delle varie modalità con traslazione

% Definisci le condizioni specifiche
specificConditions = {'LowAss','MediumAss','HighAss','ConstRes','InteractionRes'};
start=0;
 entrato=0;
% Iterate through each angle
for a = 1:length(movements)
    movement = movements{a};
    
    % Create a new figure for each angle
    figure('Name', sprintf('Comparison - %s', movement), 'NumberTitle', 'off');
  

    % Iterate through each muscle to create a subplot
    for m = 1:numMuscles
        subplot(4, 4, m);
        hold on; % Hold the plot to overlay signals from different conditions
        
        % Iterate through each specific condition
        for c = 1:length(conditions)
            condition = conditions{c};

            if strcmp(movement,'LatFlex')
                if strcmp(condition,'ImpLowAss')
                    start=32000-6000;

                elseif strcmp(condition,'ImpMediumAss')
                    start= 26000-6000;

                elseif strcmp(condition,'ImpHighAss')
                    start=36000-6000;


                elseif strcmp(condition,'ResConst')
                    start=26000-6000;

                elseif strcmp(condition,'ResInteraction')
                    start=26000-6000;

                end

            elseif strcmp(movement,'FrontFlex')

                if strcmp(condition,'ImpLowAss')
                    start=22000-6000;
                elseif strcmp(condition,'ImpMediumAss')
                    start=20000-6000;
                elseif strcmp(condition,'ImpHighAss')
                    start=20000-6000;
                elseif strcmp(condition,'ResConst')
                    start=20000-6000;
                elseif strcmp(condition,'ResInteraction')
                    start= 18000-6000;
                end
            end

            % Access the raw EMG data for this condition and angle
            emgData = allTxtData.(condition).(movement).data;
            time = allTxtData.(condition).(movement).time;
            
            % Segnale EMG grezzo
            raw_signal = emgData{:, m};
            
            % Filtraggio del segnale con passa banda Butterworth
            filtered_signal = filtfilt(b_BP, a_BP, raw_signal);
    
            % Applicazione del filtro notch
            filtered_signal = filtfilt(d_Notch, c_Notch, filtered_signal);
            
            % Full wave rectification
            rec_signal = abs(filtered_signal);
            
            % RMS envelope
            squared_signal = raw_signal.^2;
            rms_moving_avg = movmean(squared_signal, window_size_samples);
            rms_envelope = sqrt(rms_moving_avg);
    
            % Normalization to MVC
            signal_norm = rms_envelope / mvc(m);
            signal_norm=signal_norm(start:end);
            Time=linspace(0,time(end),length(signal_norm));
            % Plot RMS envelope normalized to MVC
            plot(Time, signal_norm, 'DisplayName', condition);
        end
        
        title(nameMuscles{m}, 'Interpreter', 'none');
        xlabel('Time [s]');
        ylabel('Signal []');
        legend('show');
        grid on;
        hold off; % Release the plot for the next subplot
    end
    
    % Add a main title for the entire figure
    sgtitle(sprintf('Comparison of Conditions - Movement: %s', movement));
end