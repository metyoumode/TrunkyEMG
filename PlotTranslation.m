clear
close all
clc

% setting movements parameters
opengl('software');

warning('off')

% Movemnts type
movements = {'FrontFlex','LatFLex'};

% Assistance and Resistance levels
conditions = {'ImpLowAss', 'ImpMediumAss', 'ImpHighAss', 'ResConst', 'ResInteraction'};

% Muscles engaged in EMG analysis
% nameMuscles = {'L.Thoracic Es', 'R.Rect.Abdom.Up.', 'L.Rect.Abdom.Up.', 'R.Ext.Oblique', ...
%                'L.Ext.Oblique','R.Scm', 'L.Scm','L.Lat.Dorsi', 'R.Lat.Dorsi','L.Thoracic Es'};
nameMuscles = {'L.Lat.Dorsi', 'R.Lat.Dorsi', 'R.Thoracic Es', 'L.Thoracic Es', ...
               'R.Rect.Abdom.Up.', 'L.Rect.Abdom.Up.', 'R.Ext.Oblique', 'L.Ext.Oblique', ...
               'R.Scm', 'L.Scm'};
MuscleName={'latDorsi','estensori','abs','RightOblique','Leftoblique','sternocleido'};
numMuscles = length(nameMuscles);

%% Reading .txt data
fs=2000;
%data directory
datadir = fullfile('C:\Users\Mattia\OneDrive - Fondazione Istituto Italiano Tecnologia\Desktop\poli\tesi\Aquisition\Lore_emg', 'txt');

% Structure initialization for storing all TXT data
allTxtData = struct();

% Cycle for all conditions
for c = 1:length(conditions)
    condition = conditions{c};
    
    % Cycle for all angles
    for a = 1:length(movements)
        movement = movements{a};
        
        % Creation of the file .txt name
        fname = sprintf('%s%s.txt', condition, movement);
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

        % Imposta i nomi delle variabili con i nomi dei muscoli
        data.Properties.VariableNames = nameMuscles;

        % Numero di campioni e tempo
        nsamples = height(data); % Usa height invece di length per le tabelle
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
    nsamples = height(data); % Usa height invece di length per le tabelle
    time = (0:(nsamples-1)) / fs; % fc deve essere definito, viene dai file .c3d

    if strcmp(MuscleName{i},'latDorsi')
        muscleData1 = data{:, 1};
        muscleData2 = data{:, 2};
        mvc(1) = max(muscleData1);
        mvc(2) = max(muscleData2);

    elseif strcmp(MuscleName{i},'abs')
        muscleData1 = data{:, 5};
        muscleData2 = data{:, 6};
        mvc(5) = max(muscleData1);
        mvc(6) = max(muscleData2);
    elseif strcmp(MuscleName{i},'estensori')
        muscleData1 = data{:, 3};
        muscleData2 = data{:, 4};
        mvc(3) = max(muscleData1);
        mvc(4) = max(muscleData2);
    elseif strcmp(MuscleName{i},'sternocleido')
        muscleData1 = data{:, 9};
        muscleData2 = data{:, 10};
        mvc(9) = max(muscleData1);
        mvc(10) = max(muscleData2);
    elseif strcmp(MuscleName{i},'RightOblique')
        muscleData1 = data{:, 7};
        mvc(7) = max(muscleData1);
    elseif strcmp(MuscleName{i},'LeftOblique')
        muscleData1 = data{:, 8};
        mvc(8) = max(muscleData1);
    end



end
%% RMS evelope and mvc normalization
close all
%FIR filter 
fs=2000;
% Define parameters
window_size_ms = 300;  % Window size in milliseconds
low_cutoff = 30;  % Frequenza di taglio inferiore 
high_cutoff = 700;  % Frequenza di taglio superiore
notch_freq = 50;  % Frequenza di disturbo (es. 50 Hz per l'interferenza di rete)
notch_bandwidth = 5;  % Larghezza di banda del filtro notch (5 Hz)

% Progettazione del filtro passa banda Butterworth
[b_BP, a_BP] = butter(6, [low_cutoff, high_cutoff] / (fs / 2), 'bandpass');
% [b_LP, a_LP] = butter(6, high_cutoff / (fs / 2), 'low');

% Progettazione del filtro notch
[d_Notch, c_Notch] = iirnotch(notch_freq / (fs / 2), notch_bandwidth / (fs / 2));

% Convert window size to samples
window_size_samples = round(window_size_ms * fs / 1000);

% Filtraggio e plot dei dati
figure;
%signal_norm=zeros(104032,10);
StartTime=zeros(1,10);
value=zeros(1,10);
p=0;
% Iterate through each condition
for c = 1:length(conditions)
    condition = conditions{c};

    % Iterate through each movement
    for a = 1:length(movements)
        movement = movements{a};

        % Access the raw EMG data for this condition and angle
        emgData = allTxtData.(condition).(movement).data;
        rec_signal = zeros(size(emgData,1), numMuscles);

        time = allTxtData.(condition).(movement).time;

        % Create a new figure for each condition and angle
        figure('Name', sprintf('%s - %s', condition, movement), 'NumberTitle', 'off');
        n=7; %prendo il right oblique
        raw_signal = emgData{:, n};

        % Filtraggio del segnale con passa banda Butterworth
        %         filtered_signal = filtfilt(b_LP, a_LP, raw_signal);
        filtered_signal = filtfilt(b_BP, a_BP, raw_signal);

        % Applicazione del filtro notch
        filtered_signal = filtfilt(d_Notch, c_Notch, filtered_signal);

        % Full wave rectification
        rec_signal(:,n) = abs(filtered_signal);

        % RMS envelope
        % singleRec_signal=rec_signal(:,m);
        % Step 1: Square the EMG signal
        squared_signal =rec_signal .^2;

        % Step 2: Apply the moving average on the squared signal
        rms_moving_avg = movmean(squared_signal, window_size_samples);

        % Step 3: Take the square root to get the RMS envelope
        rms_envelope = sqrt(rms_moving_avg);

        % Normalizatio to MVC
        signal_norm = (rms_envelope/mvc(n))*1e4;
        signal_norm=round(signal_norm);
        campioni=round(time*fs);
        % Plot RMS envelope normlized to MVC
        plot(campioni, signal_norm(:,n), 'Color', "#A2142F");
%         [StartTime(p),value(p)]=ginput(signal_norm(:,n));
%         close(gcf);
        xlabel('Time [s]');
        ylabel('Signal []');
        legend({'Normalized to MVC'});
        grid on;
        p=p+1;
      
    end

    % Add a main title for the entire figure
    sgtitle(sprintf('Condition: %s, movement: %s', condition, movement));
end

