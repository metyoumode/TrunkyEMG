close all
clear
clc

%% lettura dati EMG
% setting movements parameters
opengl('software');

warning('off')

% Movemnts type
movements = {'FrontFlex','LatFlex'};

% Assistance and Resistance levels
conditions = {'ImpHighAss', 'ImpMediumAss', 'ImpLowAss', 'ResConst', 'ResInteraction'};

% Muscles engaged in EMG analysis
nameMuscles = {'L.Erector Spinae', 'R.Rect.Abdom.Up.', 'L.Rect.Abdom.Up.', 'R.Ext.Oblique', ...
               'L.Ext.Oblique','R.Scm', 'L.Scm','L.Lat.Dorsi', 'R.Lat.Dorsi','R.Erector Spinae'};
% nameMuscles = {'L.Lat.Dorsi', 'R.Lat.Dorsi', 'R.Thoracic Es', 'L.Thoracic Es', ...
%                'R.Rect.Abdom.Up.', 'L.Rect.Abdom.Up.', 'R.Ext.Oblique', 'L.Ext.Oblique', ...
%                'R.Scm', 'L.Scm'};
MuscleName={'latDorsi','estensori','abs','RightOblique','LeftOblique','sternocleido'};
numMuscles = length(nameMuscles);
fs= 2000;

%% loading speedgoat load
Log=load("ModeLog.mat");

%saving the starting time instant of the different exercises
State=Log.logsout{11}.Values;

 % 1 è High imp, 2 medium, 3 low
StartFrontFlex=[0 0 0];
EndFrontFlex=[ 0 0 0];
StartLatFlex=[ 0 0 0];
EndLatFlex=[ 0 0 0];

repetition=0;
j=1;
i=1;
startExercise=1;
entrato=0;

for t=1:length(State.Data)
    %se ho iniziato una flessione e sono sotto il numero di ripetizioni
    if (State.Data(t) == 18 && startExercise==1 )
       StartFrontFlex(j)=t-3/0.002;  
        EndFrontFlex(j)=t+((12+2)/0.002);
        j=j+1;
        startExercise=0;
        repetition=repetition+1;
    end
    %se ho iniziato una flessione laterale e sono sotto il numero di ripetizioni
    if (State.Data(t) == 19 && startExercise==1 )
        StartLatFlex(i)=t-3/0.002;
        EndLatFlex(i)=t+((12+2)/0.002);
        i=i+1;
        startExercise=0;
        repetition=repetition+1;
    end
     % if the exercise is finished we enable the starting of the ext one
    if State.Data(t)== 22 
        startExercise=1;
        entrato=1;
    end

end
%% reading EMG data
%data directory
datadir = fullfile('C:\Users\Mattia\OneDrive - Fondazione Istituto Italiano Tecnologia\Desktop\poli\tesi\Aquisition\Mode_emg', 'txt');

% Structure initialization for storing all TXT data
allTxtData = struct();


% Cycle for all conditions
for c = 1:length(conditions)
    condition = conditions{c};
    
    % Cycle for all angles
    for a = 1:length(movements)
        movement = movements{a};
        
        % Creation of the file .txt name
        fname = sprintf('%s%s_mode.txt', condition, movement);
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

%% plot dati 
   % definisco i filtri
close all
clc
low_cutoff = 25;  % Frequenza di taglio inferiore 
high_cutoff = 400;  % Frequenza di taglio superiore
notch_freq = 50;  % Frequenza di disturbo (es. 50 Hz per l'interferenza di rete)
notch_bandwidth = 5;  % Larghezza di banda del filtro notch (5 Hz)

% Progettazione del filtro passa banda Butterworth
[b_BP, a_BP] = butter(6, [low_cutoff, high_cutoff] / (fs / 2), 'bandpass');
% [b_LP, a_LP] = butter(6, high_cutoff / (fs / 2), 'low');

% Progettazione del filtro notch
[d_Notch, c_Notch] = iirnotch(notch_freq / (fs / 2), notch_bandwidth / (fs / 2));

% Define parameters
window_size_ms = 300;  % Window size in milliseconds

% Convert window size to samples
window_size_samples = round(window_size_ms * fs / 1000);
%variabili per i giunti
joints=[1 2 3 4 5];
figure;

%iterate through each joint
for j=4:length(joints)
    joint=joints(j);
    
   pos=Log.logsout{17}.Values.(sprintf('PDO_M%d',joint)).(sprintf('pos_actual_%d',joint));
    %torque=Log.logsout{17}.Values.PDO_M1.torque_actual_1;

  
        % Iterate through each angle
        for a = 1:length(movements)
            movement = movements{a};

            % Access the raw EMG data for this condition and angle
            emgData = allTxtData.(condition).(movement).data;
            % inizialize the matrix to have one column of data for each
            % muscle
            rec_signal = zeros(size(emgData,1), numMuscles);

            %prendo il giusto start per ogni movimento
        
             
          %  time = allTxtData.(condition).(movement).time;

            % Create a new figure for each condition and angle
            figure('Name', sprintf('%s - %s', condition, movement), 'NumberTitle', 'off');

            % Plot each muscle signal in a subplot
            for m = 1:numMuscles
                
                        subplot(4,4,m)
                        hold on
                        grid on
                %Iterate through each condition
                for c = 1:3
                    condition = conditions{c};
                    if strcmp(movement,'LatFlex')
                        if strcmp(condition,'ImpLowAss')
                            start=26000-6000;
                        elseif strcmp(condition,'ImpMediumAss')
                            start= 32000-6000+2000;
                        elseif strcmp(condition,'ImpHighAss')
                            start=28000-6000+2000;
                        end

                    elseif strcmp(movement,'FrontFlex')

                        if strcmp(condition,'ImpLowAss')
                            start=17000-6000;
                        elseif strcmp(condition,'ImpMediumAss')
                            start=17600-6000;
                        elseif strcmp(condition,'ImpHighAss')
                            start=13700-6000-4000;
                        end
                    end
                    % Segnale EMG grezzo for one muscle
                    raw_signal = emgData{:, m};

                    if m==2 || m==3 %condition for filtering abs
                        % Progettazione del filtro passa banda Butterworth
                        [b_BP, a_BP] = butter(6, [35, 200] / (fs / 2), 'bandpass');

                        % Progettazione del filtro notch
                        [d_Notch, c_Notch] = iirnotch(notch_freq / (fs / 2), 5 / (fs / 2));
                        [num,den] = iirnotch(1.4 / (fs / 2), 20 / (fs / 2));

                        filtered_signal = filtfilt(b_BP, a_BP, raw_signal);

                        % Applicazione del filtro notch
                        filtered_signal = filtfilt(d_Notch, c_Notch, filtered_signal);
                        %second notch filter
                        filtered_signal=filtfilt(num,den,filtered_signal);

                    else
                        % Filtraggio del segnale con passa banda Butterworth
                        filtered_signal = filtfilt(b_BP, a_BP, raw_signal);

                        % Applicazione del filtro notch
                        filtered_signal = filtfilt(d_Notch, c_Notch, filtered_signal);

                    end
                    % Full wave rectification
                    rec_signal(:,m) = abs(filtered_signal);

                    % RMS envelope
                    % Step 1: Square the EMG signal
                    squared_signal =rec_signal .^2;

                    % Step 2: Apply the moving average on the squared signal
                    rms_moving_avg = movmean(squared_signal, window_size_samples);

                    % Step 3: Take the square root to get the RMS envelope
                    rms_envelope = sqrt(rms_moving_avg);

                    % Normalizatio to MVC
                    signal_norm = rms_envelope/mvc(m);
                    %take only useful signal
                    signal_norm=signal_norm(start:end,:);

                    %different plot depending on the movement
                    if a==1

                        movementSample=(EndFrontFlex(c)-StartFrontFlex(c));
                        scalingFactor=floor(length(signal_norm)/movementSample);
                        signal_norm_reduce=downsample(signal_norm,scalingFactor);
                        %converting inc into degree
                        posDeg=(pos.Data(StartFrontFlex(c):EndFrontFlex(c)))/5461;

                        % Plot RMS envelope normlized to MVC
                        
                        plot(posDeg, signal_norm_reduce(1:(EndFrontFlex(c)-StartFrontFlex(c)+1),m), 'DisplayName', condition);
                       hold on
                    else

                        movementSample=(EndLatFlex(c)-StartLatFlex(c));
                        scalingFactor=floor(length(signal_norm)/movementSample);
                        signal_norm_reduce=downsample(signal_norm,scalingFactor);
                        posDeg=(pos.Data(StartLatFlex(c):EndLatFlex(c)))/5461;
                        subplot(4,4,m)

                        % Plot RMS envelope normlized to MVC
                        
                        plot(posDeg, signal_norm_reduce(1:(EndLatFlex(c)-StartLatFlex(c)+1),m),  'DisplayName', condition);
hold on

                        %legend({'Normalized to MVC'});
                    end
                end
                title(nameMuscles{m}, 'Interpreter', 'none');
                 xlabel('Position [inc]', 'FontSize', 12);
                    ylabel('Amplitude', 'FontSize', 12);
                   set(gca, 'FontSize', 12);  % Set font size for tick labels
                hold off
           end
            % Add a main title for the entire figure
            sgtitle(sprintf('Movement: %s Joint: %d', movement,joint));
            %legend('show', 'Location', 'northeastoutside');
            % Save the figure as an SVG file
            %filename=sprintf('%s%s_envelope.svg', condition, movement);
            %filename = 'plot_image.svg'; % Name of the file
        end

end