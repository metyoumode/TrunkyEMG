clc
close all
clear


%sampling frequency
fs=2000;
% Parametri dei filtri
low_cutoff = 25;   % Frequenza di taglio bassa del filtro passa-banda
high_cutoff = 400; % Frequenza di taglio alta del filtro passa-banda


% Creazione del filtro passa-banda Butterworth (ordine 6)
[b_BP, a_BP] = butter(6, [low_cutoff, high_cutoff] / (fs / 2), 'bandpass');

% Convertiamo i coefficienti del filtro in oggetti di sistema per poter usare bode
sys_BP = tf(b_BP, a_BP, 1/fs); % Filtro passa banda Butterworth

% Plot del Bode per il filtro Butterworth con asse x in Hz
figure('Name', 'Bode Plot - Butterworth Bandpass Filter', 'NumberTitle', 'off');
[mag_BP, phase_BP, wout_BP] = bode(sys_BP);  % Ottenere dati Bode

% Convertire la frequenza da rad/s a Hz
freq_Hz_BP = wout_BP / (2 * pi);

% Grafico in Hz
subplot(2,1,1);
semilogx(freq_Hz_BP, squeeze(20*log10(mag_BP))); % Magnitude in dB
title('Magnitude');
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
grid on;

subplot(2,1,2);
semilogx(freq_Hz_BP, squeeze(phase_BP)); % Fase in gradi
title('Phase');
xlabel('Frequency [Hz]');
ylabel('Phase [°]');
grid on;
sgtitle('Butterworth Band-Pass Filter (25-400 Hz)')

%% notch 50 HZ

notch_freq = 50;  % Frequenza del disturbo (es. 50 Hz per il rumore della linea elettrica)
notch_bandwidth = 5;  % Larghezza di banda del filtro Notch (es. 5 Hz)


% Creazione del filtro Notch
[d_Notch, c_Notch] = iirnotch(notch_freq / (fs / 2), notch_bandwidth / (fs / 2));

% Convertiamo i coefficienti del filtro in oggetti di sistema per poter usare bode
sys_Notch = tf(d_Notch, c_Notch, 1/fs); % Filtro Notch


% Plot del Bode per il filtro Notch con asse x in Hz
figure('Name', 'Bode Plot - Notch Filter', 'NumberTitle', 'off');
[mag_Notch, phase_Notch, wout_Notch] = bode(sys_Notch);  % Ottenere dati Bode

% Convertire la frequenza da rad/s a Hz
freq_Hz_Notch = wout_Notch / (2 * pi);

% Grafico in Hz
subplot(2,1,1);
semilogx(freq_Hz_Notch, squeeze(20*log10(mag_Notch))); % Magnitudo in dB
title('Magnitude');
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
grid on;

subplot(2,1,2);
semilogx(freq_Hz_Notch, squeeze(phase_Notch)); % Fase in gradi
title('Phase');
xlabel('Frequency [Hz]');
ylabel('Phase [°]');
grid on;
sgtitle('Notch Filter (50 Hz, Bandwidth 5 Hz)')

%% notch 1.25 Hz
notch_freq = 1.25;  % Frequenza del disturbo (es. 50 Hz per il rumore della linea elettrica)
notch_bandwidth = 20;  % Larghezza di banda del filtro Notch (es. 5 Hz)


% Creazione del filtro Notch
[f_Notch, b_Notch] = iirnotch(notch_freq / (fs / 2), notch_bandwidth / (fs / 2));

% Convertiamo i coefficienti del filtro in oggetti di sistema per poter usare bode
sys_Notch = tf(f_Notch, b_Notch, 1/fs); % Filtro Notch


% Plot del Bode per il filtro Notch con asse x in Hz
figure('Name', 'Bode Plot - Notch Filter', 'NumberTitle', 'off');
[mag_Notch, phase_Notch, wout_Notch] = bode(sys_Notch);  % Ottenere dati Bode

% Convertire la frequenza da rad/s a Hz
freq_Hz_Notch = wout_Notch / (2 * pi);

% Grafico in Hz
subplot(2,1,1);
semilogx(freq_Hz_Notch, squeeze(20*log10(mag_Notch))); % Magnitudo in dB
title('Magnitude');
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
grid on;

subplot(2,1,2);
semilogx(freq_Hz_Notch, squeeze(phase_Notch)); % Fase in gradi
title('Phase');
xlabel('Frequency [Hz]');
ylabel('Phase [°]');
grid on;
sgtitle('Notch Filter (1.25 Hz, Bandwidth 20 Hz)')