%% 3.1A

% 1. Read the audio file
[fileData, fs] = audioread('oddity.wav');

% 2. Verify Sampling Rate
if fs == 44100
    disp('The sampling rate is 44.1 kHz');
else
    disp(['The sampling rate is ', num2str(fs/1000), ' kHz']);
end

% 3. Extract Left Track
m_t = fileData(:, 1);

% 4. Plot the Signal
time = (1:length(m_t))/fs; % Create a time vector
figure;
plot(time, m_t);
xlabel('Time (s)');
ylabel('Amplitude');
title('3.1 Part A');

% 5. Play the Audio Clip
soundsc(m_t, fs);

% 6. Compute Total Energy
energy = sum(abs(m_t).^2) / fs;
disp(['Total energy in the track: ', num2str(energy)]);


%% 3.1B


A_c = 10; % Carrier gain
f_s = 44100; % Original sampling rate
f_c = 1310e3; % Carrier frequency in Hz
f_h = 10e6; % Simulation sampling rate

% Call the upconvert function
 upsampleFactor = round(f_h / f_s);
    MessageSignalRF = zeros(upsampleFactor * length(m_t), 1);
    MessageSignalRF(1:upsampleFactor:end) = m_t;
    Th=1/f_h;
    t = (1:length( MessageSignalRF))/f_h;
    c = A_c * cos(2 * pi * f_c * t);
    s_t =  MessageSignalRF.* c;
time2 = 0:(length(MessageSignalRF)-1) *Th;
figure;
plot(t,m_t);
xlabel('Time (s)');
ylabel('Amplitude');
title('m(t) sampled at fh');
energy_m = sum(abs(MessageSignalRF).^2) / f_s;
disp(['Energy of upsampled m(t): ', num2str(energy_m)]);

figure;
T = length(s_t)/f_h;
t = 0:1/f_h:T-1/f_h;
plot(t, s_t);
title(' Signal s(t)');
xlabel('Time (s)');
ylabel('Amplitude');

energy_s =  sum(abs(s_t).^2) / f_s;
disp(['Energy of s(t): ', num2str(energy_s)]);


%% 3.1C

% Assuming s_t is your modulated signal
s_t = s_t(:,1);

% Define the values of sigma^2
sigma_squared_values = [0.01, 0.1, 1.0];

% Create a new figure for the subplots
figure;

for i = 1:length(sigma_squared_values)
    sigma_squared = sigma_squared_values(i);
    
    % Generate Gaussian noise for the given sigma^2
    sigma = sqrt(sigma_squared);
    n = length(s_t);
    w_t = sigma * randn(1, n);
    w_t = w_t';

    % Compute the received signal y(t)
    y_t = w_t + s_t;

    % Create subplots
    subplot(length(sigma_squared_values), 2, 2 * (i - 1) + 1);
    plot(t, s_t);
    title(['Modulated Signal s(t) for \sigma^2 = ' num2str(sigma_squared)]);
    xlabel('Time (s)');
    ylabel('Amplitude');
    
    subplot(length(sigma_squared_values), 2, 2 * i);
    plot(t, y_t);
    title(['Received Signal y(t) for \sigma^2 = ' num2str(sigma_squared)]);
    xlabel('Time (s)');
    ylabel('Amplitude');
end






