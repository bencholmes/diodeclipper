clear;

%% Simulation settings

fs = 48e3;          % Sampling frequency
T  = 1/fs;          % Sampling period
dur = 2;          % Signal duration
Ns = dur*fs;        % Number of samples in the signal
tv = (0:Ns-1)./fs;  % Time vector

%% Simulate

effect = clipperEffect(fs);

parameters = [0.8*ones(1,Ns);...    % Distortion
              0.5*ones(1,Ns);...    % Brightness
              linspace(0,1,Ns)];% Warmth

% Input signal
u = 2*sin(2*pi*1e2*tv);

tic
y = effect.simulate(u,parameters);
simTime = toc;

fprintf('Simulation time: %g\n\n',simTime);

%% Plot

figure(1);
clf;
plot(tv,y);
xlabel('Time (s)');
ylabel('Output Voltage (V)');
legend('Variable','Fixed');

%% Play

sound(y, fs)