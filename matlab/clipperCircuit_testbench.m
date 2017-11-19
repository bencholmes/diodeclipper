clear;

%% Simulation settings

fs = 48e3;          % Sampling frequency
T  = 1/fs;          % Sampling period
dur = 0.1;          % Signal duration
Ns = dur*fs;        % Number of samples in the signal
tv = (0:Ns-1)./fs;  % Time vector

%% Simulate

dcModel = clipperCircuit(fs);

% Input signal
u = 3*sin(2*pi*100*tv);

y = dcModel.simulate(u);

%% Plot

figure(1);
clf;
plot(tv,u,tv,y);
xlabel('Time (s)');
ylabel('Voltage');
legend('Input','Output');