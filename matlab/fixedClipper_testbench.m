clear;

%% Simulation settings

fs = 48e3;          % Sampling frequency
T  = 1/fs;          % Sampling period
dur = 0.02;          % Signal duration
Ns = dur*fs;        % Number of samples in the signal
tv = (0:Ns-1)./fs;  % Time vector

%% Simulate

varModel = clipperCircuit(fs);
fixModel = fixedClipperCircuit(fs);

% Input signal
u = 10*sin(2*pi*1e3*tv);

tic
yV = varModel.simulate(u);
varTime = toc;

tic
yF = fixModel.simulate(u);
fixTime = toc;

fprintf('Variable time: %g, Fixed time: %g\n\n',varTime,fixTime);

%% Plot

figure(1);
clf;
subplot(211);
plot(tv,yV,tv,yF);
xlabel('Time (s)');
ylabel('Output Voltage (V)');
legend('Variable','Fixed');
subplot(212);
plot(tv,yV-yF);
xlabel('Time (s)');
ylabel('Error');