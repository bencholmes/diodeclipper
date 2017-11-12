%%% Diode clipper - Virtual analogue model of a diode clipper
% Written by Ben Holmes
% 2017/11/12
% GPL

% This script is written to clearly identify the different parts required
% to simulate a circuit.

% The diode clipper is a 2 node distortion circuit, featuring a resistor
% between nodes 1 and 2, and a capacitor and two antiparallel diodes
% between 2 and ground.

% This model is created assuming that the diodes are exactly equivalent.
% The trapezoidal rule is used to discretise the system. Newton interation
% is used to solve the transcendental equation.

% Continuous time representation:
% (Vi-Vo)/R - C dv/dt - 2*Is*sinh(V/(NVt) = 0

% Discrete time representation
% (Vi(n)-Vo(n))/R - C T/2 (Vo(n-1) + Vo(n)) - 2*Is*sinh(Vo(n)/(NVt) = 0

%% Simulation settings

fs = 48e3;          % Sampling frequency
T  = 1/fs;          % Sampling period
dur = 0.01;         % Signal duration
Ns = dur*fs;        % Number of samples in the signal
tv = (0:Ns-1)./fs;  % Time vector

%% Model parameters

R = 100;        % Resistance
C = 0.1e-6;     % Capacitance
Is = 1e-8;      % Saturation current
N = 1.2;        % Ideality factor
Vt = 25e-3;     % Thermal voltage

%% Input/output signals

% Input signal
u = 3*sin(2*pi*500*tv);

% Output defined as vector of zeros to preallocate memory
y = zeros(size(u));

%% Processing loop
for nn=1:Ns
    
    % If at the start of the loop give the memory a zero value
    % This is valid as the circuit has no DC offset
    if nn==1, yn1=0; else, yn1 = y(nn-1); end
    
    % This value sets the initial iterate, i.e. the first guess to the
    % solution. the previous solution is usually pretty good.
    yc = yn1;       % Current value of the output (being solved for)
    
    res = 10;       % set residual above limit of while loop
    iters = 0;      % set the number of iterations to 0
    
    % Newton solver
    % This loop will exit when either the step taken is very small (i.e.
    % the newton algorithm has converged) or too many iterations have been
    % performed.
    while res>1e-10 && iters < 100
        
        % Calculate the function value
        f = (u(nn)-yc)/R - (C*T/2)*(yc+yn1) - 2*Is*sinh(yc/(N*Vt));
        
        % Calculate it's derivative
        j = -1/R - (C*T/2) - (2*Is/(N*Vt))*cosh(yc/(N*Vt));
        
        % Newton step is defined by the function divided by the derivative.
        step = -f/j;
        
        % Only magnitude is usefuly for defining a residual
        res = abs(step);
        
        % Take a Newton step towards the solution (potentially)
        yc = yc + step;
        
        % Increment iterations
        iters = iters + 1;
    end
    
    % Place solution into output vector
    y(nn) = yc;
end

%% Plot

figure(1);
clf;
plot(tv,u,tv,y);
xlabel('Time (s)');
ylabel('Voltage');
legend('Input','Output');

%% Play

% Uncomment to hear the output
% soundsc(y,fs)