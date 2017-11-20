classdef clipperCircuit < handle
    %ClipperCircuit: A circuit focussed representation of the diode
    %clipper.
    
    % Simulation model properites
    properties
        fs  = 48e3;      % Sampling frequency
        x   = 0;         % Model state
        yn1 = 0;         % Last output (for iterative method)
    end
    
    % Circuit model properties
    properties
        R = 2200;       % Resistance
        C = 0.1e-6;     % Capacitance
        Is = 1e-10;     % Saturation current
        N = 1.2;        % Ideality factor
        Vt = 25e-3;     % Thermal voltage
    end
    
    methods
        function m = clipperCircuit(fs)
            m.fs = fs;
        end
        
        function y = simulate(m,u)
            
            if (size(u,1) > 1) && (size(u,2) > 1)
                error('Input must be 1D.');
            end
            
            Ns = length(u);
            y = zeros(1,Ns);
            
            for nn=1:Ns                
                % This value sets the initial iterate, i.e. the first guess to the
                % solution. the previous solution is usually pretty good.
                yc = m.yn1;       % Current value of the output (being solved for)
                
                res = 10;       % set residual above limit of while loop
                iters = 0;      % set the number of iterations to 0
                
                % Newton solver
                % This loop will exit when either the step taken is very small (i.e.
                % the newton algorithm has converged) or too many iterations have been
                % performed.
                while res>1e-10 && iters < 100
                    
                    % Calculate the function value
                    f = (u(nn)-yc)/m.R - ((2*m.C*m.fs)*yc) + m.x - 2*m.Is*sinh(yc/(m.N*m.Vt));
                    
                    % Calculate it's derivative
                    j = -1/m.R - (2*m.C*m.fs) - (2*m.Is/(m.N*m.Vt))*cosh(yc/(m.N*m.Vt));
                    
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
                
                % Update last solution
                m.yn1 = yc;
                
                % Calculate state
                m.x = (4*m.C*m.fs)*yc - m.x;
            end
        end
    end
    
end

