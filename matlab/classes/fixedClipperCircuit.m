classdef fixedClipperCircuit < clipperCircuit
    %FixedClipperCircuit: A circuit focussed representation of the diode
    %clipper, using a fixed number of iterations to remove branch
    %operations.
    
    methods
        function m = fixedClipperCircuit(fs)
            m = m@clipperCircuit(fs);
        end
        
        function y = simulate(m,u)
            
            if (size(u,1) > 1) && (size(u,2) > 1)
                error('Input must be 1D.');
            end
            
            Ns = length(u);
            y = zeros(1,Ns);
            
            for nn=1:Ns
                
                % If at the start of the loop give the memory a zero value
                % This is valid as the circuit has no DC offset
                if nn==1, yn1=0; else, yn1 = y(nn-1); end
                
                % This value sets the initial iterate, i.e. the first guess to the
                % solution. the previous solution is usually pretty good.
                yc = yn1;       % Current value of the output (being solved for)
                
                % Newton solver
                % This loop will exit when either the step taken is very small (i.e.
                % the newton algorithm has converged) or too many iterations have been
                % performed.
                for mm=1:7
                    
                    % Calculate the function value
                    f = (u(nn)-yc)/m.R - ((2*m.C*m.fs)*yc) + m.x - 2*m.Is*sinh(yc/(m.N*m.Vt));
                    
                    % Calculate it's derivative
                    j = -1/m.R - (2*m.C*m.fs) - (2*m.Is/(m.N*m.Vt))*cosh(yc/(m.N*m.Vt));
                    
                    % Newton step is defined by the function divided by the derivative.
                    step = -f/j;
                    
                    % Take a Newton step towards the solution (potentially)
                    yc = yc + step;
                end
                
                % Place solution into output vector
                y(nn) = yc;
                
                % Calculate state
                m.x = (4*m.C*m.fs)*yc - m.x;
            end
        end
    end
    
end

