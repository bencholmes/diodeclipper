classdef clipperEffect < handle
    %CLIPPEREFFECT A musical effect using the diode clipper model
    
    properties(SetAccess = private)
        fs = 48e3;
        dc = clipperCircuit(48e3);
        
        warmth = 0.5
        offset = 0;
        Vforward = 0.7
    end
    
    methods
        function m = clipperEffect(fs)
            m.fs = fs;
            m.dc = clipperCircuit(fs);
            
            m.setDistortion(0.5);
            m.setBrightness(0.5);
            m.setWarmth(0.5);
        end
        
        function setDistortion(m, distortion)
            m.dc.N = 2-((1.25*80^distortion - 1.25)*1.9e-2);
            m.Vforward = m.dc.N*m.dc.Vt*log((1e-3/(2*m.dc.Is) + 1));
            
            m.setWarmth(m.warmth);
        end
        
        function setBrightness(m, brightness)
            fc = 10e2*(1.25*80^brightness - 1.25);
            m.dc.C = 1/(2*pi*m.dc.R*fc);
        end
        
        function setWarmth(m, warmth)
            m.offset = warmth*(1-m.Vforward);
        end
        
        
        function y = simulate(m, u, parameters)
            % Parameters: distortion, brightness, warmth
            if nargin < 2
                error('Not enough input arguments.');
            end
            
            if size(u,1) > 1 && size(u,2) > 1
                error('Input signal must be 1D.');
            end
            
            Ns = length(u);
            
            if nargin == 3
                if size(parameters,1) ~= 3 || size(parameters,2) ~= Ns
                    error('Parameters must be 3xNs.');
                end
            elseif nargin ==2
                parameters = 0;
            else
                error('Too many input arguments.');
            end
            
            y = zeros(1,Ns);
            
            if size(parameters,1) == 1
                u = u + m.offset;
                y = m.dc.simulate(u);
            else
                for nn=1:Ns
                    m.setDistortion(parameters(1,nn));
                    m.setBrightness(parameters(2,nn));
                    m.setWarmth(parameters(3,nn));
                    
                    u(nn) = u(nn) + m.offset;
                    
                    y(nn) = m.dc.simulate(u(nn));
                end
                y = y./m.Vforward;
            end
            
        end
    end
    
end

