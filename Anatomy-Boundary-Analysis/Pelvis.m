 classdef Pelvis
    % Pelvis class for the pelvis attributes
    
    % Developed by C.Micheler,
    % Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
    % Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich
    
    
    properties 
        import = Import(0);     % Import: loaded stl data  
        boundaries = Boundary;  % Object array for class boundary: Boundary body around pelvis
    end
    
    methods
        %% Constructer: generate object
        function obj = Pelvis()
            disp('class pelvis initialized')
        end  
    end
    
 end