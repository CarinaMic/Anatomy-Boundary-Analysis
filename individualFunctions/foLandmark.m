% Calculation of the landmark foramen obturatum (fo) 

% Input 'fo_p', 'fo_i', 'fo_a', 'fo_ip': Coordinates of four points within foramen obturatum (fo)
%       The four points describe the foramen obturatum by approximating it as an ellipse and ...
%       describing the major and minor axes accordingly.

% Output fo: Coordinates of the landmark fo (centre point of fo)

% Developed by C.Micheler, 
% Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich

function fo = foLandmark(fo_p, fo_i, fo_a, fo_ip)
    
    % Calculation of landmark fo (centre of foramen obturatum):
    % half distance between two skewed straight lines
    
    % Straight line equation with two points (P,Q): x = P + lambda*(P-Q)
    % Line fo-p to fo-i (main axis fo elipse): x = fo_p + lambda*u
    u_main = fo_p - fo_i;
    % Line fo-a to fo-ip (minor axis fo elipse): x = fo_a + phi*u
    u_minor = fo_a - fo_ip;
    
    % Vector perpendicular to both direction vectors of the straight line:
    % cross product of u_main and u_minor
    n = cross(u_main, u_minor);
    
    % Auxiliary plane: x = fo_p + lambda*u_main + omega*n
    % Plump bob point: equalise auxiliary plane and line fo-a to fo-ip
    % fo_p + lambda*u_main + omega*n = fo_a + phi*u_minor
    % System of equations: lampda*u_main + omega*n - phi*u_minor = fo_a-fo_p
    % gr = inv(A)*(fo_a-fo_p); gr = [lambda; omega; phi]
    A = [u_main', n', (-u_minor)'];
    b = (fo_a - fo_p)';
    gr = A\b; % system of equations
    
    % Plump bob points with straight line equations
    S1 = fo_p + gr(1)*u_main;
    S2 = fo_a + gr(3)*u_minor;
    
    % Centre point of vector S1-S2 (half distance)
    fo = 0.5 * (S1 + S2);

end

