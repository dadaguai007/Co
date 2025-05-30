function [Ex, Ey] = pbs(E, theta)
    % Polarization beam splitter (PBS).
    % 
    % Parameters
    % ----------
    % E : (N,2) Input pol. multiplexed optical field.
    % theta : Rotation angle of input field in radians. The default is 0.
    % 
    % Returns
    % -------
    % Ex : (N,1) Ex output single pol. field.
    % Ey : (N,1) Ey output single pol. field.
    
    if nargin < 2
        theta = 0;  % Default rotation angle
    end
    
    % Check and reshape E if necessary
    if isrow(E)
        E=E.';
    end
    % N×2
    if size(E, 2) == 1
        E = [E, zeros(size(E))];  % Assume Y component is 0
    end
    
    % Create the rotation matrix
    rot = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    
    % Rotate the input optical field
    E = E * rot;
    
    % Extract X and Y polarization components
    Ex = E(:, 1);
    Ey = E(:, 2);
end
