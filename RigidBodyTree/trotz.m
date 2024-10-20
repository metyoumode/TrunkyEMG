function R = trotz(theta)
    % Matrice di rotazione attorno all'asse Z
    R = [cos(theta) -sin(theta) 0 0;
         sin(theta) cos(theta) 0 0;
         0 0 1 0;
         0 0 0 1];
end