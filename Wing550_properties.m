% Drone with 550mm Wing

%% Mass Properties 
m = 1.88;   % kg
Jn = 0.074741;  % kgm^2
Je = 0.074757;  % kgm^2
Jd =  0.144779;  % kgm^2

%% Wing Characteristics (From xflr5)
Iw = 0.000403;    % Inertia for Wing at rotation axis
l     = 0.47;         % Arm Length (m)
b     = l * 2;        % Wingspan (m)
AR    = 7.932;        % Aspect ratio
tr    = 1.37;         % Taper ratio
S     = (pi*l^2)/3; % 0.116/2;      % Wing Area (m^2)
y_MAC = (b / 6) * ((1 + 2 * tr) / (1 + tr)); % Wing location of MAC
MAC   = 0.124;        % Mean Aerodynamic Chord