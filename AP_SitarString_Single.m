% Author: Abhiroop Paul
% Date: 29th July 2026
%
% TITLE: FD Simulation of a Single Sitar String Colliding with a Distributed Bridge
%
% DESCRIPTION:
% FD simulation of a stiff string with internal losses (frequency dependent
% and independent) that collides with a distributed rigid bridge on one
% end.
%
% The output can be convolved with a recorded sitar IR.
% The user must input the midi value of the note they desire to generate
% (for example 60 for the note C4).
%
% DEPENDENCIES:
% The user must ensure the following files are in the folder path for the
% scheme to work:
%
% AP_midi2f.m
% AP_D2x.m
% AP_D4x.m
% AP_myfastconv.m
% AP_fadeOut.m
% IR_Contact_1.wav / IR_Condenser_1.wav
% AP_midi2note.m
%
% REFERENCE
% [1] Bilbao, S. (2009). Numerical sound synthesis: finite difference schemes and
% simulation in musical acoustics. John Wiley & Sons.
%
% =======================================================================
clear all; close all; clc;

name    =   "Bridge_Single_";
% ====================================
% GLOBAL INPUT PARAMETERS
% ====================================
% Sampling Related
Fs      =   48000;                          % Sample Rate (Hz) (48k RECOMMENDED)
Tf      =   10.0;                            % Simulation duration (seconds)
tFo     =   0.025;                          % Fadeout time (seconds)

% String System Related
midi    =   55;                             % Midi note value
Bxn     =   1e-7;                           % Inharmonicity Parameter

% Losses (set by trial and error to achieve empirical T60 values)
sig0    =   8.2e-1;                         % Frequency Independent loss param
sig1    =   9e-4;                           % Freq. dependent loss param

% Bridge Related
lB      =   0.06;                           % Fractional length of bridge (0.04 - 0.06)
aB      =   5e-9;                           % Coefficient of x^2 for parabolic bridge
kB      =   6e17;                           % Bridge Stiffness parameter
alpha   =   3.2;                            % alpha parameter for potential (>3.5 for stability)
tol     =   1e-18;                          % Tolerance for 0    
NRtol   =   1e-12;                           % Threshold for Newton Raphson
maxiter =   20;                             % NR Iteration Limit

% Boundary Condition
BC      =   1;                              % 0 = Clamped, 1 = Simply Supported

% Initial Condition
exc     =   1;                              % 0 = raised cos, 1 = triangle                                  
x_pluck =   0.25;                           % Normalised position of pluck
w_pluck =   5e-1;                           % Width of raised cosine
h_pluck =   5e-5;                           % Height of pluck for triangular IC
u0      =   1;                              % Scaling of initial displacement
v0      =   0;                              % Initial velocity

% Readout Related
rp      =   0.25;                           % Normalised position of readout

% Output Related
animate =   false;                          % Turn on to generate animation (recommend making Tf < 1sec)
velo    =   true;                           % Output reads velocity if true, disp. if false
bodyIR  =   1;                              % Body filter: 0 = no conv., 1 = contact mic IR, 2 = cond. mic IR
listen  =   true;                           % Turn on for sound output 
write   =   false;                          % Turn on to write to wav


% ====================================
% GLOBAL CALCULATED PARAMETERS
% ====================================
% Sampling Related
Nf      =   floor(Fs * Tf);                 % Duration in samples
T       =   1/Fs;                           % Time step (seconds)

% String System Related
f0      =   AP_midi2f(midi);                % Calculation of fundamental from midi
gamma   =   2 * f0;                         % Gamma Parameter
k       =   (gamma/pi)*sqrt(Bxn);           % Stiffness Parameter "k"

% Grid Related
h       =   sqrt((gamma*T)^2 + sqrt((((gamma*T)^4 + (16*(k*T)^2))))/2);
N       =   floor(1/h);                     % Number of grid points
h       =   1/N;                            % Grid spacing final

% Bridge Related
lN      =   ceil(N*lB);                    % Number of points in bridge
al1     =   alpha + 1;                      % alpha + 1 in phi calculation

% Scheme Related
lamhsq  =   (T * gamma)^2;                  % Lambda squared * h = 1 (for stability)
muhsq   =   (T * k)^2;                      % Mu squared * h

% Initial Condition
N_pluck =   floor(N*x_pluck);               % Grid index of pluck position


% ====================================
% SYSTEM MATRICES
% ====================================
D2x     =   AP_D2x(N,h,BC);             % Spatial 2nd order diffOp dep. on AP_D2x.m
D4x     =   AP_D4x(N,h,BC);             % Spatial 4th order diffOp dep. on AP_D4x.m

A   =   sparse((1+(sig0*T))*eye(N+1) - (sig1*T)*D2x);
B   =   sparse(2*eye(N+1) + lamhsq*D2x - muhsq*D4x);
C   =   sparse((1-(sig0*T))*eye(N+1) + (sig1*T)*D2x);

Q   =   A(1:lN,1:lN);


% ====================================
% OUTPUT INTERPOLATION (1st Order)
% ====================================
lr      =   1 + floor(N*rp);
rFrac   =   N*rp - floor(N*rp);
R       =   zeros(1,N+1);
R(lr)   =   N*(1-rFrac);
R(lr+1) =   N*rFrac;


% ====================================
% INITIAL CONDITIONS
% ====================================
if floor(mod(exc,2)) == 0  % If exc ID is even:
    % Raised cosine IC [1]
    xax     =   [1:N+1]'*h;
    ind     =   sign(max(-(xax-x_pluck-w_pluck/2).*(xax-x_pluck+w_pluck/2),0));
    ic      =   h_pluck*0.5*ind.*(1+cos(2*pi*(xax-x_pluck)/w_pluck));
else
    % Triangular Pluck IC
    up      =   [0:1:N_pluck]' * (h_pluck / N_pluck);
    down    =   flip([0:1:N - N_pluck - 1]') * (h_pluck / N_pluck);
    ic      =   [up;down];
end


% ====================================
% INITIALISING VECTORS
% ====================================
x       =   [0:N]'/N;                       % x vector for bridge

% Bridge Construction
b       =   zeros(N+1,1);
b(1:lN) =   -aB*(x(1:lN).^2);
b(lN+1:end)=   -(3.0*h_pluck);

u       =   zeros(N+1,1);                   % Empty vector of grid-size
u_temp  =   zeros(N+1,1);

et      =   zeros(N+1,1);
et_temp =   zeros(N+1,1);

Ffull   =   zeros(N+1,1);                   % Force on string from bridge

out     =   zeros(Nf,1);                    % Empty output vector


% ====================================
% SET INITIAL CONDITIONS
% ====================================
u2      =   u0*ic;                          % First Time Step 
u1      =   (u0+T*v0)*ic;                   % Second Time Step

et2     =   u2 - b;
et1     =   u1 - b;


% ====================================
% ANIMATION STORAGE
% ====================================
FPS = 20;

nFrames = floor(Nf/FPS);
frames = zeros(N+1,nFrames);

iFrame = 1;

% Timer start for runtime
tic
% ====================================
% MAIN LOOP (INCLUDES NR AT EACH STEP)
% ====================================
for n   =   3:Nf
    % Vectorised uncolliding update
    u_temp  = A \ (B*u1 - C*u2);

    % Temporary Penetration
    et_temp = u_temp - b;

    % eta calculation
    a       =   et1(1:lN);
    c       =   et_temp(1:lN) - et1(1:lN);

    r       =   zeros(lN,1);                % First guess for NR

    F       =   zeros(lN,1);                % Vector to hold force values
    delF    =   zeros(lN,1);                % Vector to hold force derivatives


    G       =   ones(lN,1);
    dG      =   inf(lN,1);

    iter    =   0;

    % -------------------------------
    % NEWTON RAPHSON
    % -------------------------------
    while norm(dG) > NRtol && iter < maxiter

        ps1 = max(a,0);
        ps2 = max(a+r,0);

        % Subsets of r where r = 0 (or small in practical case)
        zr = abs(r) < tol;
        nz = ~zr;

        % Initialising F(r) and F'(r) vectors
        F = zeros(lN,1);
        delF = zeros(lN,1);

        % Elements of F(r) & corresponding F'(r) where r = 0
        F(zr) = (kB*(ps1(zr).^alpha))/2;
        delF(zr) = 0.5*kB*alpha*(ps1(zr).^(alpha-1));

        % Elements of F(r) & corresponding F'(r) for all other values of r
        F(nz) = (kB/al1) .* ((ps2(nz).^al1)-(ps1(nz).^al1)) ./ r(nz);
        delF(nz) = (r(nz).*kB.*(ps2(nz).^alpha) - (kB/al1).*((ps2(nz).^al1)-(ps1(nz).^al1))) ./ (r(nz).^2);

        % Jacobian
        JF = sparse(diag(delF));
        JG = sparse(eye(lN) + T^2*(Q\JF));

        % Function to G(r) that needs to converge
        G = r + T^2*(Q\F) + c;

        dG = -JG\G;

        r = r + dG;

        iter = iter + 1;
    end

    % disp(iter);

    % Recalculation of F after NR
    ps1 = max(a,0);
    ps2 = max(a+r,0);

    zr = abs(r) < tol;                      % subset of r where ri = 0
    nz = ~zr;                               % Remaining r

    F(zr) = kB*(ps1(zr).^alpha);            % subset of r where ri = 0
    F(nz) = (kB/al1) .* ((ps2(nz).^al1)-(ps1(nz).^al1)) ./ r(nz);

    Ffull   =   [F;zeros(N+1-lN,1)];        % Full size force vector

    % Real State Calculation
    u = u_temp + T^2*(A\Ffull);             % String position
    v   =   (u - u1)/T;                     % String velocity

    % Current Penetration
    et = b - u;
    
    % Readout
    if velo == true
        out(n)      =   R*v;
    else
        out(n)      =   R*u;
    end

    % Animation
    if mod(n,FPS)==0
        frames(:,iFrame) = u;
        iFrame = iFrame + 1;
    end

    % Update
    u2  =   u1;     et2 =   et1;
    u1  =   u;      et1 =   et;

end

% ====================================
% ANIMATION
% ====================================
if animate == true
    figure
    for j = 1:size(frames,2)

        plot(x, b, "k", x,frames(:,j), "r")
        ylim([-2*h_pluck 2*h_pluck])
        xlim([0 1])

        xlabel("String Axis")
        ylabel("Displacement")

        drawnow
    end
end


% ====================================
% OUTPUT SECTION
% ====================================
out     =   out - mean(out);                    % Remove DC Offset
out     =   out / max(abs(out));                % Normalisation

% Convolving with body IR
if floor(mod(bodyIR,3)) == 1
    hIR     =   audioread("IR_Contact_1.wav");
    body    =   "ContactIR";  
elseif floor(mod(bodyIR,3)) == 2
    hIR     =   audioread("IR_Condenser_1.wav");
    body    =   "CondIR";
else
    body    =   "NoBody";
    hIR     =   [1;zeros(length(out),1)];
end
out     =   AP_myfastconv(hIR,out);

% FADEOUT
out     =   AP_fadeOut(out,tFo,Fs);
out     =   out/max(abs(out));       % Renormalise

% Playing Sound
if listen == true
    soundsc(out,Fs);
end

% Timer end for runtime
toc

% WRITE TO FILE
if write == true
    % Audiowrite
    note    =   AP_midi2note(midi(1));
    
    audiowrite(name+note+".wav",out,Fs);
end

% ==================================END====================================