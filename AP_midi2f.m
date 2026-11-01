function f = AP_midi2f(midi)
% Author: Abhiroop Paul
% Date: 2nd Aug 2026
%
% Title: Midi to Frequency Converter
% 
% Description:
% Function to convert midi note value into frequency. Input needs only one
% argument, which is the midi note value.
% =======================================================================
Afreq   =   440;                % Centre of tuning convention (eg. A440, A432 etc)

oct     =   floor(midi/12) - 5; % Keyboard octave the note is in
note    =   mod(midi,12);       % Note number in chromatic scale (0 indexing)

% Frequency calc according to Afreq
f       =   Afreq*(2^oct)*(2^((note - 9)/12));

end

