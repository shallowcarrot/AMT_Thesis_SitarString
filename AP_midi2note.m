function name = AP_midi2note(midi)
% Author: Abhiroop Paul
% Date: 2nd Aug 2026
%
% Title: Midi to Note Converter
% 
% Description:
% Function to convert midi note value into note name.
% =======================================================================

oct     =   num2str(floor(midi/12) - 1);    % Keyboard octave the note is in
note    =   mod(midi,12);                   % Note number in chromatic scale (0 indexing)

% Note Name assignment
if note == 0;
    name = 'C';
elseif note == 1
    name = 'C#';
elseif note == 2
    name = 'D';
elseif note == 3
    name = 'D#';
elseif note == 4
    name = 'E';
elseif note == 5
    name = 'F';
elseif note == 6
    name = 'F#';
elseif note == 7
    name = 'G';
elseif note == 8
    name = 'G#';
elseif note == 9
    name = 'A';
elseif note == 10
    name = 'A#';
elseif note == 11
    name = 'B';
end

name    =   strcat(name,oct);

end