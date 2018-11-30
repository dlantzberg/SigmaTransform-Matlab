function [y,Fs] = bat()
%BAT returns bat signal
%   USAGE: [y,Fs] = bat()
%	INPUT:
%
%	OUTPUT:
%		y	: Bat signal
%		Fs	: sampling frequency
%
%	AUTHOR:	Daniel Lantzberg, Nov. 2017
%
%   Acknowledgement: The author wishes to thank Curtis Condon, Ken White, and Al Feng
%   of the Beckman Institute of the University of Illinois for the bat data and
%   for permission to use it in this example
%   http://dsp.rice.edu:80/software/bat-echolocation-chirp

    % load signal
    y   = load('-ascii','bat.asc');
    % set Fs
    Fs  = 143000;

end
