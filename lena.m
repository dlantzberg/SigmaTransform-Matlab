function [L,Fs2D] = lena()
%LENA returns LEN(N)A in 128 x 128 and greyscale
%   USAGE: [L,Fs2D] = lena()
%	INPUT:
%
%	OUTPUT:
%		L       : lena "signal"
%		Fs2D    : assuming spatial domain [0,1)x[0,1),
%                 gives the "spatial frequencies" 128 x 128
%
%	AUTHOR:	Daniel Lantzberg, Nov. 2017

    % load image
    L = rgb2gray( imread('lena.png') );
    L = imresize(L,1/4);
    Fs2D = size(L);
end
