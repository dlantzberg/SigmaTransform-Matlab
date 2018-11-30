function [ t ] = freq2time( w , centertime )
%FREQ2TIME gives time axis from freq axis
%   USAGE: [ w ] = time2freq( t )
%
%	INPUT:
%       w         : freq
%       centertime: optional parameter, indicating
%                   that the time-axis should be symmetric
%	OUTPUT:
%       t         : time
%
%	AUTHOR:	D Lantzberg, Nov. 2015

	N	= length(w);
	dt	= -1/w(1)/2;

	t	= linspace(-N/2*dt , N/2*dt , N+1*(~mod(N,2)));
	t	= t(1:end-~mod(N,2));

    if~exist('centertime','var')
        t = t - t(1);
    end;
end
