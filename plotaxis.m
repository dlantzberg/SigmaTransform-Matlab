function [  ] = plotaxis( n,cola,LineWidth )
%PLOTAXIS plots an axis through the origin
%   USAGE: [  ] = plotaxis( n,cola )
%	INPUT:
%		n		: distance from axis to boundary in percent
%       cola    : colour and line specs
%	OUTPUT:
%
%	AUTHOR:	D Lantzberg, Jan. 2014
    
	if~exist('LineWidth')
		LineWidth = 1;
	end;

	if~exist('cola')
		cola = 'k';
	end;
	
	if~exist('n')
		n = 5;
	end;
    N = 2^2;
    reax(n,n);
    ax = axis;
    hor = linspace( ax(1), ax(2) , N);
    ver = linspace( ax(3), ax(4) , N);
    zero = zeros(1,N);
    hold on,
    plot( hor,zero, cola, zero, ver, cola );
    hold off;
    reax(n,n);
end

