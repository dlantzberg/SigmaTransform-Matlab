function [ nf ] = norm1( f , p )
%NORM1 "p"-normalizes f
%   USAGE: [ nf ] = norm1( f , p )
%	INPUT:
%		f		: vector
%		p		: p-norm
%	OUTPUT:
%		nf		: normalized vector
%
%	AUTHOR:	D Lantzberg, Nov. 2014

    if~exist('p')
        p = 2;
    end;
    nf = f / norm(f,p);
end
