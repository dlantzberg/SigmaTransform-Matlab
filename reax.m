function [  ] = reax( hor, ver )
%REAX resizes axes symmetrically
%   USAGE: [  ] = reax( hor, ver , depth )
%	INPUT:
%		hor		: x resize
%		ver		: y resize
%	OUTPUT:
%
%	AUTHOR:	D Lantzberg, Jan. 2014
    
    if~exist('hor')
        hor = 10;
    end;
    if~exist('ver')
        ver = 10;
    end;
    hor = hor/100/2;
    ver = ver/100/2;
    ax = axis; 
    deltah = ax(2) - ax(1);
    deltav = ax(4) - ax(3);
    ax(1) = ax(1) - deltah*hor; 
    ax(2) = ax(2) + deltah*hor; 
    ax(3) = ax(3) - deltav*ver;
    ax(4) = ax(4) + deltav*ver; 
    axis(ax);

end

