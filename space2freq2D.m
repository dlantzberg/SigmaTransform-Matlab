function [ Wx, Wy ] = space2freq2D( Tx , Ty )
%SPACE2FREQ2D gives Frequency-Plane for spatial plane
%   USAGE: [ Wx, Wy ] = space2freq2D( Tx , Ty )
%   INPUT:
%       Tx, Ty : Time-Domain-Matrices, as given by "meshgrid"
%   OUTPUT:
%       Wx, Wy : Corresponding Frequency-Domain-Matrices
%
%	AUTHOR:	D Lantzberg, Nov. 2016


	[ height , width ] = size(Tx);

	dth = Tx(1,2)-Tx(1,1);
	dtv = Ty(2,1)-Ty(1,1);

	fh_max = 1/dth/2;
	fv_max = 1/dtv/2;

	dwh = 2*fh_max/width;
	dwv = 2*fv_max/height;

	[ Wx , Wy ] = meshgrid( 0 : dwh : fh_max , 0 : dwv : fv_max ) ;

	Wx = [ -fliplr(Wx(:,2:end)) , Wx  ];
	Wx = [  Wx ;  flipud(Wx(2:end,:)) ];

	Wy = [ -flipud(Wy(2:end,:)) ; Wy ];
	Wy = [  fliplr(Wy(:,2:end)) , Wy ];

	% delete highest freq, if even
	if( ~mod(height,2) )
		Wy=Wy(1:end-1,:);
		Wx=Wx(1:end-1,:);
	end;
	if( ~mod(width,2) )
		Wx=Wx(:,1:end-1);
		Wy=Wy(:,1:end-1);
	end;
end
