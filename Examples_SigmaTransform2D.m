%   Some usage examples for SigmaTransform2D()
%
%   AUTHOR: D Lantzberg, 2017 - 2018

% is this octave or matlab?
if(exist('OCTAVE_VERSION', 'builtin'))
    disp('Octave detected, trying to load image package ..'); 
    disp('.. if it is not installed, you may do so with: ');
    disp('    pkg install image -forge');
    pkg load image;
end;

close all;

%% 2D-STFT, Lena 128x128, rectangular window
sigmaX = @(x,y) x;
sigmaY = @(x,y) y;

%actX     = @(x,y,xp,yp) x - xp;
%actY     = @(x,y,xp,yp) y - yp;

% load image of lena
[Lena,Fs2D]    = lena();

% define time- and frequency planes
t = linspace( 0, 1 , 128+1 ); t=t(1:end-1);
[X,Y] = meshgrid(t);
[WX,WY] = space2freq2D( X , Y );
WX = fftshift(WX); WY = fftshift(WY);

% define rectangular window (of "width" 5 x 5 )
win = @(X,Y)    ( X >= -2.5 ) .* ( X < 2.5 ) ...
            .*  ( Y >= -2.5 ) .* ( Y < 2.5 );

% make channels
[Xstep,Ystep] = meshgrid( -sigmaX(Fs2D(2)/2,Fs2D(1)/2) : 5 : sigmaX(Fs2D(2)/2,Fs2D(1)/2) , ...
                          -sigmaY(Fs2D(2)/2,Fs2D(1)/2) : 5 : sigmaY(Fs2D(2)/2,Fs2D(1)/2) );

% the transform
Lena2DSTFT = SigmaTransform2D(  ...
    Lena , ...                  % the image
    win  , ...                  % the window
    Xstep(:), Ystep(:) , ...    % the channels
    sigmaX ,  sigmaY , ...      % the diffeomorphism (x,y) |-> (sigmaX(x,y),sigmaY(x,y))
    WX,  WY ...                 % the Fourier-domain
);


%%% Alternative %%%

% define some function handle for a "named" transform like the "STFT2D" ...
STFT2D = @( signal , Xchannels, Ychannels , Wx, Wy ) ...
    SigmaTransform2D( signal , win , Xchannels, Ychannels , ...
                      @(x,y) x, @(x,y) y, Wx , Wy );

% ... and apply it
Lena2DSTFT = STFT2D(        ...
    Lena ,                  ... % the image
    Xstep(:), Ystep(:) ,    ... % the channels
    WX,  WY                 ... % the Fourier-domain
);

%%%%%%%%%%%%%%%%%%%

% Save Reconstruction
recLena2DSTFT = Lena2DSTFT.reconstruct('dual');

% Plot emerging reconstruction
Lena2DSTFT.emergingRec( 0.001 , 2 );
pause;

% plot full Reconstruction
imagesc( real(recLena2DSTFT )), axis ij; grid off, axis off;
title('Reconstruction, using dual frame');
pause;

%% 2D Wavelet, Lena 128x128, rectangular window
sigmaX = @(x,y) log2(abs(x)+eps);
sigmaY = @(x,y) log2(abs(y)+eps);

%actX     = @(x,y,xp,yp)    x - xp;
%actY     = @(x,y,xp,yp) y - yp;

% load image of lena
[Lena,~]    = lena();

% define time- and frequency planes
t = linspace( 0, 1 , 128+1 ); t=t(1:end-1);
[X,Y] = meshgrid(t);
[WX,WY] = space2freq2D( X , Y );
WX = fftshift(WX); WY = fftshift(WY);

% define rectangular window ( of "width" 1x1 )
win = @(X,Y)   (X <= 1 ) .* ( X > 0 ) ...
            .* (Y <= 1 ) .* ( Y > 0 );

% make channels
[Xstep,Ystep] = meshgrid(   linspace(log2(1) , log2(32) , 6 ),...
                            linspace(log2(1) , log2(32) , 6 ) );

% the transform
Lena2DWT = SigmaTransform2D(  ...
    Lena , ...                  % the image
    win  , ...                  % the window
    Xstep(:), Ystep(:) , ...    % the channels
    sigmaX ,  sigmaY , ...      % the diffeomorphism (x,y) |-> (sigmaX(x,y),sigmaY(x,y))
    WX,  WY ...                 % the Fourier-domain
);

% Save Reconstruction
recLena2DWT = Lena2DWT.reconstruct('resid');

% Plot emerging reconstruction
Lena2DWT.emergingRec( .1 , 1 );
pause;

% plot full Reconstruction
imagesc( real(recLena2DWT )), axis ij; grid off, axis off;
title('Reconstruction, using dual frame and residuum');
pause;

%% (nonParabolic) SIM(2)-Let, Lena 128x128, rectangular window
sigmaX   = @(x,y) log2(x.^2+y.^2 + eps)/2;
sigmaY   = @(x,y) atan(y./x);

%actX     = @(x,y,xp,yp)  x - xp;
%actY     = @(x,y,xp,yp)  y - yp;

% load image of lena
[Lena,~]    = lena();

% define time- and frequency planes
t = linspace( 0, 1 , 128+1 ); t=t(1:end-1);
[X,Y] = meshgrid(t);
[WX,WY] = space2freq2D( X , Y );
WX = fftshift(WX); WY = fftshift(WY);

% define rectangular window (of "width" 1 x (pi/16) )
win = @(X,Y)   (Y < 2*pi/32 ) .* ( Y >= 0) ...
            .* ( X > 0 ) .* ( X <= 1);

% make channels
[Xstep,Ystep] = meshgrid(   log2(1) : 1 : log2(64) , ...
                            -pi/2 : 2*pi/32 : pi/2 );

% the transform
LenaSIM2T = SigmaTransform2D(  ...
    Lena , ...                  % the image
    win  , ...                  % the window
    Xstep(:), Ystep(:) , ...    % the channels
    sigmaX ,  sigmaY , ...      % the diffeomorphism (x,y) |-> (sigmaX(x,y),sigmaY(x,y))
    WX,  WY ...                 % the Fourier-domain
);

%%% Alternative %%%

% define some function handle for a "named" transform like the "SIM2Transform" ...
SIM2Transform = @( signal , Xchannels, Ychannels , Wx, Wy ) ...
    SigmaTransform2D( signal , win , Xchannels, Ychannels , ...
                      sigmaX, sigmaY, Wx , Wy );

% ... and apply it
LenaSIM2T = SIM2Transform(        ...
    Lena ,                  ... % the image
    Xstep(:), Ystep(:) ,    ... % the channels
    WX,  WY                 ... % the Fourier-domain
);

%%%%%%%%%%%%%%%%%%%

% Save Reconstruction
recLenaSIM2T = LenaSIM2T.reconstruct('resid');

% Plot reconstruction
LenaSIM2T.emergingRec( 0.1 , 1 );
pause;

% plot full Reconstruction
imagesc( real(recLenaSIM2T )), axis ij; grid off, axis off;
title('Reconstruction, using dual frame and residuum');
pause;

%% 2D-NonParabolic SHEAR, Lena 128x128, rectangular window
sigmaX  = @(x,y) log2(abs(x)+eps);
sigmaY  = @(x,y) y./x;

%actX    = @(x,y,xp,yp) x - xp;
%actY    = @(x,y,xp,yp) y - yp ;

% load image of lena
[Lena,~]    = lena();

% define time- and frequency planes
t = linspace( 0, 1 , 128+1 ); t=t(1:end-1);
[X,Y] = meshgrid(t);
[WX,WY] = space2freq2D( X , Y );
WX = fftshift(WX); WY = fftshift(WY);

% define rectangular window  (of "width" 2 x 1 )
win = @(X,Y)   ( Y <= 1/2) .* ( Y > -1/2) ...
            .* ( X >   0 ) .* ( X <=  2 );

% make channels
[Xstep,Ystep] = meshgrid( (0: 2 : 8) , -4 : 1 : 4 );

% the Transform
LenaNPShear = SigmaTransform2D(  ...
    Lena , ...                  % the image
    win  , ...                  % the window
    Xstep(:), Ystep(:) , ...    % the channels
    sigmaX ,  sigmaY , ...      % the diffeomorphism (x,y) |-> (sigmaX(x,y),sigmaY(x,y))
    WX,  WY ...                 % the Fourier-domain
);

% Save Reconstruction
recLenaNPShear = LenaNPShear.reconstruct('resid');

% plot reconstruction
LenaNPShear.emergingRec( 0.1 , 1 );
pause;

% plot full Reconstruction
imagesc( real(recLenaNPShear )), axis ij; grid off, axis off;
title('Reconstruction, using dual frame and residuum');
pause;

%% Curvelet / (Parabolic SIM(2)-let)  Lena 128x128, rectangular window
sigmaX = @(x,y) log2(x.^2+y.^2+eps)/2;
sigmaY = @(x,y) atan(y./x);
%detSIGMA = @(x,y) exp(-2.*x);

actX    = @(x,y,xp,yp)                x - xp;
actY    = @(x,y,xp,yp) 2.^(-xp*.5) .*( y - yp );

% load image of lena
[Lena,~]    = lena();

% define time- and frequency planes
t = linspace( 0, 1 , 128+1 ); t=t(1:end-1);
[X,Y] = meshgrid(t);
[WX,WY] = space2freq2D( X , Y );
WX = fftshift(WX); WY = fftshift(WY);

% define rectangular window (of "width" 2 x (pi/32) )
win = @(X,Y)   ( Y >= 0 ) .* ( Y <  2*pi/64 ) ...
            .* ( X >  0 ) .* ( X <= 2 );

% make channels
[Xstep,Ystep] = meshgrid( ( 0 : 2 : 6 ) , pi/2*(-1:1/16:1) );
% adjust step-width to parabolic scaling
Ystep = Ystep .* 2.^(Xstep/2);

% the transform
LenaCL = SigmaTransform2D(  ...
    Lena , ...                  % the image
    win  , ...                  % the window
    Xstep(:), Ystep(:) , ...    % the channels
    sigmaX ,  sigmaY , ...      % the diffeomorphism (x,y) |-> (sigmaX(x,y),sigmaY(x,y))
    WX,  WY ,...                % the Fourier-domain
    actX , actY ...             % the translational action (if non-abelian)
);

% Save Reconstruction
recLenaCL = LenaCL.reconstruct('resid');

% Plot reconstruction
LenaCL.emergingRec( 0.1 , 1 );
pause;

% plot full Reconstruction
imagesc( real(recLenaCL )), axis ij; grid off, axis off;
title('Reconstruction, using dual frame and residuum');
pause;

%% 2D-SHEAR, Lena 128x128, rectangular window
sigmaX  = @(x,y) log2(abs(x) + eps);
sigmaY  = @(x,y) y./x;

actX    = @(x,y,xp,yp)                x - xp;
actY    = @(x,y,xp,yp) 2.^(-xp/2) .*( y - yp );

% load image of lena
[Lena,~]    = lena();

% define time- and frequency planes
t = linspace( 0, 1 , 128+1 ); t=t(1:end-1);
[X,Y] = meshgrid(t);
[WX,WY] = space2freq2D( X , Y );
WX = fftshift(WX); WY = fftshift(WY);

% define rectangular window  (of "width" 2 x 1 )
win = @(X,Y)   ( Y <= 1/2) .* ( Y > -1/2) ...
            .* ( X >   0 ) .* ( X <=  2 );

% make channels
 xstep = (0: 2 : 4);
 Xstep = []; Ystep = [];
 for k = 1 : length(xstep),
     ystep = -2*2.^(xstep(k)/2) : 1 : 2*2.^(xstep(k)/2);
     Xstep = [ Xstep ; repmat( xstep(k) , [ length( ystep ) , 1 ] ) ];
     Ystep = [ Ystep ; ystep(:) ];
 end;
% adjust step-width to parabolic scaling
Ystep = Ystep .* 2.^(Xstep/2);

% the transform
Lena2DShear = SigmaTransform2D(  ...
    Lena , ...                  % the image
    win  , ...                  % the window
    Xstep(:), Ystep(:) , ...    % the channels
    sigmaX ,  sigmaY , ...      % the diffeomorphism (x,y) |-> (sigmaX(x,y),sigmaY(x,y))
    WX,  WY ,...                % the Fourier-domain
    actX , actY ...             % the translational action (if non-abelian)
);

% Save Reconstruction
recLena2DShear = Lena2DShear.reconstruct('resid');

% plot reconstruction
Lena2DShear.emergingRec( .1 , 1 );
pause;

% plot full Reconstruction
imagesc( real(recLena2DShear )), axis ij; grid off, axis off;
title('Reconstruction, using dual frame and residuum');
pause;
