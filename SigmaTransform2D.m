function [ out , PSIs ] = SigmaTransform2D( f , psi , Xsteps , Ysteps , sigmaX , sigmaY , wFsx , wFsy, actX, actY , detinvsigma )
%SigmaTransform2D 2D Continuous SigmaTransform, 
%   USAGE: [ out , PSIs ] = SigmaTransform2D( f , psi , Xsteps , Ysteps , sigmaX , sigmaY , wFsx , wFsy, actX, actY , detinvsigma )
%   INPUT:
%       f           : sig
%       psi         : wavelet; same size as sig; centerfreq & time at first index
%       X/Ystep     : 1-D vectors of steps in warped domain (or numsteps)
%       sigmaX/Y    : diffeomorphism as 2Dfunctions
%       wFs_x/y     : Samplingfrequency, if known
%       actX/Y      : optional, non-abelian shifting in warped domain
%       detinvsigma : determinant of jacobian of inverse of 2Ddiffeomorphism
%   OUTPUT:
%       out         : struct, containing W_psi f
%       [PSIs]      : [optional mat of fourier transforms of "scaled wavelets" ]
%
%   AUTHOR: D Lantzberg, Okt. 2016

    % config
    lensteps = length( Xsteps );
    
    if( lensteps ~= length(Ysteps) ),
        error('Stepvectors have different length');
    end;
    
    [height,width] = size( f );
    
    % Fourier domain provided or shall it be created?
    if(~exist('wFsx','var') || ~exist('wFsy','var') )
        wFsx = width;
        wFsy = height;
    end;
    if( length(wFsx) == 1 || length(wFsy) ==1 )
        warning('creating axis..');
        x = FourierAxis( wFsx , width );
        y = FourierAxis( wFsy , height );
        [WX,WY] = meshgrid( x , y );
    else
        WX = wFsx;
        WY = wFsy;
    end;

    % enough memory?
    if( log2(lensteps*height*width*8*2) >= 30 )
        warning('More than 1 GiByte RAM will be used. Sure? CTRL+C to cancel!');
        pause;
    end;
    if( log2(lensteps*height*width*8*2) >= 31 )
        warning('More than 2 GiByte RAM will be used. Sure?? CTRL+C to cancel!');
        pause;
    end;
    if( log2(lensteps*height*width*8*2) >= 32 )
        warning('More than 4 GiByte RAM will be used. Sure??? CTRL+C to cancel!');
        pause;
    end;
    if( log2(lensteps*height*width*8*2) >= 33 )
        warning('More than 8 GiByte RAM will be used. Sure???? CTRL+C to cancel!');
        pause;
    end;
    if( log2(lensteps*height*width*8*2) >= 34 )
        warning('More than 16 GiByte RAM would be used. Canceling.');
        return;
    end;
    % shall we scale the windows?
    if~exist('detinvsigma','var')
        scale = 0;
        detinvsigma = @(x,y) 1;
    else 
        scale = 1;
    end;
    
    % if no action given -> use abelian structure
    if~exist('actX','var')
        actX = @(x,y,xp,yp) x - xp;
    end;
    if~exist('actY','var')
        actY = @(x,y,xp,yp) y - yp;
    end;
        
    % get spectrum of signal
    f( isnan(f) )   = 0;
    F               = fft2( f );
        
    % diffeomorph Fourierdomain
    sigmaWX     = sigmaX(WX,WY); sigmaWX = repmat(sigmaWX,[1,1,lensteps]);
    sigmaWY     = sigmaY(WX,WY); sigmaWY = repmat(sigmaWY,[1,1,lensteps]);
    
    % make steps
    Xsteps      = reshape(Xsteps,1,1,lensteps);
    Ysteps      = reshape(Ysteps,1,1,lensteps);
    sigmaXSteps = repmat( Xsteps , [height , width , 1] );
    sigmaYSteps = repmat( Ysteps , [height , width , 1] );

    % evaluate window
    PSIs = psi( actX( sigmaWX, sigmaWY, sigmaXSteps , sigmaYSteps ) , ...
                actY( sigmaWX, sigmaWY, sigmaXSteps , sigmaYSteps ) );
    
    % handle NaNs
    PSIs(isnan(PSIs)) = 0;
    
    % Very Slow
    if( scale == 1 ),
        deter           = detinvsigma(2-Xsteps , 2-Ysteps )./detinvsigma(2,2);
        PSIs            = bsxfun( @times , deter.^.5 , PSIs );
        Mask            = sum( bsxfun( @times , deter.^-1 , abs(PSIs).^2 ) , 3 );
    else
        deter           = ones(size(Xsteps));
        Mask            = sum( abs(PSIs).^2 , 3 );
    end;
    
    % transform
    coeff = ifft2( bsxfun( @times, F , conj(PSIs) ) );
    
    % save residuum, if steps do not make up a frame
    resid = ifft2( F .* ( 1 - (Mask>eps) ) );
    
    % save and return stuff
    out = struct( ...
        'coeff'         , coeff  , ...      % transform coefficients
        'PSIs'          , PSIs , ...        % spectra of the windows
        'psi_hat'       , psi,  ...         % function handle of window in warped Fourier domain
        'omega_x'       , WX , ...          % X-components of the Fourier domain mesh-grid
        'omega_y'       , WY , ...          % Y-components of the Fourier domain mesh-grid
        'steps_x'       , Xsteps , ...      % used X-steps in warped domain
        'steps_y'       , Ysteps , ...      % used Y-steps in warped domain
        'sigma_x'       , sigmaX, ...       % function handle of first component of the spectral diffeomorphism
        'sigma_y'       , sigmaY, ...       % function handle of second component of the spectral diffeomorphism
        'warpFactors'   , deter, ...        % weighting factors
        'FourierMask'   , Mask , ...        % sum of the squared spectra of the windows
        'residuum'      , resid,  ...       % residuum (part of the signal which is lost due to missing windows)
        'reconstruct'   , @(modus) [], ...  % function handle for reconstruction
        'emergingRec'   , @(t,fig) t  ...   % function handle for plotting an emerging reconstruction
    );
    out.emergingRec = @(t,fig) plotReconstruction(out,t,fig);
    out.reconstruct = @(modus) invSigmaTransform2D(out,modus);
end

%% Plots the emerging Reconstruction
function t = plotReconstruction( c , t , fig )
    if~exist('t','var')
        % 50 ms by default
        t = 0.05;
    end;
    if~exist('fig','var')
        % default to figure 1
        fig = 1;
    end;
    % setup
    EmergingImage = 0;
    EmergingFourier = 0;
    figure(fig),colormap gray;
    emergImage    = subplot(121);
    title('Emerging Reconstruction','parent', emergImage );
    axis image; axis ij; axis off;
    set(gca,'NextPlot','replacechildren') ;
    emergFourier  = subplot(122);
    title('Emerging Fourier Plane','parent', emergFourier );
    axis image;
    xlabel('Wavenumber k_x \rightarrow','parent', emergFourier ); 
    ylabel('Wavenumber k_y \rightarrow','parent', emergFourier ); 
    set(gca,'NextPlot','replacechildren') ;
    shg;
    for k = 1 : 1 : length( c.steps_x(:) ),% : -1 : 1,
        % update images
        curr = c.PSIs(:,:,k);
        EmergingImage   = EmergingImage   + ifft2(fft2(c.coeff(:,:,k)).*curr);
        
        curr = fftshift(curr).^2;
        EmergingFourier = EmergingFourier + curr;
        
        % show updated images
        imagesc( abs(EmergingImage), 'parent', emergImage ), 
        imagesc( c.omega_x(1,:), c.omega_y(:,1), EmergingFourier + curr,'parent',emergFourier ), 
        
        % sleep for "t" seconds
        pause( t ); 
    end; 
end

%% inverse sigma transform
function [ recf ] = invSigmaTransform2D( C , dual )
%invSigmaTransform2D INVERSE Continuous 2D-sigmaTransform,  

    % reconstruct, using Frame
    if~exist('dual','var')
        disp('reconstructing using frame.');
        recf = ifft2( sum( bsxfun( @times, fft2( C.coeff ) .* C.PSIs, C.warpFactors.^-1 ) , 3 ) );
    else
        if( strcmp(dual,'dual') )
            disp('reconstructing using dualframe.');
            recf =  ifft2(sum( bsxfun( @times, fft2( C.coeff ) .* C.PSIs, C.warpFactors.^-1 ), 3 )...
                 .* StableInverse( ifftshift( C.FourierMask ) ) );
        % reconstruct, using dualFrame and "unresolvable" Residuum
        elseif( strcmp(dual,'resid') )
            disp('reconstructing using dualframe + residuum.');
            recf =  ifft2(sum( bsxfun( @times, fft2( C.coeff ) .* C.PSIs, C.warpFactors.^-1 ), 3 )...
                 .* StableInverse( ifftshift( C.FourierMask ) ) ) + C.residuum;
        else
            disp('reconstructing using frame.');
            recf = ifft2( sum( bsxfun( @times, fft2( C.coeff ) .* C.PSIs, C.warpFactors.^-1 ) , 3 ) );
        end 
    end
end

function [ f ] = StableInverse( f )
%StableInverse sets small values to 1 (for inversion)
    f( abs(f) < 100*eps ) = 1;
    f = 1./f;
end
