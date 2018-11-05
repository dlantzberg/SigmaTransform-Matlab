function [ out , PSIs ] = SigmaTransform1D( f , psi , steps , sigma , wFs , action , detinvsigma )
%SigmaTransform1D 1D Continuous SigmaTransform	
%   USAGE: [ out , PSIs ] = SigmaTransform1D( f , psi , steps , sigma , wFs , action, detinvsigma )
%	INPUT:
%		f       : sig - or "0", to indicate calculation of windows, only
%		psi     : waveletfunc, in the warped domain
%		steps   : vector of steps to analyze in warped domain, or numsteps
%		sigma   : diffeomorphism as a function
%		wFs     : Samplingfrequency, if known, or FourierAxis
%		action  : translation action (optional, uses abelian by default)
%		detinv. : determinant of jacobian of inverse of diffeomorphism
%	OUTPUT:
%		out     : struct, containing W_psi f
%		[PSIs]	: [optional mat of fourier transforms of "scaled wavelets" ]
%
%	AUTHOR:	Daniel Lantzberg, Okt. 2016

	% config
	lensteps	= length( steps );
	lenf        = length( f );
	
	% ERRORCHECKING
	if~exist('detinvsigma','var')
		scale = 0;
        detinvsigma = @(x) 1;
    else 
        scale = 1;
	end;
	
	if~exist('wFs' , 'var')
		%warning('no axis given - using axis -1:1');
        wFs = lenf/2;
	end;
	    
    if~exist('action','var')
        % if no action given -> use abelian structure
		%warning('no action given - using abelian action');
		action = @(x,xp) x - xp;
    end;

	% signal given?
    if( lenf == 1 ) % get PSIs only		    
        
        if( length(wFs) == 1 )
            error('without signal, axis is needed');
        end;

        w       = reshape( wFs    , 1  , [] );	
        sigmaw  = sigma(w);

        % number of sample-points
        if( lensteps == 1 ) 
            lensteps = steps;
            steps = linspace( min(sigmaw) , max(sigmaw) , lensteps );
        end;	
        steps   = reshape( steps , [] , 1  );	
        
        if(~isa(psi,'function_handle') )
            % if no function handle is given: use warped Gaussian of
            % "width psi"
            width = psi;
            psi = @(x)  exp(  -pi * ( x/width * (lensteps/(steps(end)-steps(1)))  ).^2 );
        end;
        
        vars    = bsxfun( action , sigmaw , steps );
        PSIs    = psi( vars );    
        PSIs(isnan(PSIs)) = 0;    

        %scaling? -> takes time
        if( scale == 1 )
            deter   = detinvsigma(w);
            %deter  = bsxfun(@rdivide,detinvsigma(w),detinvsigma(vars));
            PSIs    =      bsxfun( @times , deter.^.5 , PSIs );
            Mask    = sum( bsxfun( @times , deter.^-1 , abs(PSIs).^2 ) , 1);
        else
            deter   = ones(size(steps(:)));
            Mask    = sum( abs(PSIs).^2 , 1);
        end;
        
        % irrelevant
        FF = 0;
        resid = 0;
    else
        % make axis, if wFs is sampling frequency (in timedomain)
        if( length(wFs) == 1 )
            wFs    = FourierAxis( wFs , lenf );
        end;
        if(lenf ~= length(wFs) )
            error('signal length doesnt match axis length');
        end;

        % handle nans
        f(isnan(f)) = 0;
        
        % make axes
        w       = reshape( wFs    , 1  , [] );	
        sigmaw  = sigma(w);
        % make infinity huge, but finite
        sigmaw( isinf(sigmaw) ) = .01/eps;

        % number of sample-points
        if( lensteps == 1 ) 
            lensteps = steps;
            steps    = linspace( min(sigmaw) , max(sigmaw) , lensteps-4 );
            stepsize = steps(2)-steps(1);
            steps    = [ steps(1)-stepsize*2, steps(1)-stepsize , steps , ...
                         steps(end)+stepsize , steps(end)+stepsize*2 ];
        end;		
        steps   = reshape( steps , [] , 1  );	
        
        if(~isa(psi,'function_handle') )
            % if no function handle is given: use warped Gaussian of
            % "width psi"
            width = psi;
            psi = @(x)  exp(  -pi * ( x/width * (lensteps/(steps(end)-steps(1)))  ).^2 );
        end;
        
        % make windows
        vars    = bsxfun( action , sigmaw , steps );		
        PSIs    = psi( vars );    
        PSIs(isnan(PSIs)) = 0;  
        
        % get spectrum of signal
        F     = reshape( fft( f ) , 1 , [] );

        % Very slow
        if( scale == 1 )
            deter       = bsxfun(@rdivide,detinvsigma(w),detinvsigma(vars));
            PSIs		=      bsxfun( @times , deter.^.5 , PSIs );
            Mask        = sum( bsxfun( @times , deter.^-1 , abs(PSIs).^2 ) , 1);
        else
            deter		= eye(size(vars));
            Mask        = sum( abs(PSIs).^2 , 1);
        end;

        % transform
        FF = ifft( bsxfun( @times , F , conj(PSIs) ) , [] , 2 );
        
        % save residuum, if steps do not make up a frame
        resid = ifft( F .* ( 1 - (Mask>eps) ) );
    end;
	out  = struct( ...
        'coeff'			    , FF, ...               % transform coefficients
        'PSIs'			    , PSIs , ...            % spectra of the windows
        'psi_hat'		    , psi,	...             % function handle of window in warped Fourier domain
        'omega'			    , w , ...               % the Fourier domain vector
        'steps'			    , steps , ...           % used steps in warped domain
        'sigma'			    , sigma, ...            % function handle of the spectral diffeomorphism
        'warpFactors'   	, deter, ...            % weighting factors
        'FourierMask'   	, Mask , ...            % sum of the squared spectra of the windows
        'residuum'      	, resid,  ...           % residuum (part of the signal which is lost due to missing windows)
        'reconstruct'   	, @( modus ) [], ...    % function handle for reconstruction
        'emergingRec'   	, @(t,fig) [], ...      % function handle for plotting an emerging reconstruction
        'plotFrameogram'	, @( titlestr ) [], ... % function handle for plotting frame-o-gram
        'plotWindows'   	, @( titlestr ) [] ...  % function handle for plotting the windows
	);
    % use "this" structure
    out.emergingRec     = @(t,fig)      plotReconstruction(out,f,t,fig);
    out.plotFrameogram  = @(titlestr)   plotFrameogram(out,titlestr);
    out.plotWindows     = @(titlestr)   plotWindows(out,titlestr);
    out.reconstruct     = @(modus)      invSigmaTransform1D(out,modus);
end

%% Plots the Frame-o-gram
function [] = plotReconstruction( c , f , t , fig )
    if~exist('t','var')
        % 50 ms by default
        t = 0.05;
    end;
    if~exist('fig','var')
        % default to figure 1
        fig = 1;
    end;
    
    % setup
    x = freq2time( ifftshift(c.omega) );
    w = fftshift(c.omega);
    EmergingSignal = 0;
    EmergingFourier = 0;
    figure(fig);shg;
    pause(.1);
    
    emergSignal  = subplot(211); pause(.01);
    plot( x , norm1(real(EmergingSignal),inf) ); axis tight; plotaxis;
    axis( [ x(1), x(end) , -1 , 1 ] );
    title('Emerging Reconstruction, using frame in blue; original in dashed-red','parent', emergSignal );grid on;
    xlabel('time t \rightarrow','parent', emergSignal ); 
    ylabel('f','parent', emergSignal );
    plotaxis; set(gca,'NextPlot','replacechildren') ;
    pause(.1);
    
    emergFourier = subplot(212); pause(.01);
    plot( w , norm1(fftshift(EmergingFourier),inf), 'k' ); axis tight;  plotaxis;
    axis( [ w(1), w(end) , 0 , 1 ] );
    title('Emerging Fourier Axis','parent', emergFourier ); grid on;
    xlabel('Freq \omega \rightarrow','parent', emergFourier ); 
    ylabel('|SPEC(\psi)|^2','parent', emergFourier );
    set(gca,'NextPlot','replacechildren') ;
    pause(.1);
    

    for k = length( c.steps ) : -1 : 1,
        % update images
        curr = c.PSIs(k,:);
        EmergingSignal  = EmergingSignal  + ifft(fft(c.coeff(k,:)).*curr);
        EmergingFourier = EmergingFourier + curr.^2;
        
        % show updated plots
        plot( x , norm1(real(EmergingSignal),inf), ...
              x , norm1(real(f),inf), 'r--','parent', emergSignal, 'LineWidth', 1  );
        
        plot( w, norm1(fftshift(EmergingFourier),inf), 'k' , ...
              w, norm1(fftshift(curr),inf), 'r--' , 'parent', emergFourier, 'LineWidth',1 );
        
        % sleep for "t" seconds
        pause( t ); 
    end; 
end

%% Plots the Frameogram
function [] = plotFrameogram( c , titlestr )
    imagesc( freq2time( ifftshift(c.omega) ) , c.steps , abs( c.coeff ).^2 )
    %axis square, 
    axis xy,axis tight,plotaxis(0),colormap jet;
    xlabel('t \rightarrow');
    ylabel('\sigma(\omega) \rightarrow');
    title( titlestr ); grid on;
    shg;
end

%% Plots the windows
function [] = plotWindows( c , titlestr )
    imagesc( ifftshift(c.omega) , c.steps , fftshift(abs( c.PSIs ),2) )
    %axis square, 
    axis xy,axis tight,plotaxis(0),colormap jet;
    xlabel('\omega \rightarrow');
    ylabel('\sigma(\omega) \rightarrow');
    title(titlestr); grid on;
    % plot the diffeomorphism
    hold on,
    plot( ifftshift(c.omega) , c.sigma( ifftshift(c.omega) ) , 'w--' , 'LineWidth' , 1.1);
    hold off;
    shg;
end

%% reconstructs from coefficients
function [ recf ] = invSigmaTransform1D( C  , dual )
	% reconstruct, using Frame
	if(~exist('dual','var'))
        disp('reconstructing using frame.');
        recf = ifft( sum(C.warpFactors^-1 * ( fft(C.coeff,[],2) .* C.PSIs ) , 1 ) );		
    else	
        % reconstruct, using dualFrame
        if( strcmp(dual,'dual') )
            disp('reconstructing using dualframe.');
            recf =  ifft( sum(C.warpFactors^-1 * ( fft(C.coeff,[],2) .* C.PSIs ) , 1 ) ...
                 .* StableInverse( C.FourierMask ) );		 
        % reconstruct, using dualFrame and "unresolvable" Residuum
        elseif( strcmp(dual,'resid'))
            disp('reconstructing using dualframe and residuum.');
            recf =  ifft( sum(C.warpFactors^-1 * ( fft(C.coeff,[],2) .* C.PSIs ) , 1 ) ...
                 .* StableInverse( C.FourierMask ) ) + C.residuum;		
        else
            disp('reconstructing using frame.');
            recf = ifft( sum(C.warpFactors^-1 * ( fft(C.coeff,[],2) .* C.PSIs ) , 1 ) );	
        end;
    end;
end

function [ f ] = StableInverse( f )
%StableInverse sets small values to 1 (for inversion)
    f( abs(f) < 100*eps ) = 1;
    f = 1./f;
end
