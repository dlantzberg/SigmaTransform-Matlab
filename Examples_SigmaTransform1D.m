%   Some usage examples for SigmaTransform1D()
%
%   AUTHOR: D Lantzberg, 2017 - 2018

% is this octave or matlab?
if(exist('OCTAVE_VERSION', 'builtin'))
    disp('Octave detected, trying to load signal package ..'); 
    disp('.. if it is not installed, you may do so with: ');
    disp('    pkg install signal -forge');
    pkg load signal;
end;

close all;

%% 1D STFT; Bat-Signal, Warped Gaussian window

%%% MAKE DATA %%%

% define diffeomorphism
sigma   = @(x) x;
action  = @(x,xp) x - xp;

% load bat-signal
[bat_signal,Fs] = bat();

% 400 linear spaced sampling points in warped domain
% steps = linspace( sigma(min_freq) , sigma(max_freq) , num_steps );
num_steps = 400;
domwidth = Fs;

% width of the window
win_width = 8; % in sampling steps

% make the window
warpedG   = @(x)  exp(  -pi * ( x/win_width * (num_steps/domwidth)  ).^2 );

% the transform
batSTFT = SigmaTransform1D( ...
    bat_signal, ... % the signal
    warpedG,    ... % the window - alternatively: the win_width for warped Gaussian as a number
    num_steps,  ... % the number of channels
    sigma,      ... % the diffeomorphism
    Fs,         ... % the sampling frequency
    action      ... % the "group action" (optional)
);

% save reconstruction
recbatSTFT = batSTFT.reconstruct('dual');

%%% PLOT DATA %%%

figure(1),shg;
% plot coeffs
subplot(121);
batSTFT.plotFrameogram('Frameogram for \sigma(x) = x (Spectrogram)');
pause(.1);

% plot Windows
subplot(122);
batSTFT.plotWindows('Windows in Fourier domain and diffeomorphism (dashed) for \sigma(x) = x');
pause(.1);pause();pause(.1);

figure(2),shg; pause(.1);
% plot emerging reconstruction
batSTFT.emergingRec( .001 , 2 );
pause;

% plot perfect reconstruction
tt = 0 : 1/Fs : (400-1)/Fs;
plot(tt,bat_signal,'b-',tt,real(recbatSTFT),'r--','LineWidth',1);
title('Reconstruction, using dual-frame in blue; original in dashed-red');
xlabel('time t \rightarrow' ); ylabel('f' );
axis tight; plotaxis; grid on;
pause;

%% 1D Wavelet; Bat-Signal, Warped Gaussian window

%%% MAKE DATA %%%

% define diffeomorphism, its inverse and jacobian
sigma   = @(x) log2(abs(x)+eps);
action  = @(x,xp) x - xp;

% load bat-signal
[bat_signal,Fs] = bat();

% 400 linear spaced sampling points in warped domain
num_steps = 400;
chan = linspace( sigma(Fs*0.005) , sigma(Fs/2*1.1) , num_steps );
domwidth = chan(end) - chan(1);

% width of the window
win_width = 8; % in sampling steps

% make the window
warpedG   = @(x)  exp(  -pi * ( x/win_width * (num_steps/domwidth)  ).^2 );

% the transform
batWT = SigmaTransform1D(   ...
    hilbert(bat_signal),    ... % the signal
    warpedG,                ... % the window
    chan,                   ... % the channels
    sigma,                  ... % the diffeomorphism
    Fs,                     ... % the sampling frequency
    action                  ... % the "group action" (optional)
);


%%% Alternative %%%

% define some function handle for a "named" transform like the "WaveletTransform" ...
WaveletTransform = @( signal , channels , Fs ) ...
    SigmaTransform1D( signal , warpedG , channels, @(x) log2(abs(x)+eps) , Fs );

% ... and apply it
batWT = WaveletTransform(   ...
    hilbert(bat_signal),    ... % the signal
    chan,                   ... % the channels
    Fs                      ... % the sampling frequency
);

%%%%%%%%%%%%%%%%%%%


% save reconstruction
recbatWT = batWT.reconstruct( 'resid' );

%%% PLOT DATA %%%

figure(1),shg;
% plot coeffs
subplot(121);batWT.plotFrameogram('Frameogram for \sigma(x) = log_2|x| (Scalogram)');

% plot Windows
subplot(122);
batWT.plotWindows('Windows in Fourier domain and diffeomorphism (dashed) for \sigma(x) = log_2|x|');
pause;

% plot emerging reconstruction
batWT.emergingRec( .001 , 2 );
pause;

% plot perfect reconstruction
tt = 0 : 1/Fs : (400-1)/Fs;
plot(tt,bat_signal,'b-',tt,real(recbatWT),'r--','LineWidth',1);
title('Reconstruction, using dual-frame and residuum in blue; original in dashed-red');
xlabel('time t \rightarrow' ); ylabel('f' );
axis tight; plotaxis; grid on;
pause;

%% 1D ConstantQ-Transform; Bat-Signal, Warped Gaussian window

%%% MAKE DATA %%%

f_0 = 256; % CenterFrequency in Hz
B   = 24;   % Bins per octave

sigma  = @(x) log2( abs(x)/f_0 + eps) * B;
action  = @(x,xp) x - xp;

% load bat-signal
[bat_signal,Fs] = bat();

% 400 linear spaced sampling points in warped domain
num_steps = 400;
chan = linspace( sigma(Fs*0.01) , sigma(Fs/2*1.1) , num_steps );

% the transform
batCQ = SigmaTransform1D( ...
    hilbert(bat_signal),    ... % the signal
    8,                      ... % the window
    chan,                   ... % the channels
    sigma,                  ... % the diffeomorphism
    Fs,                     ... % the sampling frequency
    action                  ... % the "group action" (optional)
);

% save reconstruction
recbatCQ = batCQ.reconstruct( 'resid' );

%%% PLOT DATA %%%

figure(1),shg;
% plot coeffs
subplot(121);
batCQ.plotFrameogram('Frameogram for \sigma_{CQ} (CQ-o-gram?)');

% plot Windows
subplot(122);
batCQ.plotWindows('Windows in Fourier domain and diffeomorphism (dashed) for \sigma_{CQ}');
pause;

% plot emerging reconstruction
batCQ.emergingRec( .001 , 2 );
pause;

% plot perfect reconstruction
tt = 0 : 1/Fs : (400-1)/Fs;
plot(tt,bat_signal,'b-',tt,real(recbatCQ),'r--','LineWidth',1);
title('Reconstruction, using dual-frame and residuum in blue; original in dashed-red');
xlabel('time t \rightarrow' ); ylabel('f' );
axis tight; plotaxis; grid on;
pause;

%% 1D ERBLet-Transform; Bat-Signal, Warped Gaussian window

%%% MAKE DATA %%%

a = 9.265;
b = 1;
c = 228.8445;

sigma   = @(x) a * log(b + c*abs(x));
action  = @(x,xp) x - xp;

% load bat-signal
[bat_signal,Fs] = bat();

% 400 linear spaced sampling points in warped domain
num_steps = 400;
chan = linspace( sigma(Fs*0.01) , sigma(Fs/2*1.1) , num_steps );
domwidth = chan(end) - chan(1);

% width of the window
win_width = 12; % in sampling steps

% make the window
warpedG   = @(x)  exp(  -pi * ( x/win_width * (num_steps/domwidth)  ).^2 );

% the transform
batERB = SigmaTransform1D( ...
    hilbert(bat_signal),    ... % the signal
    warpedG,                ... % the window
    chan,                   ... % the channels
    sigma,                  ... % the diffeomorphism
    Fs,                     ... % the sampling frequency
    action                  ... % the "group action" (optional)
);

% save reconstruction
recbatERB = batERB.reconstruct('resid');

%%% PLOT DATA %%%

figure(1),shg;
% plot coeffs
subplot(121);
batERB.plotFrameogram('Frameogram for \sigma_{ERB} (ERB-o-gram?)');

% plot Windows
subplot(122);
batERB.plotWindows('Windows in Fourier domain and diffeomorphism (dashed) for \sigma_{ERB}');
pause;

% plot emerging reconstruction
batERB.emergingRec( .001 , 2 );
pause;

% plot perfect reconstruction
tt = 0 : 1/Fs : (400-1)/Fs;
plot(tt,bat_signal,'b-',tt,real(recbatERB),'r--','LineWidth',1);
title('Reconstruction, using dual-frame and residuum in blue; original in dashed-red');
xlabel('time t \rightarrow' ); ylabel('f' );
axis tight; plotaxis; grid on;
pause;

%% 1D sigma_1(x) ~ 2*x + sin(x); Bat-Signal, Warped Gaussian window

%%% MAKE DATA %%%

% load bat-signal
[bat_signal,Fs] = bat();

% the transform
batMisc1 = SigmaTransform1D( ...
    hilbert(bat_signal),                    ... % the signal
    16,                                     ... % the window_width
    400,                                    ... % the number of channels
    @(x) 2*x + (Fs/26)*sin(2*pi*x/Fs*8),    ... % the diffeomorphism
    Fs                                      ... % the sampling frequency
);

% save reconstruction
recbatMisc1 = batMisc1.reconstruct( 'resid' );

%%% PLOT DATA %%%

figure(1),shg;
% plot coeffs
subplot(121);
batMisc1.plotFrameogram('Frameogram for \sigma_!');

% plot Windows
subplot(122);
batMisc1.plotWindows('Windows in Fourier domain and diffeomorphism (dashed) for \sigma_1');
pause;

% plot emerging reconstruction
batMisc1.emergingRec( .01 , 2 );
pause;

% plot perfect reconstruction
tt = 0 : 1/Fs : (400-1)/Fs;
plot(tt,bat_signal,'b-',tt,real(recbatMisc1),'r--','LineWidth',1);
title('Reconstruction, using dual-frame and residuum in blue; original in dashed-red');
xlabel('time t \rightarrow' ); ylabel('f' );
axis tight; plotaxis; grid on;
pause;

%% 1D sigma_2(x) ~ 1.25*x + cos(x); Bat-Signal, Warped Gaussian window

%%% MAKE DATA %%%

% load bat-signal
[bat_signal,Fs] = bat();

% define some transform, with fixed window(of width 16) and diffeomorphism, as a function handle...
SomeTransform = @( sig , numchan , Fs ) ...
    SigmaTransform1D( sig , 16 , numchan , @(x) 1.25*x +(Fs/20.5)*cos(2*pi*x/Fs*4) , Fs );

% Apply...
batMisc2 = SomeTransform(   ... % ...some transform
    hilbert(bat_signal),    ... % the signal
    400,                    ... % the number of channels
    Fs                      ... % the sampling frequency
);

% save reconstruction
recbatMisc2 = batMisc2.reconstruct('dual');

%%% PLOT DATA %%%

figure(1),shg;
% plot coeffs
subplot(121);
batMisc2.plotFrameogram('Frameogram for \sigma_"');

% plot Windows
subplot(122);
batMisc2.plotWindows('Windows in Fourier domain and diffeomorphism (dashed) for \sigma_2');
pause;

% plot emerging reconstruction
batMisc2.emergingRec( .01 , 2 );
pause;

% plot reconstruction
tt = 0 : 1/Fs : (400-1)/Fs;
plot(tt,bat_signal,'b-',tt,real(recbatMisc2),'r--','LineWidth',1);
title('Reconstruction, using dual-frame in blue; original in dashed-red');
xlabel('time t \rightarrow' ); ylabel('f' );
axis tight; plotaxis; grid on;
pause;
