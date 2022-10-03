function [y] = delayfm(t,x,delay)
% Delay a matrix of time series data x(Nt,Ne) by a vector of time delays (Ne)
% M Gray, updated 08 Sep 2022

%% Check data matrix size and shape
t = t(:);                   % Time vector (s)
fs = 1/mean(diff(t));       % Sample rate (Hz)
Nt = length(t);             % Length of time vector
Ne = length(delay);         % Length of delays/data channels

if length(x(:,1))~=Nt
    swp = 1;
    x = transpose(x);
else
    swp = 0;
end


%% Apply delays in frequency domain, then invert back to time. Positive values are delays, negative values are leads
fx = fft(x,[],1); 

if ~rem(Nt,2)  % even length series

    iwt = repmat(-2i*pi*(0:Nt/2)'*fs/Nt,[1 Ne]).*repmat(delay(:)',[(Nt/2+1) 1]);    % argument for positive frequencies
    tmp = fx(1:(Nt/2+1),:).*exp(iwt);                                               % delayed positive frequencies
    tmp((Nt/2+2):Nt,:) = conj(tmp(Nt/2:-1:2,:));                                    % delayed negative frequencies
                       
else           % odd length series
    
    lt2 = floor(Nt/2);
    iwt = repmat(-2i*pi*(0:lt2)'*fs/Nt,[1 Ne]).*repmat(delay(:)',[(lt2+1) 1]);
    
    tmp = fx(1:(lt2+1),:).*exp(iwt);
    tmp((lt2+2):Nt,:) = conj(tmp((lt2+1):-1:2,:));
                    
end

y = real(ifft(tmp,[],1));                                                           % return to time domain


%% Return data to original input orientation if needed
if swp
    y = transpose(y);
end