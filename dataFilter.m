function [filtered] = dataFilter(sig,varargin)
%filter all data for a patient
TR=2.6; %sec
tipo='';
filtered = zeros(size(sig,2),size(sig(1).tac,1));
if(nargin == 1)
   tipo = 'ellip'; 
else
    if(strcmp(varargin{1},'ellip') || strcmp(varargin{1},'butter'))
        tipo = varargin{1};
    else 
        err('Il tipo di filtro deve essere o "ellip" o "butter" ');
    end
end
% Butterworth Highpass filter designed using FDESIGN.HIGHPASS.
% All frequency values are in Hz.
if( strcmp(tipo,'ellip'))
    

    Fs = 1/TR;  % Sampling Frequency

    Fstop = 0.005;   % Stopband Frequency
    Fpass = 0.006;   % Passband Frequency
    Astop = 80;      % Stopband Attenuation (dB)
    Apass = 1;       % Passband Ripple (dB)
    match = 'both';  % Band to match exactly

    % Construct an FDESIGN object and call its ELLIP method.
    h  = fdesign.highpass(Fstop, Fpass, Astop, Apass, Fs);
    Hd = design(h, 'ellip', 'MatchExactly', match);
    parfor i=1:1:size(sig,1)
        filtered(i,:) = filter(Hd,sig(i).tac)';
    end
    
elseif  strcmp(tipo,'butter') %%use butterworth
    Fs=1/TR;
    Wp=0.007/(Fs/2);
    Ws=0.004/(Fs/2);
    Rp=0.5; %dB                   
    Rs=25; %dB
    [ord, Wn]=buttord(Wp,Ws,Rp,Rs);
    [B,A]=butter(ord,Wn,'high');
    parfor i=1:1:size(sig,2)
        filtered(i,:)=filtfilt(B,A,sig(i).tac)';
    end
else
    err('Il tipo di filtro deve essere o "ellip" o "butter" ');
end


end
