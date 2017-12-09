function [filtered] = dataFilter(sig)
%filter all data for a patient

% Butterworth Highpass filter designed using FDESIGN.HIGHPASS.

% All frequency values are in Hz.
Fs = 0.38461538462;  % Sampling Frequency

Fstop = 0.005;   % Stopband Frequency
Fpass = 0.006;   % Passband Frequency
Astop = 80;      % Stopband Attenuation (dB)
Apass = 1;       % Passband Ripple (dB)
match = 'both';  % Band to match exactly

% Construct an FDESIGN object and call its ELLIP method.
h  = fdesign.highpass(Fstop, Fpass, Astop, Apass, Fs);
Hd = design(h, 'ellip', 'MatchExactly', match);

% filtered = zeros(size(sig,1),size(sig,1),1);
for i=1:1:size(sig,1)
    tac = sig(i).tac;
    filtered(i) = filter(Hd,tac);
end
end
