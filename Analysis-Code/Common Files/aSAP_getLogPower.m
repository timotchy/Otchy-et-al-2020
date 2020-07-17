function [sigLogPow, filteredAudio] = aSAP_getLogPower(sig, fs)

%historic code used variable name audio instead of sig.
audio = sig;

%demean the audio.
audio = audio - mean(audio);

%include only very relevant power:
%below 8000Hz
filt2.order = 50; %sufficient for 44100Hz of lower
filt2.win = hann(filt2.order+1);
filt2.cutoff = 6000; %Hz
filt2.fs = fs;
filt2.lpf = fir1(filt2.order, filt2.cutoff/(filt2.fs/2), 'low', filt2.win);
audio = filtfilt(filt2.lpf, 1, audio);

%above 860Hz
filt3.order = 50; %sufficient for 44100Hz of lower
filt3.win = hann(filt3.order+1);
filt3.cutoff = 860; %Hz
filt3.fs = fs;
filt3.hpf = fir1(filt3.order, filt3.cutoff/(filt3.fs/2), 'high', filt3.win);
audio = filtfilt(filt3.hpf, 1, audio);

%compute power
audioPow= audio.^2; 

%smooth the power, lpf:
filt.order = 100; 
filt.win = hann(filt.order+1);
filt.cutoff = 50; %Hz
filt.fs = fs;
filt.lpf = fir1(filt.order, filt.cutoff/(filt.fs/2), 'low', filt.win);
audioPow = filtfilt(filt.lpf, 1, audioPow);

%compute log Pow
audioLogPow = log(audioPow + eps);

filteredAudio = audio;
sigLogPow = audioLogPow;