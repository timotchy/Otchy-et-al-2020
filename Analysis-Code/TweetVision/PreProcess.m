function [data] = PreProcess(data, PP)
% Takes in single channel dataset, recorded at the given sampling rate, and
% uses the parameters specified in the PreProcess structure to run
% sequential processing steps and return the output.

% Sequentially do each of the stages
for i=1:length(PP.Steps)
    if PP.Steps(i) == 1; % None
        
    elseif PP.Steps(i) == 2; % Local Detrend
        data = locdetrend(data,PP.params.Fs,[PP.Detrend.window PP.Detrend.winstep]);
        
    elseif PP.Steps(i) == 3; % Denoising
        params = PP.params;
        params.fpass = [PP.Denoise.LP PP.Denoise.HP];
        data = rmlinesc(data,params,PP.Denoise.p,PP.Denoise.plt);
    
    elseif PP.Steps(i) == 4; % High Pass
        % Prep the coeffients for the filter
        fNorm = PP.HP.cutoff/(PP.params.Fs/2);
        if strcmp('butter',PP.HP.type)
            [b,a] = butter(PP.HP.order,fNorm,'high');
        elseif strcmp('cheby1',PP.HP.type)
            [b,a] = cheby1(PP.HP.order,PP.HP.rpass_ripple,fNorm,'high');
        elseif strcmp('cheby2',PP.HP.type)
            [b,a] = cheby2(PP.HP.order,PP.HP.rstop_ripple,fNorm,'high');
        elseif strcmp('ellip',PP.HP.type)
            [b,a] = ellip(PP.HP.order,PP.HP.rpass_ripple,PP.HP.rstop_ripple,fNorm,'high');
        end
        
        % Run the specified zero-phase filter
        data = filtfilt(b,a,data);

    elseif PP.Steps(i) == 5; % Low Pass
         % Prep the coeffients for the filter
        fNorm = PP.LP.cutoff/(PP.params.Fs/2);
        if strcmp('butter',PP.LP.type)
            [b,a] = butter(PP.LP.order,fNorm,'low');
        elseif strcmp('cheby1',PP.LP.type)
            [b,a] = cheby1(PP.LP.order,PP.LP.rpass_ripple,fNorm,'low');
        elseif strcmp('cheby2',PP.LP.type)
            [b,a] = cheby2(PP.LP.order,PP.LP.rstop_ripple,fNorm,'low');
        elseif strcmp('ellip',PP.LP.type)
            [b,a] = ellip(PP.LP.order,PP.LP.rpass_ripple,PP.LP.rstop_ripple,fNorm,'low');
        end
        
        % Run the specified zero-phase filter
        data = filtfilt(b,a,data);               
                
    elseif PP.Steps(i) == 6; % Amplification
        data = data*PP.Amp.gain;
        
    elseif PP.Steps(i) == 7; % Smoothing
        span = PP.Smooth.window*PP.params.Fs;
        
        if strcmp('sgolay',PP.Smooth.type)
            data = smooth(data,span,PP.Smooth.type,PP.Smooth.degree);
        else
            data = smooth(data,span,PP.Smooth.type);
        end
        
    elseif PP.Steps(i) == 8; %Down Sampling
        factor = floor(PP.params.Fs/PP.DSample.target);
        data = downsample(data,factor);
        
    end

end