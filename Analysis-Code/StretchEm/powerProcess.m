% This function processes the input power spectrum signal, by 
% 1) band pass filtering
% 2) clipping the signal to the start and end



function [Pout Fout T] = powerProcess(y,Fs)

    [S F T P] = spectrogram(y,ones(441,1)',400,513, Fs);
    

    [dum sIndex] = min(abs(F-1000));
    [dum eIndex] = min(abs(F-10000));


    Pout = P(sIndex:eIndex,:);

    Fout = F(sIndex:eIndex);
    
    
    

end