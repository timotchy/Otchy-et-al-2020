function [SAPFeats, m_spec_deriv]  = aSAP_generateASAPFeatures(audio, fs, param)

if(~exist('param'))
    param = Parameters;
end

SAPFeats.lengthAudio = length(audio);
SAPFeats.audioSampRate = fs;

%compute the features...
chunkSize = 882000;
if(length(audio) <= chunkSize)
    [m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude ,gravity_center, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight ]= deriv(audio,fs);
else
    warning('aSAP_generateASAPFeatures: Not computing features because audio is longer than 882000 samples.');
    [m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude ,gravity_center, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight ] = deal([]);
    %for(nChunk = 1:ceil(length(audio)/chunkSize))
    %    startNdx = (nChunk-1) * chunkSize + 1;
    %    endNdx = min((nChunk-1) * chuckSize + chunkSize, length(audio))
    %    [t_spec_deriv , t_AM, t_FM ,t_Entropy , t_amplitude ,t_gravity_center, t_PitchGoodness , t_Pitch , t_Pitch_chose , t_Pitch_weight ]= deriv(audio(startNdx:endNdx),fs);
    %    append t values to m...
    %end
end
%Put the features into the struct
SAPFeats.param = param;
SAPFeats.m_AM = single(m_AM);
SAPFeats.m_FM = single(m_FM);
SAPFeats.m_Entropy = single(m_Entropy);
SAPFeats.m_amplitude = single(m_amplitude);
SAPFeats.gravity_center = single(gravity_center);
SAPFeats.m_PitchGoodness = single(m_PitchGoodness);
SAPFeats.m_Pitch = single(m_Pitch);
SAPFeats.Pitch_chose = single(Pitch_chose);
SAPFeats.Pitch_weight = single(Pitch_weight); 
SAPFeats.time = 0; %To be set elsewhere