function [uNoise, uSound, sdNoise, sdSound] = aSAP_estimateTwoMeans(audioLogPow)    

%Run EM algorithm on mixture of two gaussian model:
    
%set initial conditions
l = length(audioLogPow);
len = 1/l;
m = sort(audioLogPow);
uNoise = mean(m(1:round(end/2)));
uSound = mean(m(round(end/2):end));
sdNoise = .5;
sdSound = 2;

%compute estimated log likeilyhood given these initial conditions...
prob = zeros(2,l);
prob(1,:) = (exp(-(audioLogPow - uNoise).^2 / (2*sdNoise^2)))./sdNoise;
prob(2,:) = (exp(-(audioLogPow - uSound).^2 / (2*sdSound^2)))./sdSound;
[estProb, class] = max(prob);
logEstLike = sum(log(estProb)) * len;        
logOldEstLike = -Inf;
    
%maximize using Estimation Maximization
while(abs(logEstLike-logOldEstLike) > .005)
    logOldEstLike = logEstLike;

    %Which samples are noise and which are sound.
    nndx = find(class==1);
    sndx = find(class==2);

    %Maximize based on this classification.
    uNoise = mean(audioLogPow(nndx));
    sdNoise = std(audioLogPow(nndx));
    uSound = mean(audioLogPow(sndx));
    sdSound = std(audioLogPow(sndx));

    %Given new parameters, recompute log liklyhood.
    prob(1,:) = (exp(-(audioLogPow - uNoise).^2 / (2*sdNoise^2)))./sdNoise;
    prob(2,:) = (exp(-(audioLogPow - uSound).^2 / (2*sdSound^2)))./sdSound;
    [estProb, class] = max(prob);
    logEstLike = sum(log(estProb)) * len;       
end

