% Choose which figure to use based on the type of computer being used
if strcmp(computer, "PCWIN64") % machine is a PC
    TweetVision_PC
else % computer is a mac (or Linux machine; don't have special settings for that)
    TweetVision
end