function patterns = PatternFile_test
%This script defines the set of stimulation patterns used in a Nanoclip
%experiment (2018). This protoype file (mainly for testing) is based on the
%LR82 experiments (4/7/2018).
%
% Nanoclip Pad Layout:
%     _______   _______   _______ 
%    |   1   | |   2   | |   6   |
%    |_______| |_______| |_______|
%     _______   _______   _______ 
%    |   3   | |   5   | |   4   |
%    |_______| |_______| |_______|
% 
% Modified by TMO 04/18/2018

%Define vars
patterns = [];
patterns.I1 = [];
patterns.I2 = [];
patterns.T1 = [];
patterns.T2 = [];
patterns.TI = [];
patterns.Reps = [];
patterns.Freq = [];

default = [];
default.I1 = [];
default.I2 = [];
default.T1 = [200, 200, 200, 200, 200, 200];
default.T2 = [200, 200, 200, 200, 200, 200];
default.TI = [5, 5, 5, 5, 5, 5];
default.Reps = [100, 100, 100, 100, 100, 100];
default.Freq = [1000, 1000, 1000, 1000, 1000, 1000];

%Pattern 1
i = 1;
patterns(i).I1 = [50, 50, 50, 50, 50, 50];
patterns(i).I2 = [-50, -50, -50, -50, -50, -50];
patterns(i).T1 = default.T1;
patterns(i).T2 = default.T2;
patterns(i).TI = default.TI;
patterns(i).Reps = default.Reps;
patterns(i).Freq = default.Freq;

%Pattern 2
i = 2;
patterns(i).I1 = [-50, -50, -50, -50, -50, -50];
patterns(i).I2 = [50, 50, 50, 50, 50, 50];
patterns(i).T1 = default.T1;
patterns(i).T2 = default.T2;
patterns(i).TI = default.TI;
patterns(i).Reps = default.Reps;
patterns(i).Freq = default.Freq;

%Pattern 3
i = 3;
patterns(i).I1 = [50, -100, 50, 50, -100, 50];
patterns(i).I2 = [-50, 100, -50, -50, 100, -50];
patterns(i).T1 = default.T1;
patterns(i).T2 = default.T2;
patterns(i).TI = default.TI;
patterns(i).Reps = default.Reps;
patterns(i).Freq = default.Freq;

%Pattern 4
i = 4;
patterns(i).I1 = [-50, 100, -50, -50, 100, -50];
patterns(i).I2 = [50, -100, 50, 50, -100, 50];
patterns(i).T1 = default.T1;
patterns(i).T2 = default.T2;
patterns(i).TI = default.TI;
patterns(i).Reps = default.Reps;
patterns(i).Freq = default.Freq;

%Pattern 5
i = 5;
patterns(i).I1 = [100, 0, 100, -100, 0, -100];
patterns(i).I2 = [-100, 0, -100, 100, 0, 100];
patterns(i).T1 = default.T1;
patterns(i).T2 = default.T2;
patterns(i).TI = default.TI;
patterns(i).Reps = default.Reps;
patterns(i).Freq = default.Freq;

%Pattern 6
i = 6;
patterns(i).I1 = [-100, 0, -100, 100, 0, 100];
patterns(i).I2 = [100, 0, 100, -100, 0, -100];
patterns(i).T1 = default.T1;
patterns(i).T2 = default.T2;
patterns(i).TI = default.TI;
patterns(i).Reps = default.Reps;
patterns(i).Freq = default.Freq;

%Pattern 7
i = 7;
patterns(i).I1 = [0, 100, 0, -100, 100, -100];
patterns(i).I2 = [0, -100, 0, 100, -100, 100];
patterns(i).T1 = default.T1;
patterns(i).T2 = default.T2;
patterns(i).TI = default.TI;
patterns(i).Reps = default.Reps;
patterns(i).Freq = default.Freq;

%Pattern 8
i = 8;
patterns(i).I1 = [-100, 100, -100, 0, 100, 0];
patterns(i).I2 = [100, -100, 100, 0, -100, 0];
patterns(i).T1 = default.T1;
patterns(i).T2 = default.T2;
patterns(i).TI = default.TI;
patterns(i).Reps = default.Reps;
patterns(i).Freq = default.Freq;

%Pattern 9
i = 9;
patterns(i).I1 = [0, 0, 150, 0, 0, -150];
patterns(i).I2 = [0, 0, -150, 0, 0, 150];
patterns(i).T1 = default.T1;
patterns(i).T2 = default.T2;
patterns(i).TI = default.TI;
patterns(i).Reps = default.Reps;
patterns(i).Freq = default.Freq;

%Pattern 10
i = 10;
patterns(i).I1 = [150, 0, 0, -150, 0, 0];
patterns(i).I2 = [-150, 0, 0, 150, 0, 0];
patterns(i).T1 = default.T1;
patterns(i).T2 = default.T2;
patterns(i).TI = default.TI;
patterns(i).Reps = default.Reps;
patterns(i).Freq = default.Freq;










