function output = DTAreader(FILENAME)
% Reads the .DTA files from the Gamry 600 potentiostat and output test
% information and data in a standard format for later processing.
% 
% Inputs:
% FILENAME - the absolute path for the file to be read
%
% Outputs:
% The output is a standardized structure:
%   filename -- name of the file read
%   testType -- identifier for the type of eChem test performed
%   date     -- date of the test
%	time     -- time of the test
%   notes    -- notes entered into the notes field
%   ocvcurve -- data from the pre-test
%   cvcurve  -- cyclic voltammetry data
%   eis      -- electrode impedance spectroscopy data
%
%Inspired by Yarden's github:
%https://github.com/gardner-lab/lab-operations/blob/master/utilities/Electrochemistry%20measurements/DTAread.m

%Set default input values
delim = '\t';

%Establish basic output structure
output = struct('filename', [], 'testType', [], 'date', [], 'time', [], 'notes', [], 'ocvcurve', [], 'cvcurve', [], 'eis', []);
    
%Parse out the filename
t = regexp(FILENAME, filesep, 'split');
output.filename = t{end};

%Open the file for binary read access
fid=fopen(FILENAME,'r');

%Read out and discard the first line
fgetl(fid);

%Retrieve the test type from second column
curLine=fgetl(fid);
t=regexp(curLine,delim,'split');
output.testType = t{2}; % Str: 'CV', 'EISOT', ?

%Read through all lines until the end of file flag is set
testNum = 1;
while ~feof(fid)
    %Read a single line
    curLine=fgetl(fid);
    
    %Parse out the tab-delimited column tokens
    t=regexp(curLine,delim,'split');
    
    %Process the line by switching between first-colum cases
    switch t{1}
        %Test date
        case 'DATE'
            output.date = t{3};
            
            %Test time of day
        case 'TIME'
            output.time = t{3};
            
            %Test notes
        case 'NOTES'
            if ~strcmp(t{4}, '&Notes...')
                output.notes = t{4};
            else
                output.notes = [];
            end
            
            %The marker for the preroutine data eChem runs
        case 'OCVCURVE'
            %Read the data block
            [data, fid] = readBlock(fid, delim);
            
            %Copy the data to output structure
            output.ocvcurve.time = data(:,2);
            output.ocvcurve.Vf = data(:,3);
            output.ocvcurve.Vm = data(:,4);
            
            %Marker for the start of a EIS datablock
        case 'ZCURVE'
            %Read the data block
            [data, fid] = readBlock(fid, delim);
            
            %Copy the data to output structure
            output.eis.time = data(:,2);
            output.eis.freq = data(:,3);
            output.eis.Zreal = data(:,4);
            output.eis.Zmod = data(:,7);
            output.eis.Zph = data(:,8);
            
            %The first column doesn't (by itself) contain unique identifying
            %data... check a little further
        otherwise
            %Slice off the first 5 characters and check the tag
            if numel(t{1})>=5
                r = t{1}(1:5);
            end
            
            if strcmp(r, 'CURVE')
                
                %Read the data block
                [data, fid] = readBlock(fid, delim);
                
                %Copy the data to output structure
                output.cvcurve(testNum).time = data(:,2);
                output.cvcurve(testNum).Vf = data(:,3);
                output.cvcurve(testNum).Im = data(:,4);
                
                %Update the testNum pointer
                testNum = testNum + 1;

            end
    end
end

%Close the file
fclose(fid);

function [data, fid] = readBlock(fid, delim)
%Read out the whole block of text

%Trash two lines
fgetl(fid); fgetl(fid);

%Flag for continue to read lines
bGo = true;
data = [];
while bGo && ~feof(fid)
    %Read in purported first dataline
    t=regexp(fgetl(fid),delim,'split');
    
    %If the first column is empty, we're in a data block
    if strcmp(t{1}, '')
        %Fix for avoiding conversion of an empty string to number
        %Tact the current line tokens onto the end of the datablock
        if isempty(str2num(t{end}))
            data=[data; cellfun(@str2num,t(2:end-1))];
        else
            data=[data; cellfun(@str2num,t(2:end))];
        end
    else
        %Exit the while loop
        bGo = false;
    end
end
