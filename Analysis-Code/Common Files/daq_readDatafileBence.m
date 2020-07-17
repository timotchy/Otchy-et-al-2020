function [HWChannels, data, time, startSamp]  = daq_readDatafileBence(filename)
%Read a data file created by a trigger.  

%As of file format '-2', if two saving sessions are appended into the 
%same file, this function will return the first session.

%filename:  string specifing the filename of the data file.

%HWChannels:  the hardware channel numbers of the columns of the data
%matrix.

%data:  matrix in which each column contains data from the specified HW
%channel number in HWChannels.

%time:  In file format -1 and -2, column vector containing the time relative to daq_Start each row of
%                                 data was taken.
%       In file format -3: a vector containing
%       [acquisitionStartTime (6 element datevec), approxFileCreatedTime (6 element datevec), startSampTime (seconds), endSampTime (seconds), endSampNumber (integer)]

%%%%NOTE:  DO NOT CHANGE THE FILE FORMAT, 
%%%%       WITHOUT INCREMENTING THE TRIGGER FILE FORMAT ID NUMBER!!!! 
startSamp=[];
time=[];
if(~exist(filename,'file'))
    error('File does not exist');
end

%Open the file
fid = fopen(filename);
%Read first element in the file which is the file format
trigFileFormat = fread(fid, 1, 'float64');

if(trigFileFormat == -3)
    
    %Read absolute clock times in datevec format, 6 number [Y,M,D,H,M,S]. 
    acquisitionStartTime = fread(fid, 6, 'float64'); %Time the daq was started.
    approxFileCreatedTime = fread(fid, 6, 'float64'); %Approximate time this file was writtin.
    
    %Read the number of hardware channels
	numHWChannels = fread(fid, 1, 'float64');

    %Read the hardware indices of the recorded channels.
	HWChannels = fread(fid, numHWChannels, 'float64');
    
    %Read information for converting data from native format to doubles
    %double data= (native data)*(native scaling constant) + (native offset constant)
    for(nChan = 1:numHWChannels)
        nativeScale(nChan) =  fread(fid, 1, 'float64');
        nativeOffset(nChan) =  fread(fid, 1, 'float64');
    end
    
    %Read the name of the native data format.
    nLetter = 1;
    while(true)
        c = fread(fid, 1, 'float64');
        if(c == trigFileFormat)
            break;
        end
        nativeDataType(nLetter) = c;
        nLetter = nLetter + 1;
    end
    nativeDataType = char(nativeDataType);
    
    %Read internal daq timing for first sample in file.
    startSampNumber = fread(fid, 1, 'float64'); %The sample number.  First sample after daq was started is sample #1.
    startSampTime = fread(fid, 1, 'float64'); %The time the sample number in seconds from first sample. The function getdata return a time for each sample.

    %Read the remainder of the file in native data type
    data = fread(fid, inf, nativeDataType);
    
    %fileFormatID is used to mark the beginning and end of each recording
    %session.  Use this mark to locate the end of the data block.
    mark = find(data == -3);
    
	%Reshape the datafile so that each column represents a channel.
	data = reshape(data(1:mark(end)-1), numHWChannels, length(data(1:mark(end)-1))/(numHWChannels))';
	 
    %Convert the daq native format to a double
    for(nChan = 1:numHWChannels)
        data(:,nChan) = data(:,nChan)*nativeScale(nChan) + nativeOffset(nChan);
    end
    
    %back up from end of file by two double.
    fseek(fid, -16, 'eof');
    
    %Read internal daq timing for last sample in file.
    endSampNumber = fread(fid, 1, 'float64');
    endSampTime = fread(fid, 1, 'float64');

    %close the file.
    fclose(fid);    
    
    %Create time array:
    %time = [acquisitionStartTime', approxFileCreatedTime', startSampTime, endSampTime, endSampNumber];
elseif(trigFileFormat == -2)
    %read remainder of file.
    data = fread(fid, inf, 'float64');
    fclose(fid);
        
    %fileFormatID is used to mark the beginning and end of each recording
    %session.  Use this mark to locate the end of the data block.
    mark = find(data == -2);
    
 	%Parse the first float64 as the number of channels
	startSampNumber = data(1);   
    
	%Parse the second float64 as the number of channels
	numHWChannels = data(2);
	
	%Read the hardware indices of the recorded channels.
	HWChannels = data(3:3+numHWChannels-1);
	
	%Reshape the datafile so that each column represents a channel.  (Time
	%is stored as the last channel, so add one).
	data = reshape(data(3+numHWChannels:mark(end)-1), numHWChannels+1, length(data(3+numHWChannels:mark(end)-1))/(numHWChannels+1))';
	
	%Break off the last column as the time stamp for each sample.
	%time = data(:,end);
	data = data(:,1:end-1);    
elseif(trigFileFormat == -1)
    %read remainder of file.
    data = fread(fid, inf, 'float64');
    fclose(fid);    
    
 	%Parse the second float64 as startSample number.  Useful for orienting
 	%across files.
	startSampNumber = data(1);   
    
	%Parse the second float64 as the number of channels
	numHWChannels = data(2);
	
	%Read the hardware indices of the recorded channels.
	HWChannels = data(3:3+numHWChannels-1);
	
	%Reshape the datafile so that each column represents a channel.  (Time
	%is stored as the last channel, so add one).
	data = reshape(data(3+numHWChannels:end), numHWChannels+1, length(data(3+numHWChannels:end))/(numHWChannels+1))';
	
	%Break off the last column as the time stamp for each sample.
	%time = data(:,end);
	data = data(:,1:end-1);
    
elseif(trigFileFormat == -4)
    %read remainder of file.
    data = fread(fid, inf, 'float32');
    fclose(fid);    
    
 	%Parse the second float64 as startSample number.  Useful for orienting
 	%across files.
	startSampNumber = data(1);   
    
	%Parse the second float64 as the number of channels
	%numHWChannels = data(2);
	
	%Read the hardware indices of the recorded channels.
	HWChannels = 1; %data(3:3+numHWChannels-1);
	
	%Reshape the datafile so that each column represents a channel.  (Time
	%is stored as the last channel, so add one).
	%data = reshape(data(3+numHWChannels:end), numHWChannels+1, length(data(3+numHWChannels:end))/(numHWChannels+1))';
	
	%Break off the last column as the time stamp for each sample.
	%time = data(:,end);
	%data = data(:,1:end-1);
else
    error(['Unknown trigger-file format:', num2str(trigFileFormat)]);
end
    

