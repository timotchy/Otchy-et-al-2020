function time = getFileTime(filename)

strparts = regexp(filename,'_', 'split');

y = str2double(strparts{3});
m = str2double(strparts{4});
d = str2double(strparts{5});
th = str2double(strparts{6});
tm = str2double(strparts{7});
ts = strparts{8};
ts =  str2double(ts(1:end-4));

time = datenum(y,m,d,th,tm,ts);
