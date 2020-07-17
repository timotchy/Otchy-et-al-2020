function ls = plotTimeSeriesQuick(varargin)
%Auther: Aaron Andalman 2006
%The function is useful for plotting time series of uniformly spaced
%samples.  It accepts the same inputs as plot, but subsamples your
%timeseries based on the number of pixels available for display.  This
%increases the speed of plotting signicantly.

%Currently the function can only accept 1 timeseries at a time.  matricies
%of data are not supported.

%determine if first element in varargin in the axis...
v1 = varargin{1};
if(ishandle(v1))
    ax = varargin{1};
    varargin = varargin(2:end);
else
    ax = gca;
end

%find x and y in varargin
if((length(varargin)>1) && length(varargin{2}) == length(varargin{1}))
    ud.x = varargin{1};
    ud.y = varargin{2};
    ud.varargin = varargin(3:end);
else
    ud.y = varargin{1};
    ud.x = [1:length(ud.y)];
    ud.varargin = varargin(2:end);
end
ud.startndx = 1;
ud.endndx = length(ud.x);
ud.ax = ax;

set(ud.ax, 'UserData', ud);
set(ud.ax, 'ButtonDownFcn', @buttondown_plotquick);

ls = helper_plotquick(ud);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ls = helper_plotquick(ud)
%determine size of axis relative to size of x.
set(ud.ax,'Units','pixels')
pixSize = get(ud.ax,'Position');
ratio = (ud.endndx - ud.startndx) / pixSize(3);

if(ratio < 8)
    %if x is not much longer than the axis, then just plot it...  
    axes(ud.ax);
    ls = plot(ud.x(ud.startndx:ud.endndx), ud.y(ud.startndx:ud.endndx), ud.varargin{:});
    axis tight;
    set(ls, 'HitTest', 'off');
else
    %otherwise decimate signal before plotting...
    ratio = floor(ratio/4);
    
    %Old way was to simply down sample...
    %ssx = downsample(ud.x(ud.startndx:ud.endndx),ratio);
    %ssy = downsample(ud.y(ud.startndx:ud.endndx),ratio); 
    
    %Now do peak detect... truncates the remainder... would be better to append NaN
    x = ud.x(ud.startndx:ud.endndx);
    boxx = reshape(x(1:end-mod(length(x),ratio)), ratio, []);
    minx = min(boxx);
    maxx = max(boxx);
    y = ud.y(ud.startndx:ud.endndx);
    boxy = reshape(y(1:end-mod(length(y),ratio)), ratio, []);
    miny = min(boxy);
    maxy = max(boxy);
    axes(ud.ax);
    ls(1) = plot(minx, miny, ud.varargin{:}); hold on;
    ls(2) = plot(maxx, maxy, ud.varargin{:}); hold off;
    axis tight;
    set(ls, 'HitTest', 'off');
end
set(ud.ax,'Units','normalized')
set(ud.ax, 'UserData', ud);
set(ud.ax, 'ButtonDownFcn', @buttondown_plotquick);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function buttondown_plotquick(src, evnt)
ud = get(src, 'UserData');
axes(ud.ax);
mouseMode = get(gcf, 'SelectionType');
clickLocation = get(ud.ax, 'CurrentPoint');

if(strcmp(mouseMode, 'alt'))
    %right click to pan
    %NOT YET IMPLEMENTED
elseif(strcmp(mouseMode, 'open'))
    %double click to zoom out
    ud.startndx = 1;
    ud.endndx = length(ud.x);
elseif(strcmp(mouseMode, 'normal'))
    %left click to zoom in.
    rect = rbbox;
    endPoint = get(gca,'CurrentPoint'); 
    point1 = clickLocation(1,1:2);              % extract x and y
    point2 = endPoint(1,1:2);
    p1 = min(point1,point2);             % calculate locations
    offset = abs(point1-point2);         % and dimensions
    if(offset(1)/diff(xlim) < .001)
        quarter = (ud.endndx - ud.startndx) / 4;
        ud.startndx = max(1,floor(p1(1) - quarter));
        ud.endndx = min(length(ud.x),ceil(p1(1) + quarter));
    else
        ud.startndx = min(find(ud.x >= p1(1)))
        ud.endndx = max(find(ud.x <= p1(1) + offset(1)));
    end
end
set(ud.ax,'UserData',ud);
helper_plotquick(ud);

     