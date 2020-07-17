p = p;
q = q;

if length(p) ~= length(q)
    return
end

songs = length(p);

for i = 1:songs
    %Create path vector
    path = [p{i};q{i}]; 
    %Get the lengths of master, slave, and path vectors
    s_length = max(path(2,:));
    m_length = max(path(1,:));
    path_length = sqrt(m_length^2+s_length^2);
    %Find angle to rotate vector so it finishes horizontal
    theta = -1*atan(s_length/m_length);
    %Transform vector
    trans = [cos(theta),-1*sin(theta);sin(theta), cos(theta)];
    rot_path = trans*path;
    %Linearly scale the "x" component to the length of master
    sRot_path = rot_path(1,:)*(m_length/path_length);
    %Interpolate to find evenly spaced time points
    stretch(i,:) = interp1(sRot_path(1,:),rot_path(2,:),1:floor(m_length),'linear');
end

%Take the time derivative to find areas of change
diff_stretch = diff(stretch')';

%Plot scaled image matrix
imagesc(mat2gray(diff_stretch));