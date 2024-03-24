%% Basic solution
load mandrill
padding = 100;

% pad all sides
X = [ones(size(X,1),padding) X zeros(size(X,1), padding) ]*2;
X = [ones(padding, size(X,2)); X; zeros(padding, size(X,2)) ]*2;
figure;
[x,y,z] = sphere;
Xhalf = [ones(size(X,1),375)*max(max(X)), X, ones(size(X,1),125)*max(max(X))]*2;
Xhalf(Xhalf(:) == 0) = 2;
surface(x,y,z, 'FaceColor','texturemap','EdgeColor','none','Cdata',flipud(Xhalf))
colormap(map)
view(112, 4)
axis equal
axis off

%% project onto sphere
clear xx yy zz;
load mandrill
padding = 100;

% pad all sides
X = [ones(size(X,1),padding) X zeros(size(X,1), padding) ]*2;
X = [ones(padding, size(X,2)); X; zeros(padding, size(X,2)) ]*2;

x = [-0.5 0 0.5];
y = [-0.5 0 0.5];
[x,y] = meshgrid(x,y);
z = ones(3,3);
x, y, z

for k = 1:length(x(:))
    dist = abs(x(k)+j*y(k));
    ang  = angle(x(k)+j*y(k));

    dplane = sin(dist);
    zz(k)  = cos(dist);

    realpos = dplane*exp(j*ang);
    xx(k) = real(realpos);
    yy(k) = imag(realpos);

    % project back to x, y and z
    % top of sphere 0,0,1
end
xx = reshape(xx, size(x))
yy = reshape(yy, size(x))
zz = reshape(zz, size(x))

figure;
surface(xx,yy,zz, 'FaceColor','texturemap','EdgeColor','none','Cdata',flipud(X))
colormap(map)
view(112, 4)
axis equal

%% Same but use the image coordinates
load mandrill
padding = 0;
X = [ones(size(X,1),padding) X zeros(size(X,1), padding) ]*2;
X = [ones(padding, size(X,2)); X; zeros(padding, size(X,2)) ]*2;

% Create a basic 3D surface for the example (e.g., a sphere)
[u, v] = meshgrid(linspace(0, 2 * pi, size(X, 2)), linspace(0, pi, size(X, 1)));
x = cos(u) .* sin(v);
y = sin(u) .* sin(v);
z = cos(v);

% Display the surface and apply the texture
surf(x, y, z, 'CData', X, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
axis equal  % Adjust the axis for a better view
colormap(map)


%% preserve geometry
load mandrill
clear xx yy zz;
x = linspace(-1,1, size(X,1));
y = linspace(-1,1, size(X,2));
[x,y] = meshgrid(x,y);
z = ones(size(x));

for k = 1:length(x(:))
    dist = abs(x(k)+j*y(k));
    ang  = angle(x(k)+j*y(k));

    dplane = sin(dist);
    zz(k)  = cos(dist);

    realpos = dplane*exp(j*ang);
    xx(k) = real(realpos);
    yy(k) = imag(realpos);

    % project back to x, y and z
    % top of sphere 0,0,1
end
xx = reshape(xx, size(x));
yy = reshape(yy, size(x));
zz = reshape(zz, size(x));

figure;
surface(xx(1:10:end,1:10:end),yy(1:10:end,1:10:end),zz(1:10:end,1:10:end), 'FaceColor','texturemap','EdgeColor','none','Cdata',flipud(X))
colormap(map)
view(112, 4)
axis equal

%% using the patch function
% Load the texture image
load mandrill % X contains the data
X = X(1:4:end, 1:4:end);
[u, v] = meshgrid(linspace(0, 2 * pi, size(X, 2)), linspace(0, pi, size(X, 1)));
x = cos(u) .* sin(v);
y = sin(u) .* sin(v);
z = cos(v);
mesh = surf2patch(x, y, z);

% Create a figure and render the textured mesh
colors = zeros(numel(X),3);
for iX = 1:length(X(:))
    colors(iX,:) = map(X(iX),:);
end

figure;
patch('Vertices', mesh.vertices, 'Faces', mesh.faces, 'FaceVertexCData', X(:), 'FaceColor', 'interp', 'EdgeColor', 'none');
%patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'interp', 'EdgeColor', 'none');
view(3);
axis equal;
colormap(map)

%% preserve geometry and use the patch function
load mandrill
X = X(1:4:end, 1:4:end);
clear xx yy zz;
x = linspace(-1,1, size(X,1));
y = linspace(-1,1, size(X,2));
[x,y] = meshgrid(x,y);
z = ones(size(x));

for k = 1:length(x(:))
    dist = abs(x(k)+j*y(k));
    ang  = angle(x(k)+j*y(k));

    dplane = sin(dist);
    zz(k)  = cos(dist);

    realpos = dplane*exp(j*ang);
    xx(k) = real(realpos);
    yy(k) = imag(realpos);

    % project back to x, y and z
    % top of sphere 0,0,1
end
xx = reshape(xx, size(x));
yy = reshape(yy, size(x));
zz = reshape(zz, size(x));

mesh = surf2patch(xx, zz, yy);

% Create a figure and render the textured mesh
X = X';
colors = zeros(numel(X),3);
for iX = 1:length(X(:))
    colors(iX,:) = map(X(iX),:);
end

figure;
patch('Vertices', mesh.vertices, 'Faces', mesh.faces, 'FaceVertexCData', X(:), 'FaceColor', 'interp', 'EdgeColor', 'none');
%patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'interp', 'EdgeColor', 'none');
view(3);
axis equal;
colormap(map)