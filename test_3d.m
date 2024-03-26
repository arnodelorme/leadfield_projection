stlData = stlread('ear.stl');

vertices = stlData.Points;
faces = stlData.ConnectivityList;

% Get the minimum and maximum z-coordinates
zMin = min(vertices(:, 3));
zMax = max(vertices(:, 3));

% Create a color gradient based on z-coordinates
colorMap = [linspace(0, 1, 256)'; linspace(0, 0, 256)'; linspace(1, 0, 256)'];
textureImage = imread('ngc6543a.jpg');
cdata = textureImage;

[u, v] = meshgrid(linspace(0, 1, size(textureImage, 2)), ...
                  linspace(0, 1, size(textureImage, 1)));
u = u(:);
v = v(:);
texCoords = [u, v];

% Update the patch object with texture coordinates
set(gcf, 'Renderer', 'opengl');  % Switch to OpenGL renderer

[X,Y] = meshgrid([1:size(textureImage,2)],[1:size(textureImage,1)]);
colors = interp2(X, Y, mean(textureImage,3), vertices(:, 2), vertices(:, 3))*255;

% Normalize the z values to the range [0, 1] for color mapping
z = vertices(:,3); % Extract the z-coordinate of each vertex
z_normalized = (z - min(z)) / (max(z) - min(z));

% Generate a color gradient from blue (at the bottom) to red (at the top)
colors = [z_normalized, zeros(size(z_normalized)), 1-z_normalized];

% Create a patch object for visualization
figure;
patch('Vertices', vertices, 'Faces', faces, ...
      'FaceColor', 'interp', 'FaceVertexCData', colors, ...
      'EdgeColor', 'none');

daspect([1 1 1]);
view(3);
axis equal;
