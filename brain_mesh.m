clear

load mesh.mat

figure('position', [1440 1 1484 1237]);
pairwiseDist = ones(size(vertices,1),4);
colors = pairwiseDist(:,4)*2;
patch('Vertices', vertices, 'Faces', faces, ...
      'FaceColor', 'interp', 'FaceVertexCData', colors, ...
      'EdgeColor', 'none');

daspect([1 1 1]);
view(3);
axis equal;
colorbar;

camlight left; camlight right; lighting phong; axis equal; axis off;

%% make a plot shere
if 0
    hold on
    [x,y,z] = sphere;
    surface(x*100,y*100-20,z*100-10, 'FaceColor','none','EdgeColor','k')
    view(112, 4)
    axis equal
    axis off
end

%%
hold on;
load mandrill
clear xx yy zz;
X = X(1:4:end,1:4:end);
x = linspace(-0.25,0.25, size(X,1));
y = linspace(-0.25,0.25, size(X,2));
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

    % project back to 
    % x, y and z
    % top of sphere 0,0,1
end

[xx,yy,zz] = sphere;

ht1 = traditionaldipfit([ 0 -20 -10     0    0  0  100 100 100]);
%ht2 = traditionaldipfit([ 0  0     0 -pi/4 -pi/4 0  1   1   1]);
coords = ht1*[ xx(:) yy(:) zz(:) ones(length(xx),1)]';
coords = coords';

xx = reshape(coords(:,1), size(x));
yy = reshape(coords(:,2), size(x));
zz = reshape(coords(:,3), size(x));

% xx = reshape(xx, size(x))*100;
% yy = reshape(yy, size(x))*100-20;
% zz = reshape(zz, size(x))*100-10;

mesh = surf2patch(xx, yy, zz);

% Create a figure and render the textured mesh
colors = zeros(numel(X),3);
X = X';
for iX = 1:length(X(:))
    colors(iX,:) = map(X(iX),:);
    % colors{iX} = [ '#' dec2hex(floor(255*map(X(iX),1)),2) dec2hex(floor(255*map(X(iX),2)),2) dec2hex(floor(255*map(X(iX),3)),2)  ];
end

%patch('Vertices', mesh.vertices, 'Faces', mesh.faces, 'FaceVertexCData', colors, 'FaceColor', 'interp', 'EdgeColor', 'none');
patch('Vertices', mesh.vertices, 'Faces', mesh.faces, 'FaceColor', 'none', 'EdgeColor', 'k');

%surface(xx(1:10:end,1:10:end),yy(1:10:end,1:10:end),zz(1:10:end,1:10:end), 'FaceColor','texturemap','EdgeColor','none','Cdata',flipud(X))
%colormap(map)
%view(148, 90)
view(3);
axis equal


