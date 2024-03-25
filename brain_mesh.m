clear

plotTexture = true;

load refined_mesh.mat
vertices(:,2) = vertices(:,2)+20;
vertices(:,3) = vertices(:,3)+10;
pairwiseDist = ones(size(vertices,1),4);
colors = repmat(pairwiseDist(:,4)*2, [1 3]);

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
X = flipud(X);
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

% ht1 = traditionaldipfit([ 0 -20 -10 pi/4+0.1 0 pi/4+0.1 100 100 100]);
ht1 = traditionaldipfit([ 0 0 0 pi/4+0.1 0 pi/4+0.1 100 100 100]);

if 1
    ht0 = traditionaldipfit([ 0 0 0 0 0 pi/4 1 1 1]);
    coords = ht0*[ xx(:) yy(:) zz(:) ones(length(xx(:)),1)]';
    xx = coords(1,:); yy = coords(2,:); zz = coords(3,:); 
    
    tmp = xx; xx = zz; zz = tmp;
    coords = ht1*[ xx(:) yy(:) zz(:) ones(length(xx(:)),1)]';
    
    xx = reshape(coords(1,:), size(x));
    yy = reshape(coords(2,:), size(x));
    zz = reshape(coords(3,:), size(x));
    mesh = surf2patch(xx, yy, zz);
    
    % Create a figure and render the textured mesh
    colorsM = zeros(numel(X),3);
    X = X';
    for iX = 1:length(X(:))
        colorsM(iX,:) = map(X(iX),:);
        % colors{iX} = [ '#' dec2hex(floor(255*map(X(iX),1)),2) dec2hex(floor(255*map(X(iX),2)),2) dec2hex(floor(255*map(X(iX),3)),2)  ];
    end
    
    if ~plotTexture
        patch('Vertices', mesh.vertices, 'Faces', mesh.faces, 'FaceVertexCData', colorsM, 'FaceColor', 'interp', 'EdgeColor', 'none');
        %patch('Vertices', mesh.vertices, 'Faces', mesh.faces, 'FaceColor', 'none', 'EdgeColor', 'k');
        %surface(xx(1:10:end,1:10:end),yy(1:10:end,1:10:end),zz(1:10:end,1:10:end), 'FaceColor','texturemap','EdgeColor','none','Cdata',flipud(X))
        %colormap(map)
        %view(148, 90)
        view(70, 26);
        axis equal
   end
 
    % scan vertices and find the closest vertices
    diff = zeros(1,size(vertices,1));
    coords = [xx(:) yy(:) zz(:)];
    listInd = [];
    for iVert = 1:size(vertices,1)
        if vertices(iVert,1) > 20 && vertices(iVert,2) < -20 % remove quandrants to speed up
            vert = vertices(iVert,:);
            vert = vert./norm(vert)*100;
    
            % pairwise distance with xx, yy, and zz
            difftmp = bsxfun(@minus, coords, vert);
            [diff(iVert),ind] = min(vecnorm(difftmp'));
            %[diff(iVert),ind] = min(sum(abs(difftmp')));
    
            % TO DO ********************
            if diff(iVert) < 1
                colors(iVert,:) = colorsM(ind,:);
                listInd = [ listInd iVert ];
            end
        end
    end

    % vertices(listInd);
    keepInd = [];
    keepInd = zeros(1, size(faces,1), 'logical');
    parfor iFace = 1:size(faces,1)
        if any(faces(iFace,1) == listInd) || any(faces(iFace,2) == listInd) || any(faces(iFace,3) == listInd)
            % keepInd(end+1) = iFace;
            keepInd(iFace) = true;
        end
    end
    keepInd = find(keepInd);
    faces = faces(keepInd,:);

    if 1
        % plot low res mesh
        tmp = load('mesh.mat');
        verticesLR = tmp.vertices;
        verticesLR(:,2) = verticesLR(:,2)+20;
        verticesLR(:,3) = verticesLR(:,3)+10;
        pairwiseDist = ones(size(verticesLR,1),4);
        facesLR = tmp.faces;
        colorsLR = repmat(pairwiseDist(:,4)*2, [1 3]);

        % make the vertices close to the image deeper.
        for iVert = 1:size(verticesLR,1)
            if verticesLR(iVert,1) > 20 && verticesLR(iVert,2) < -20 % remove quandrants to speed up
                vert = verticesLR(iVert,:);
                vert = vert./norm(vert)*100;
    
                % pairwise distance with xx, yy, and zz
                difftmp = bsxfun(@minus, coords, vert);
                [diff(iVert),ind] = min(vecnorm(difftmp'));
                %[diff(iVert),ind] = min(sum(abs(difftmp')));
    
                % TO DO ********************
                if diff(iVert) < 5
                    verticesLR(iVert,:) = 0.9*verticesLR(iVert,:);
                end
            end
        end

        figure('position', [1440 1 1484 1237]);
        patch('Vertices', verticesLR, 'Faces', facesLR, ...
              'FaceColor', 'interp', 'FaceVertexCData', colorsLR, ...
              'EdgeColor', 'none');
    
        daspect([1 1 1]);
        view(3);
        axis equal;
        camlight headlight; lighting phong; axis equal; axis off;
    end
    
    if plotTexture
        patch('Vertices', vertices, 'Faces', faces, ...
              'FaceColor', 'interp', 'FaceVertexCData', colors, ...
              'EdgeColor', 'none');
    
        daspect([1 1 1]);
        view(70, 26);
        axis equal;
        camlight right; lighting phong; axis equal; axis off;
    end

else
    %% create sphere
    [k,l,m] = sphere;
    tmp = k; k = m; m = tmp;
    coords = ht1*[ k(:) l(:) m(:) ones(length(k(:)),1)]';
    surf(reshape(coords(1,:), size(k)), reshape(coords(2,:), size(k)), reshape(coords(3,:), size(k)), 'facecolor', 'none');
end

% next steps
% - refine the vertices and faces
%    - remove vertices which do not have a face
%    - reindex the faces
% - transform back to original space (y-20 and z-10) and save
% - use mesh to compute leadfield on 370 electrodes
% - propagate R G and B separately, add the 3 to get the color
% - add more electrodes for higher scalp resolution
