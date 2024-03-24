% Create a sphere
[n, m] = deal(50); % Resolution of the sphere
[x, y, z] = sphere(n);

% Make the sphere bumpy by perturbing the vertices
perturbation_strength = 0.01; % Adjust this value to increase or decrease the bumpiness
noise = perturbation_strength * randn(size(x));
bumpy_x = x + noise .* x;
bumpy_y = y + noise .* y;
bumpy_z = z + noise .* z;

% Visualize the bumpy sphere
%surf(bumpy_x, bumpy_y, bumpy_z, 'FaceColor', 'interp', 'EdgeColor', 'none');
surf(bumpy_x, bumpy_y, bumpy_z, 'CData', X, 'FaceColor', 'texturemap', 'EdgeColor', 'none');

% Improve the lighting to enhance the bumpy effect
camlight left; lighting phong; axis equal; axis off;