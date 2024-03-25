clear
[x,y,z] = sphere;

angle_that_create_issue = pi/4;
ht1 = traditionaldipfit([ 0 -20 -10     -pi/4    angle_that_create_issue  pi/4  1 1 1]);
coords = ht1*[ x(:) y(:) z(:) ones(length(x(:)),1)]';

figure; 
surf(reshape(coords(1,:), size(x)), reshape(coords(2,:), size(x)), reshape(coords(3,:), size(x)), 'facecolor', 'none');
axis equal
view(21,36);
