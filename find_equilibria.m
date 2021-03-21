%% Find equilibria of system by solving RHS equal to zero 
syms x y z
eqn1 = x*(s(1)-x-ro(1,2)*y-ro(1,3)*z) == 0;
eqn2 = y*(s(2)-y-ro(2,1)*x-ro(2,3)*z) == 0;
eqn3 = z*(s(3)-z-ro(3,1)*x-ro(3,2)*y) == 0;
sol = solve([eqn1, eqn2, eqn3], [x, y, z]);
mysol = [sol.x(:),sol.y(:),sol.z(:)]; % matrix of equilibria
 
%% find Jacobian of system for each positive equilibria
syms x y z
J=jacobian([x*(s(1)-x-ro(1,2)*y-ro(1,3)*z), y*(s(2)-y-ro(2,1)*x-ro(2,3)*z), z*(s(3)-z-ro(3,1)*x-ro(3,2)*y)], [x, y, z]);
J1 = subs(J, {x,y,z}, {0,0,0});
J2 = subs(J, {x,y,z}, {1.1,0,0});
J3 = subs(J, {x,y,z}, {0,1,0});
J4 = subs(J, {x,y,z}, {0,0,0.9});
J5 = subs(J, {x,y,z}, {4/15,7/15,4/15});

% find eigenvalues and eigenvectors of equilibria 
[V,D] = eig(J1);
[V1,D1] = eig(J2);
[V2,D2] = eig(J3);
[V3,D3] = eig(J4);
[V4,D4] = eig(J5);
