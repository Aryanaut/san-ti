syms u1 u2 u3 u4 mu mus

r1=[ u1+mu; u2];
r2=[u1-mus; u2];
d1 = sqrt(r1(1)^2+r1(2)^2)^3;  % Distance to the earth
d2 = sqrt(r2(1)^2+r2(2)^2)^3;  % Distance to the moon

F = [u3;
     u4;
     2*u4 + u1 - mus*(u1+mu)/d1 - mu*(u1-mus)/d2;
    -2*u3 + u2 - mus*u2/d1 - mu*u2/d2;];

simplify(jacobian(F,[u1,u2,u3,u4]))


%output
% can be further simplified, e.g. completing square
% the expression used in the code is after manipulation
% [                                                                                                                                                                                                      0,                                                                                                                                                       0,  1, 0]
% [                                                                                                                                                                                                      0,                                                                                                                                                       0,  0, 1]
% [(3*mu*(2*mus - 2*u1)*(mus - u1))/(2*((mus - u1)^2 + u2^2)^(5/2)) - mu/((mus - u1)^2 + u2^2)^(3/2) - mus/((mu + u1)^2 + u2^2)^(3/2) + (3*mus*(2*mu + 2*u1)*(mu + u1))/(2*((mu + u1)^2 + u2^2)^(5/2)) + 1,                                                      (3*mus*u2*(mu + u1))/((mu + u1)^2 + u2^2)^(5/2) - (3*mu*u2*(mus - u1))/((mus - u1)^2 + u2^2)^(5/2),  0, 2]
% [                                                                                             (3*mus*u2*(mu + u1))/((mu + u1)^2 + u2^2)^(5/2) - (3*mu*u2*(2*mus - 2*u1))/(2*((mus - u1)^2 + u2^2)^(5/2)), (3*mu*u2^2)/((mus - u1)^2 + u2^2)^(5/2) - mu/((mus - u1)^2 + u2^2)^(3/2) - mus/((mu + u1)^2 + u2^2)^(3/2) + (3*mus*u2^2)/((mu + u1)^2 + u2^2)^(5/2) + 1, -2, 0]
% 
