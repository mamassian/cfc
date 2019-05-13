function ABC = bivariate_any_triangle(vertices)

% Compute the integral of the independent (rho = 0) bivariate Normal
%   distribution with parameters (mx = 0, my = 0, sigma = 1)
%   over the arbitrary triangle ABC bounded by vertices A(Ax, Ay),
%   B(Bx, By) and C(Cx, Cy).
%
% Input: vertices = [Ax, Ay; Bx, By; Cx, Cy]
%
% Convention: the integral is positive if the triangle is positively
%   oriented, i.e. if its area is on the left as one travels along the
%   contour (i.e. the vertices are visited in a counter-clockwise fashion). 
%
% Uses the associated 'bivariate_simple_triangle' function.
% 
% 10-AUG-2017 - created - pascal mamassian
% 21-FEB-2019 - pm: added special case where all vertices are identical


Ax = vertices(1,1);  Ay = vertices(1,2);
Bx = vertices(2,1);  By = vertices(2,2);
Cx = vertices(3,1);  Cy = vertices(3,2);

% -> if all three points are identical, return 0
% if ((Bx == Ax) && (Cx == Ax) && (By == Ay) && (Cy == Ay))
% -> if at least two points are identical, return 0
if ((Bx == Ax) && (By == Ay)) || ((Cx == Ax) && (Cy == Ay)) || ((Bx == Cx) && (By == Cy))
    ABC = 0;
else

% -> let D be the orthogonal projection of the origin onto BC
tmp = (Bx * Cy - By * Cx) / ((Bx - Cx)^2 + (By - Cy)^2);
Dx = - tmp * (By - Cy);
Dy =   tmp * (Bx - Cx);

% -> let E be the projection of O onto AB
tmp = (Ax * By - Ay * Bx) / ((Ax - Bx)^2 + (Ay - By)^2);
Ex = - tmp * (Ay - By);
Ey =   tmp * (Ax - Bx);

% -> let F be the projection of O onto AC
tmp = (Ax * Cy - Ay * Cx) / ((Ax - Cx)^2 + (Ay - Cy)^2);
Fx = - tmp * (Ay - Cy);
Fy =   tmp * (Ax - Cx);


% -> the area under ABC is the same as that under OAB + OBC + OCA
%   introducing right triangles, we have
%   OAB = OAE + OEB
%   OBC = OBD + ODC
%   OCA = OCF + OFA
%   so ABC = OAE + OEB + OBD + ODC + OCF + OFA

OEA = signed_area(Ex, Ey, Ax, Ay);
OEB = signed_area(Ex, Ey, Bx, By);
ODB = signed_area(Dx, Dy, Bx, By);
ODC = signed_area(Dx, Dy, Cx, Cy);
OFC = signed_area(Fx, Fy, Cx, Cy);
OFA = signed_area(Fx, Fy, Ax, Ay);

OAE = - OEA;
OBD = - ODB;
OCF = - OFC;

ABC = OAE + OEB + OBD + ODC + OCF + OFA;

end

% -> nested function to compute the signed area within the right triangle
%   OMN, where O is the origin, and M is right-angle vertex
    function SA = signed_area(Mx, My, Nx, Ny)
        hh = sqrt(Mx*Mx + My*My);
        kk = sqrt((Nx-Mx)^2 + (Ny-My)^2);
        SA = bivariate_simple_triangle(hh, kk);
        
        % -> check that the triangle is positively oriented
        %   by computing the cross-product of NO by OM
        pos_orient = sign(Mx*Ny - My*Nx);
        SA = SA * pos_orient;
    end

end
% *** THE END ***