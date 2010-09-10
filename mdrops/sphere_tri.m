function [FV] = sphere_tri(shape,maxlevel,r,winding)

% sphere_tri - generate a triangle mesh approximating a sphere
% 
% Usage: FV = sphere_tri(shape,Nrecurse,r,winding)
% 
%   shape is a string, either of the following:
%   'ico'   starts with icosahedron (most even, default)
%   'oct'   starts with octahedron
%   'tetra' starts with tetrahedron (least even)
%
%   Nrecurse is int >= 0, setting the recursions (default 0)
%
%   r is the radius of the sphere (default 1)
%
%   winding is 0 for clockwise, 1 for counterclockwise (default 0).  The
%   matlab patch command gives outward surface normals for clockwise
%   order of vertices in the faces (viewed from outside the surface).
%
%   FV has fields FV.vertices and FV.faces.  The vertices 
%   are listed in clockwise order in FV.faces, as viewed 
%   from the outside in a RHS coordinate system.
% 
% The function uses recursive subdivision.  The first
% approximation is an platonic solid, either an  icosahedron,
% octahedron or a tetrahedron.  Each level of refinement 
% subdivides each triangle face by a factor of 4 (see also 
% mesh_refine).  At each refinement, the vertices are 
% projected to the sphere surface (see sphere_project).
% 
% A recursion level of 3 or 4 is a good sphere surface, if
% gouraud shading is used for rendering.
% 
% The returned struct can be used in the patch command, eg:
% 
% % create and plot, vertices: [2562x3] and faces: [5120x3]
% FV = sphere_tri('ico',4,1);
% lighting phong; shading interp; figure;
% patch('vertices',FV.vertices,'faces',FV.faces,...
%       'facecolor',[1 0 0],'edgecolor',[.2 .2 .6]);
% axis off; camlight infinite; camproj('perspective');
% 
% See also: mesh_refine, sphere_project
%



% $Revision: 1.15 $ $Date: 2004/05/20 22:28:45 $

% Licence:  GNU GPL, no implied or express warranties
% Jon Leech (leech @ cs.unc.edu) 3/24/89
% icosahedral code added by Jim Buddenhagen (jb1556@daditz.sbc.com) 5/93
% 06/2002, adapted from c to matlab by Darren.Weber_at_radiology.ucsf.edu
% 05/2004, reorder of the faces for the 'ico' surface so they are indeed
% clockwise!  Now the surface normals are directed outward.  Also reset the
% default recursions to zero, so we can get out just the platonic solids.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~exist('shape','var') || isempty(shape),
    shape = 'ico';
end
fprintf('...creating sphere tesselation based on %s\n',shape);

% default maximum subdivision level
if ~exist('maxlevel','var') || isempty(maxlevel) || maxlevel < 0,
    maxlevel = 0;
end

% default radius
if ~exist('r','var') || isempty(r),
    r = 1;
end

if ~exist('winding','var') || isempty(winding),
    winding = 0;
end


% -----------------
% define the starting shapes

shape = lower(shape);

switch shape,
case 'tetra',
    
    % Vertices of a tetrahedron
    sqrt_3 = 0.5773502692;
    
    tetra.v = [  sqrt_3,  sqrt_3,  sqrt_3 ;   % +X, +Y, +Z  - PPP
                -sqrt_3, -sqrt_3,  sqrt_3 ;   % -X, -Y, +Z  - MMP
                -sqrt_3,  sqrt_3, -sqrt_3 ;   % -X, +Y, -Z  - MPM
                 sqrt_3, -sqrt_3, -sqrt_3 ];  % +X, -Y, -Z  - PMM
	
    % Structure describing a tetrahedron
    tetra.f = [ 1, 2, 3;
                1, 4, 2;
                3, 2, 4;
                4, 1, 3 ];
    
    FV.vertices = tetra.v;
    FV.faces    = tetra.f;
    
case 'oct',
    
    % Six equidistant points lying on the unit sphere
    oct.v = [  1,  0,  0 ;  %  X
              -1,  0,  0 ; 	% -X
               0,  1,  0 ;  %  Y
               0, -1,  0 ; 	% -Y
               0,  0,  1 ; 	%  Z
               0,  0, -1 ];	% -Z
	
    % Join vertices to create a unit octahedron
    oct.f = [ 1 5 3 ;    %  X  Z  Y  -  First the top half
              3 5 2 ;    %  Y  Z -X
              2 5 4 ;    % -X  Z -Y
              4 5 1 ;    % -Y  Z  X
              1 3 6 ;    %  X  Y -Z  -  Now the bottom half
              3 2 6 ;    %  Y  Z -Z
              2 4 6 ;    % -X  Z -Z
              4 1 6 ];   % -Y  Z -Z
    
    FV.vertices = oct.v;
    FV.faces    = oct.f;
    
case 'ico',
    
    % Twelve vertices of icosahedron on unit sphere
    tau = 0.8506508084; % t=(1+sqrt(5))/2, tau=t/sqrt(1+t^2)
    one = 0.5257311121; % one=1/sqrt(1+t^2) , unit sphere
    
    ico.v( 1,:) = [  tau,  one,    0 ]; % ZA
    ico.v( 2,:) = [ -tau,  one,    0 ]; % ZB
    ico.v( 3,:) = [ -tau, -one,    0 ]; % ZC
    ico.v( 4,:) = [  tau, -one,    0 ]; % ZD
    ico.v( 5,:) = [  one,   0 ,  tau ]; % YA
    ico.v( 6,:) = [  one,   0 , -tau ]; % YB
    ico.v( 7,:) = [ -one,   0 , -tau ]; % YC
    ico.v( 8,:) = [ -one,   0 ,  tau ]; % YD
    ico.v( 9,:) = [   0 ,  tau,  one ]; % XA
    ico.v(10,:) = [   0 , -tau,  one ]; % XB
    ico.v(11,:) = [   0 , -tau, -one ]; % XC
    ico.v(12,:) = [   0 ,  tau, -one ]; % XD
    
    % Structure for unit icosahedron
    ico.f = [  5,  8,  9 ;
               5, 10,  8 ;
               6, 12,  7 ;
               6,  7, 11 ;
               1,  4,  5 ;
               1,  6,  4 ;
               3,  2,  8 ;
               3,  7,  2 ;
               9, 12,  1 ;
               9,  2, 12 ;
              10,  4, 11 ;
              10, 11,  3 ;
               9,  1,  5 ;
              12,  6,  1 ;
               5,  4, 10 ;
               6, 11,  4 ;
               8,  2,  9 ;
               7, 12,  2 ;
               8, 10,  3 ;
               7,  3, 11 ];
	
    FV.vertices = ico.v;
    FV.faces    = ico.f;
end


% -----------------
% refine the starting shapes with subdivisions
if maxlevel,
    
    % Subdivide each starting triangle (maxlevel) times
    for level = 1:maxlevel,
        
        % Subdivide each triangle and normalize the new points thus
        % generated to lie on the surface of a sphere radius r.
        FV = mesh_refine_tri4(FV);
        FV.vertices = sphere_project(FV.vertices,r);
        
        % An alternative might be to define a min distance
        % between vertices and recurse or use fminsearch
        
    end
end

if winding,
    fprintf('...returning counterclockwise vertex order (viewed from outside)\n');
    FV.faces = FV.faces(:,[1 3 2]);
else
    fprintf('...returning clockwise vertex order (viewed from outside)\n');
end

return


function [ FV ] = mesh_refine_tri4(FV)

% mesh_refine_tri4 - creates 4 triangle from each triangle of a mesh
%
% [ FV ] = mesh_refine_tri4( FV )
%
% FV.vertices   - mesh vertices (Nx3 matrix)
% FV.faces      - faces with indices into 3 rows
%                 of FV.vertices (Mx3 matrix)
% 
% For each face, 3 new vertices are created at the 
% triangle edge midpoints.  Each face is divided into 4
% faces and returned in FV.
%
%        B
%       /\
%      /  \
%    a/____\b       Construct new triangles
%    /\    /\       [A,a,c]
%   /  \  /  \      [a,B,b]
%  /____\/____\     [c,b,C]
% A	     c	   C    [a,b,c]
% 
% It is assumed that the vertices are listed in clockwise order in
% FV.faces (A,B,C above), as viewed from the outside in a RHS coordinate
% system.
% 
% See also: mesh_refine, sphere_tri, sphere_project
% 


% ---this method is not implemented, but the idea here remains...
% This can be done until some minimal distance (D) of the mean 
% distance between vertices of all triangles is achieved.  If
% no D argument is given, the function refines the mesh once.
% Alternatively, it could be done until some minimum mean 
% area of faces is achieved.  As is, it just refines once.


% $Revision: 1.12 $ $Date: 2004/05/10 21:01:55 $

% Licence:  GNU GPL, no implied or express warranties
% History:  05/2002, Darren.Weber_at_radiology.ucsf.edu, created
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
fprintf('...refining mesh (tri4)...')

% NOTE
% The centroid is located one third of the way from each vertex to 
% the midpoint of the opposite side. Each median divides the triangle 
% into two equal areas; all the medians together divide it into six 
% equal parts, and the lines from the median point to the vertices 
% divide the whole into three equivalent triangles.

% Each input triangle with vertices labelled [A,B,C] as shown
% below will be turned into four new triangles:
%
% Make new midpoints
% a = (A+B)/2
% b = (B+C)/2
% c = (C+A)/2
%
%        B
%       /\
%      /  \
%    a/____\b       Construct new triangles
%    /\    /\       [A,a,c]
%   /  \  /  \      [a,B,b]
%  /____\/____\     [c,b,C]
% A	     c	   C    [a,b,c]
%

% Initialise a new vertices and faces matrix
Nvert = size(FV.vertices,1);
Nface = size(FV.faces,1);
V2 = zeros(Nface*3,3);
F2 = zeros(Nface*4,3);

for f = 1:Nface,
    
    % Get the triangle vertex indices
    NA = FV.faces(f,1);
    NB = FV.faces(f,2);
    NC = FV.faces(f,3);
    
    % Get the triangle vertex coordinates
    A = FV.vertices(NA,:);
    B = FV.vertices(NB,:);
    C = FV.vertices(NC,:);
    
    % Now find the midpoints between vertices
    a = (A + B) ./ 2;
    b = (B + C) ./ 2;
    c = (C + A) ./ 2;
    
    % Find the length of each median
    %A2blen = sqrt ( sum( (A - b).^2, 2 ) );
    %B2clen = sqrt ( sum( (B - c).^2, 2 ) );
    %C2alen = sqrt ( sum( (C - a).^2, 2 ) );
    
    % Store the midpoint vertices, while
    % checking if midpoint vertex already exists
    [FV, Na] = mesh_find_vertex(FV,a);
    [FV, Nb] = mesh_find_vertex(FV,b);
    [FV, Nc] = mesh_find_vertex(FV,c);
    
    % Create new faces with orig vertices plus midpoints
    F2(f*4-3,:) = [ NA, Na, Nc ];
    F2(f*4-2,:) = [ Na, NB, Nb ];
    F2(f*4-1,:) = [ Nc, Nb, NC ];
    F2(f*4-0,:) = [ Na, Nb, Nc ];
    
end

% Replace the faces matrix
FV.faces = F2;

t=toc; fprintf('done (%5.2f sec)\n',t);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FV, N] = mesh_find_vertex(FV,vertex)

    Vn = size(FV.vertices,1);
    Va = repmat(vertex,Vn,1);
    Vexist = find( FV.vertices(:,1) == Va(:,1) & ...
                   FV.vertices(:,2) == Va(:,2) & ...
                   FV.vertices(:,3) == Va(:,3) );
    if Vexist,
        if size(Vexist) == [1,1],
            N = Vexist;
        else
            msg = sprintf('replicated vertices');
            error(msg);
        end
    else
        FV.vertices(end+1,:) = vertex;
        N = size(FV.vertices,1);
    end

return

function V = sphere_project(v,r,c)

% sphere_project - project point X,Y,Z to the surface of sphere radius r
% 
% V = sphere_project(v,r,c)
% 
% Cartesian inputs:
% v is the vertex matrix, Nx3 (XYZ)
% r is the sphere radius, 1x1 (default 1)
% c is the sphere centroid, 1x3 (default 0,0,0)
%
% XYZ are converted to spherical coordinates and their radius is
% adjusted according to r, from c toward XYZ (defined with theta,phi)
% 
% V is returned as Cartesian 3D coordinates
% 

% $Revision: 1.8 $ $Date: 2004/03/29 21:15:36 $

% Licence:  GNU GPL, no implied or express warranties
% History:  06/2002, Darren.Weber_at_radiology.ucsf.edu, created
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('v','var'),
    msg = sprintf('SPHERE_PROJECT: No input vertices (X,Y,Z)\n');
    error(msg);
end

% X = v(:,1);
% Y = v(:,2);
% Z = v(:,3);

if ~exist('c','var'),
%     xo = 0;
%     yo = 0;
%     zo = 0;
    c = [0,0,0];
% else
%     xo = c(1);
%     yo = c(2);
%     zo = c(3);
end

if ~exist('r','var'), r = 1; end

for i = 1:size(v,1)
   V(i,:) = (v(i,:)./normesp(v(i,:)))*r;
end


% alternate method is to use unit vector of V
% [ n = 'magnitude(V)'; unitV = V ./ n; ]
% to change the radius, multiply the unitV
% by the radius required.  This avoids the
% use of arctan functions, which have branches.



% % Convert Cartesian X,Y,Z to spherical (radians)
% theta = atan2( (Y-yo), (X-xo) );
% phi   = atan2( sqrt( (X-xo).^2 + (Y-yo).^2 ), (Z-zo) );
% % do not calc: r = sqrt( (X-xo).^2 + (Y-yo).^2 + (Z-zo).^2);
% 
% %   Recalculate X,Y,Z for constant r, given theta & phi.
% R = ones(size(phi)) * r;
% x = R .* sin(phi) .* cos(theta);
% y = R .* sin(phi) .* sin(theta);
% z = R .* cos(phi);
% 
% V = [x y z];

return