% SphereSurfaceRD.m
%
% Finds solutions to the Grey Scott reaction-diffusion model being on the surface of 
% a sphere using a vectorized surface finite element method.
%
% Copyright (C) 2015 Michael N
%
% This program is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of  MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.

% Reset for initial state
clear
close all

% Number of time steps
nt = 7500;
dt = 1.0;

% Load premade sphere mesh from file
verts = dlmread('verticies_high_res.txt'); % Read verticies
t = dlmread('indicies_high_res.txt'); % Read connection matrix
x = verts(:,1);
y = verts(:,2);
z = verts(:,3);
p = [x y z];
Np = length(p(:,1));

p = p'; 
t = t'; 

% Gray-Scott RD
F  = 0.036;
k  = 0.057;
Du = 0.02;
Dv = 0.01;

Perturbation = 0.01;
upde = (ones(Np, 1)*(1.00-Perturbation) + (Perturbation)*rand(Np, 1));
vpde = (ones(Np, 1)*(0.00)              + (Perturbation)*rand(Np, 1));

% Initial conditions on the centre square to get reaction going
BoxSize = 0.05;

% IC one
%upde((-BoxSize < p(1,:)) & (p(1,:) < BoxSize) & (p(3,:)>0)) = 0.5;
%vpde((-BoxSize < p(1,:)) & (p(1,:) < BoxSize) & (p(3,:)>0)) = 0.25;

% IC two
upde((-BoxSize < p(1,:)) & (p(1,:) < BoxSize) & (-BoxSize < p(2,:)) & (p(2,:) < BoxSize)) = 0.5;
vpde((-BoxSize < p(1,:)) & (p(1,:) < BoxSize) & (-BoxSize < p(2,:)) & (p(2,:) < BoxSize)) = 0.25;

% Get system setup
[ LHSu, LHSv, LLu, UUu, LLv, UUv, setup, TM ] = getRDSetup( p, t, Du, Dv, dt, Np );

% Main simulation loop
for i = 1:nt
    % Set reaction load vectors
    Fu = -upde.*vpde.*vpde + F*(1-upde);
    Fv =  upde.*vpde.*vpde - (F+k)*vpde;
    
    % Solve RD equation
    [upde, vpde] = StepRD( upde, vpde, Fu, Fv, dt, LHSu, LHSv, LLu, UUu, LLv, UUv, TM );
    
    % Render and write to GIF file
    if mod(i, 25) == 1
        figure(101)
        set(gcf, 'Position', get(0,'ScreenSize'));
        set(gcf, 'Renderer', 'zbuffer'); % Needed for win 7 gifs for Matlab-related reasons.
        %clf
        %trimesh(NN,u,v,upde);...   
        h = trisurf(t', x', y', z', vpde, 'FaceColor','interp', 'EdgeColor','none');

        lighting gouraud
        shading interp
        lightangle(-45,30)
        
        h.FaceLighting = 'gouraud';
        h.AmbientStrength = 0.94;
        h.DiffuseStrength = 0.9;
        h.SpecularStrength = 0.01;
        h.SpecularExponent = 2;
        h.BackFaceLighting = 'unlit';

        axis equal
        axis off
        grid off

        axis([-1.4 1.4 -1.4 1.4 -1.4 1.4])
        view(60,15);

        %caxis([0.0 0.4]);
        %colormap('HSV');
        colormap('jet');
        %view(0,90)
        drawnow
        
        % Uncomment to save to file
       %filename = 'two.gif';
       %im = getframe;
       %[imind,cm] = rgb2ind(im.cdata, 256);
       %if i == 1
       %    imwrite(imind, cm , filename, 'gif', 'Loopcount', inf);
       %else
       %    imwrite(imind, cm , filename, 'gif', 'DelayTime', 0, 'WriteMode', 'append');
       %end

    end
end

%
%  Function used to set up the reaction diffusion solver for the reset of
%  the simulation.
%
function [ LHSu, LHSv, LLu, UUu, LLv, UUv, setup, TM ] = getRDSetup( p, t, Du, Dv, dt, Np )
    p = p';
    t = t';

    % Local stiffness matrix 
    Ke = [ ...
        1,   -1/2,  -1/2; ...
        -1/2, 1/2,   0;   ...
        -1/2,   0,   1/2 ];
    
    % Local mass matrix 
    Me = 1/24 * [...
        2, 1, 1; ...
        1, 2, 1; ...
        1, 1, 2 ];

    % Get areas of all triangles
    v1 = p( t(:,1) , : );
    v2 = p( t(:,2) , : );
    v3 = p( t(:,3) , : );
    
    % using cross product
    Areas = sqrt( sum( cross(v2-v1, v1-v3).^2, 2 ) ) / 2;
    
    % Set up LHS and RHS matricies
    LHSu = sparse(t(:,1), t(:,1), Me(1,1)*Areas +dt*Du*Ke(1,1)*Areas, Np, Np)+ ...
         sparse(t(:,1), t(:,2), Me(1,2)*Areas +dt*Du*Ke(1,2)*Areas, Np, Np)+ ...
         sparse(t(:,1), t(:,3), Me(1,3)*Areas +dt*Du*Ke(1,3)*Areas, Np, Np)+ ...
         sparse(t(:,2), t(:,1), Me(2,1)*Areas +dt*Du*Ke(2,1)*Areas, Np, Np)+ ...
         sparse(t(:,2), t(:,2), Me(2,2)*Areas +dt*Du*Ke(2,2)*Areas, Np, Np)+ ...
         sparse(t(:,2), t(:,3), Me(2,3)*Areas +dt*Du*Ke(2,3)*Areas, Np, Np)+ ...
         sparse(t(:,3), t(:,1), Me(3,1)*Areas +dt*Du*Ke(3,1)*Areas, Np, Np)+ ...
         sparse(t(:,3), t(:,2), Me(3,2)*Areas +dt*Du*Ke(3,2)*Areas, Np, Np)+ ...
         sparse(t(:,3), t(:,3), Me(3,3)*Areas +dt*Du*Ke(3,3)*Areas, Np, Np);
     
     LHSv = sparse(t(:,1), t(:,1), Me(1,1)*Areas +dt*Dv*Ke(1,1)*Areas, Np, Np)+ ...
         sparse(t(:,1), t(:,2), Me(1,2)*Areas +dt*Dv*Ke(1,2)*Areas, Np, Np)+ ...
         sparse(t(:,1), t(:,3), Me(1,3)*Areas +dt*Dv*Ke(1,3)*Areas, Np, Np)+ ...
         sparse(t(:,2), t(:,1), Me(2,1)*Areas +dt*Dv*Ke(2,1)*Areas, Np, Np)+ ...
         sparse(t(:,2), t(:,2), Me(2,2)*Areas +dt*Dv*Ke(2,2)*Areas, Np, Np)+ ...
         sparse(t(:,2), t(:,3), Me(2,3)*Areas +dt*Dv*Ke(2,3)*Areas, Np, Np)+ ...
         sparse(t(:,3), t(:,1), Me(3,1)*Areas +dt*Dv*Ke(3,1)*Areas, Np, Np)+ ...
         sparse(t(:,3), t(:,2), Me(3,2)*Areas +dt*Dv*Ke(3,2)*Areas, Np, Np)+ ...
         sparse(t(:,3), t(:,3), Me(3,3)*Areas +dt*Dv*Ke(3,3)*Areas, Np, Np);
        
    % Setup the linear system solver
    setup.type     = 'ilutp';
    setup.droptol  = 1e-3;
    setup.milu     = 'off';
    setup.udiag    = 1;
    setup.thresh   = 0;
    
    % Use incomplete LU decomposition as the solver preconditioner
    [LLu, UUu]       = ilu(LHSu, setup);
    [LLv, UUv]       = ilu(LHSv, setup);

    % Generate the reaction diffusion linear system
    TM = sparse(t(:,1), t(:,1), Me(1,1)*Areas, Np, Np)+ ...
         sparse(t(:,1), t(:,2), Me(1,2)*Areas, Np, Np)+ ...
         sparse(t(:,1), t(:,3), Me(1,3)*Areas, Np, Np)+ ...
         sparse(t(:,2), t(:,1), Me(2,1)*Areas, Np, Np)+ ...
         sparse(t(:,2), t(:,2), Me(2,2)*Areas, Np, Np)+ ...
         sparse(t(:,2), t(:,3), Me(2,3)*Areas, Np, Np)+ ...
         sparse(t(:,3), t(:,1), Me(3,1)*Areas, Np, Np)+ ...
         sparse(t(:,3), t(:,2), Me(3,2)*Areas, Np, Np)+ ...
         sparse(t(:,3), t(:,3), Me(3,3)*Areas, Np, Np);
end

%
% Function to solve reaction diffusion system at each subsequent time step.
%
function [ u, v ] = StepRD( uprev, vprev, Fu, Fv, dt, LHSu, LHSv, LLu, UUu, LLv, UUv, TM )
    % Generate left hand side
    Mlu = dt*Fu + uprev;
    Mlv = dt*Fv + vprev;

    % and right hand side of the system
    RHSu = TM*Mlu;
    RHSv = TM*Mlv;
    
    % Set tolerances for iterative solver
    tol            = 1e-3;
    maxiter        = 200;
    
    % Solve linear system using bi congugate gradient preconditioned with 
    % an incomplete LU decomposition
    [u, ~, ~, ~] = bicgstab(LHSu, RHSu, tol, maxiter, LLu, UUu, uprev);
    [v, ~, ~, ~] = bicgstab(LHSv, RHSv, tol, maxiter, LLv, UUv, vprev);
end
