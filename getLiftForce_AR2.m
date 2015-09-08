function [forcex,forcey]=getLiftForce_AR2(x,y,rsph,U, rho, L)
%% INPUT:
%  Dimensions -L < x < L, -.5L < y < .5L
%  rsph is the particle radius (rsph < .5L)
%  U is the maximum velocity in the channel
%  rho is the fluid density
%  L is the shortest side length of the channel
%
%% OUTPUT
%  forcex - force in the x direction
%  forcey - force in the y direction
%
%% DEPENDENCIES
%  Need to have HoLeal_channel_AR2_Re1_mesh8_08-25-15.mat in the Matlab path
%  Uses Matlab's interp2 for 2D interpolation
%
%% NOTES
%  Results may be inaccurate for large particle radius
%  i.e. rsph > .2L
%
%  Code will return an error if particle radius (rsph) and location
%  (x,y) are chosen such that the particle hits the wall. 
%
%%

diffx = L-abs(x); diffy = .5*L-abs(y);

% if the location (x,y) is outside the channel, return an error
if x<-L || x>L || y<-L/2 || y>L/2
    error(['location outside channel: pick ',num2str(-L),' < x < ',...
        num2str(L),' and ',num2str(-L/2),' < y < ',num2str(L/2)])
% if the particle is too big and hits the wall, return an error
elseif rsph>min(diffx,diffy)
    error('particle is too close to wall')
% otherwise proceed
else
    % if the particle radius is large, return a warning
    if rsph/L > 0.23
        warning('large particle radius, results may be inaccurate')
    end
    
    x0=x/L; y0=y/L; rsph=abs(rsph);
   
%     load HoLeal_channel_aspectratio2_Re1_mesh8_12032014.mat
%     load HoLeal_channel_AR2_Re1_mesh8_05-07-15.mat
    load HoLeal_channel_AR2_Re1_mesh8_08-25-15.mat

    
    Fx = C4x + (rsph/L)*C5x;
    Fy = C4y + (rsph/L)*C5y;

    if x0<0 && y0<=0
        Fx=-Fx;
        Fy=-Fy;
        xloc=-xloc;
        yloc=-yloc;
    elseif x0<0 && y0>0
        Fx=-Fx;
        xloc=-xloc;
    elseif x0>=0 && y0<=0
        Fy=-Fy;
        yloc=-yloc;
    end

    forcex = (rho*U^2*rsph^4/(L^2))*interp2(xloc,yloc,Fx,x0,y0);
    forcey = (rho*U^2*rsph^4/(L^2))*interp2(xloc,yloc,Fy,x0,y0);
    
end

return
