function [forcex,forcey]=getLiftForceO4_AR2(x,y,rsph,U,rho,L)
%% INPUT:
%  Dimensions -L < x < L, -.5L < y < .5L
%  rsph is the particle radius (rsph < .5L)
%  U is the maximum velocity in the channel
%  rho is the fluid density
%  L is the side length of the square channel
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
%  Code will return an error if a particle radius (rsph) and location
%  (x,y) are chosen such that the particle hits the wall.
%
%%

diffx = L-abs(x); diffy = .5*L-abs(y);

% if the location (x,y) is outside the channel, return an error
if x<-L || x>L || y<-L/2 || y>L/2
    error(['location outside channel: pick ',num2str(-L/2),' < x < ',num2str(L/2),' and ',num2str(-L/2),' < y < ',num2str(L/2)])
% if the particle is too big and hits the wall, return an error
elseif rsph>min(diffx,diffy)
    error(['particle is too close to wall: choose rsph < ' num2str(.5*L-max(abs(x),abs(y)))])
% otherwise proceed
else
    % if the particle radius is large, return a warning
    if rsph/L > 0.23
        warning('large particle radius, results may be inaccurate')
    end
    
    x=x/L; y=y/L; rsph=abs(rsph);

%     load HoLeal_channel_aspectratio2_Re1_mesh8_10152014.mat
    load HoLeal_channel_AR2_Re1_mesh8_08-25-15.mat

    Fx = rho*U^2*rsph^4*C4x/(L^2);
    Fy = rho*U^2*rsph^4*C4y/(L^2);
    locx=xloc; locy=yloc;

    if x<0 && y<=0
        Fx=-Fx;
        Fy=-Fy;
        locx=-locx;
        locy=-locy;
    elseif x<0 && y>0
        Fx=-Fx;
        locx=-locx;
    elseif x>=0 && y<=0
        Fy=-Fy;
        locy=-locy;
    end


    forcey = interp2(locx,locy,Fy,x,y);
    forcex = interp2(locx,locy,Fx,x,y);
end

return
