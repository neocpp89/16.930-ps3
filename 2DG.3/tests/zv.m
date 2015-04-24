function [fx,fy] = zv(u,p,param,time)
%CONVECTIONV Calculate the volume flux for the Linear Convection equation.
%   [FX,FY]=CONVECTIONV(U,P,PARAM,TIME)
%
%      U(np):      np left (or plus) states
%      P(np,2):    np x,y coordinates
%      PARAM{1}:   Cell array containing either 
%                  - a constant velocity field [u,v] = PARAM{1}
%                  - a pointer to a funvtion that returns of velocity
%                    a velocity field as a function of P vvec = PARAM{1}(p)
%      TIME:       Not used
%      FX(np):     np fluxes in the x direction (f plus)  
%      FY(np):     np fluxes in the y direction (f plus)  
%                          
% - Written by: J. Peraire
%
fx = zeros(size(u));
fy = zeros(size(u));

