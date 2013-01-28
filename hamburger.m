function [temp_array]  =  hamburger(nx, ny, nz, Tgrill, duration, timestep)
% 
% DESCRIPTION
%     Simulate a hamburger being cooked from a frozen temperature.
% 
% PARAMETERS
% 
%     nx                  is the number of sections in the
%                             X-direction 
% 
%     ny                  is the number of sections in the
%                             Y-direction 
% 
%     nz                  is the number of sections in the
%                             Z-direction 
% 
%     Tgrill              is the temperature of the grill
%                                 (Celsius) 
% 
%     duration            is the time the hamburger is on the grill
%                                 (seconds)
% 
%     timestep            is the timestep used in calculations
%                                 (seconds)
%        
%                             
% RETURNS
% 
%     temp_array          is an array of nx-by-ny-by-nz elements,
%                             with the final temperature
%                             for each location in the patty
% 
% AUTHOR
%     Kevin Granger <kmg2728@rit.edu>
%     2013-02-01

%{
Additional Documentation:
%}

% Check starting values :
%Physical
assert( Tgrill >= 0, ...
    'The grill should probably be heated up *hotter* than the frozen meat. Please use a temperature above zero degrees celcius.');
assert( min([nx ny nz]) > 0 && max(rem([nx ny nz],1)) == 0,...
    'You must give a positive, nonzero, integer number of sections for each dimension of the burger.')
%Temporal
assert(duration >= timestep,...
    'Your duration must be greater than your timestep.');
assert(duration > 0 && duration > timestep,...
    'Your duration must be positive and nonzero.');
assert(timestep > 0, 'Your timestep must be positive and nonzero.');
    
% Constants
STFB = ; % stefan-boltzmann
CTOKEL = 273.15; % convert c to kelvin by adding
KELTOC = -273.15; % convert kelvin to c by adding

% MIGHT NOT NEED ALL THESE. LOOK OVER EQNS
DIFFUS = 
COND = 
SPECHEAT = 
DXSQRD = 
DENS = 
MASS =

% Populate Variables
temp_array = zeros(nx,ny,nz+1);
temp_array(:,:,1) = Tgrill;


PSEUDOCODE

for time
DO MIDDLE PARTS (incl inner bottom) (triple nested forloop)
% REMEMBER: (B)OTTOM IS CONDUCTIVE

INNER SIDES
    (T)OP
    (F)RONT
    B(A)CK
    (L)EFT
    (R)IGHT
EDGES
    TR
    TL
    TA
    TF
    BR %SPECIAL
    BL %SPECIAL
    BA %SPECIAL
    BF %SPECIAL
CORNERS
    ALT
    ART
    FLT
    FRT
    ALB %SPECIAL
    ARB %SPECIAL
    FLB %SPECIAL
    FRB %SPECIAL
    
    
end %for




end % function