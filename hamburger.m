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

OooooO
oxxxxo
oxxxxo
oxxxxo
oxxxxo
0xxxx0
######

O = corner
0 = Bottom Corner
o = edge
x = regular surrounded piece

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

% MIGHT NOT NEED ALL THESE. LOOK OVER EQNS
COND = 0.5;
global DIFFUS;
DIFFUS = 1.75E-7;
PIECE_SIZE = [ ( 0.15 ./ [nx ny] ) ( 0.01 / nz ) ];
PIECE_SIZE_SQ = PIECE_SIZE .^ 2;


% Populate Variables
temp_array = zeros(nx,ny,nz+1);
temp_array(:,:,1) = Tgrill;
change_in_temp = zeros(nx,ny,nz+1);

warned_yet = false;


%PSEUDOCODE

for time=0:timestep:duration
    
    old = temp_array;

%MIDDLE:

parfor i=2:nx-1 % x
    for j=2:ny-1 % y
        for k=2:nz % z
            % nz, not nz - 1, because 1 = the grill, instead of an edge
            % piece. So while for x and y the "edge" is x = nx or y = ny,
            % the edge for z is z = nz + 1
            change_in_temp(i,j,k) = ...
                ( ...
                ( ( old(i+1,j,k) + old(i-1,j,k) - 2*old(i,j,k) ) ...
                / PIECE_SIZE_SQ(1) ) ...
                + ...
                ( ( old(i,j+1,k) + old(i,j-1,k) - 2*old(i,j,k) ) ...
                / PIECE_SIZE_SQ(2) ) ...
                + ...
                ( ( old(i,j,k+1) + old(i,j,k-1) - 2*old(i,j,k) ) ...
                / PIECE_SIZE_SQ(2) ) ...
                ) * DIFFUS;
        end % z
    end % y
end % x

                



% INNER SIDES
%     (T)OP
%     (F)RONT
%     B(A)CK
%     (L)EFT
%     (R)IGHT
% EDGES
%     TR
%     TL
%     TA
%     TF
%     BR %SPECIAL
%     BL %SPECIAL
%     BA %SPECIAL
%     BF %SPECIAL
% CORNERS
%     ALT
%     ART
%     FLT
%     FRT
%     ALB %SPECIAL
%     ARB %SPECIAL
%     FLB %SPECIAL
%     FRB %SPECIAL
    
if (~warned_yet) && abs(max(change_in_temp)) > Tgrill
        warning('VALUES ARE BLOWING UP');
        warned_yet = true;
end

temp_array = temp_array + change_in_temp;

    
end % time




end % main hamburger function

function [qpera] = radiate(Tint)

STFB = 5.60E-8; % stefan-boltzmann
CTOKEL = 273.15; %C to kelvin by adding. Kelvin to C by subtracting.
TEXT = 20 + CTOKEL;
Tint = Tint + CTOKEL;

qpera = (TEXT^4 - Tint^4) * STFB;


end % radiate function

function [q] = spatial_heat_flow(Tcurr, Tcomp, distance)
DENS = 950;
SPECHEAT = 3000;
global DIFFUS;

q = - DIFFUS * DENS * SPECHEAT * ( (Tcurr-Tcomp) / distance ) ;
end % spatial heat flow function