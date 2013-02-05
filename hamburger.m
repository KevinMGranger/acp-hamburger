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

% Physical
assert( Tgrill >= 0, ...
    'The grill should probably be heated up *hotter* than the frozen meat. Please use a temperature above zero degrees celcius.');
assert( min([nx ny nz]) > 0 && max(rem([nx ny nz],1)) == 0,...
    'You must give a positive, nonzero, integer number of sections for each dimension of the burger.')
% Temporal
assert(duration >= timestep,...
    'Your duration must be greater than your timestep.');
assert(duration > 0 && duration > timestep,...
    'Your duration must be positive and nonzero.');
assert(timestep > 0, 'Your timestep must be positive and nonzero.');
   

% Constants

COND = 0.5;
global DENS;
DENS = 950;
global SPECHEAT;
SPECHEAT = 3000;
global DIFFUS;
DIFFUS = 1.75E-7;

PIECE_SIZE = [ ( 0.15 ./ [nx ny] ) ( 0.01 / nz ) ]; % also serves as dist.
PIECE_SIZE_SQ = PIECE_SIZE .^ 2;
PIECE_MASS = DENS * prod(PIECE_SIZE);


% Populate Working Variables
temp_array = zeros(nx,ny,nz+1);
temp_array(:,:,1) = Tgrill;
change_in_temp = zeros(nx,ny,nz+1);

warned_yet = false;



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
                ( ( old(i+1,j,k) + old(i-1,j,k) - 2*old(i,j,k) ) ...
                / PIECE_SIZE_SQ(1) ) ...
                + ...
                ( ( old(i,j+1,k) + old(i,j-1,k) - 2*old(i,j,k) ) ...
                / PIECE_SIZE_SQ(2) ) ...
                + ...
                ( ( old(i,j,k+1) + old(i,j,k-1) - 2*old(i,j,k) ) ...
                / PIECE_SIZE_SQ(3) );
        end % z
    end % y
end % x

change_in_temp(2:end-1,2:end-1,2:end-1) = ...
    change_in_temp(2:end-1,2:end-1,2:end-1) .* DIFFUS;
                


% INNER SIDES
% For each piece, take spatial heat flow from 5 meat-touching sides, then
% the 1 radiative flow from the exposed piece

% (T)OP (K fixed at end/nz+1, +K is exposed)
k = size(old,3); % habit and cleaner code.
parfor i=2:nx-1 % x
    for j=2:ny-1 % y

        change_in_temp(i,j,k) = ...
            ( ... % overall heat added for segment
            ...
            ...
            ... % X direction
            ( prod(PIECE_SIZE(2:3)) * ...
            ...
            ( spatial_heat_flow(old(i,j,k),old(i-1,j,k),PIECE_SIZE(1)) ...
            + spatial_heat_flow(old(i,j,k),old(i+1,j,k),PIECE_SIZE(1)))...
            )...
            ...
            + ... % Y dir
            ( prod(PIECE_SIZE([1 3])) * ...
            ...
            ( spatial_heat_flow(old(i,j,k),old(i,j-1,k),PIECE_SIZE(2)) ...
            + spatial_heat_flow(old(i,j,k),old(i,j+1,k),PIECE_SIZE(2)))...
            )...
            ...
            + ... % Z dir (THIS ONE HAS RADIATION)
            ( prod(PIECE_SIZE(1:2)) * ...
            ...
            ( spatial_heat_flow(old(i,j,k),old(i,j,k-1),PIECE_SIZE(3)) ...
            + radiate(old(i,j,k)))...
            )...
            );  
    end
end
        
change_in_temp(2:end-1,2:end-1,k) = ...
    change_in_temp(2:end-1,2:end-1,k) ./ ( SPECHEAT * PIECE_MASS );


% (F)RONT (I fixed at end, +I exposed)
i = size(old,2); % habit and cleaner code.
parfor j=2:ny-1 % y
    for k=3:(size(old,3)-1) % z
        change_in_temp(i,j,k) = ...
            ( ... % overall heat added for segment
            ...
            ...
            ... % X direction
            ( prod(PIECE_SIZE(2:3)) * ...
            ...
            ( spatial_heat_flow(old(i,j,k),old(i-1,j,k),PIECE_SIZE(1)) ...
            + radiate(old(i,j,k)))...
            )...
            ...
            + ... % Y dir
            ( prod(PIECE_SIZE([1 3])) * ...
            ...
            ( spatial_heat_flow(old(i,j,k),old(i,j-1,k),PIECE_SIZE(2)) ...
            + spatial_heat_flow(old(i,j,k),old(i,j+1,k),PIECE_SIZE(2)))...
            )...
            ...
            + ... % Z dir (THIS ONE HAS RADIATION)
            ( prod(PIECE_SIZE(1:2)) * ...
            ...
            ( spatial_heat_flow(old(i,j,k),old(i,j,k-1),PIECE_SIZE(3)) ...
            + spatial_heat_flow(old(i,j,k),old(i,j,k+1),PIECE_SIZE(3)))...
            )...
            );  
    end
end

change_in_temp(i,2:end-1,3:end-1) = ...
    change_in_temp(i,2:end-1,3:end-1) ./ ( SPECHEAT * PIECE_MASS );


%B(A)CK (I fixed at 1, -I exposed)

i = 1;
parfor j=2:ny-1 % y
    for k=3:(size(old,3)-1) % z
        change_in_temp(i,j,k) = ...
            ( ... % overall heat added for segment
            ...
            ...
            ... % X direction
            ( prod(PIECE_SIZE(2:3)) * ...
            ...
            ( spatial_heat_flow(old(i,j,k),old(i+1,j,k),PIECE_SIZE(1)) ...
            + radiate(old(i,j,k)))...
            )...
            ...
            + ... % Y dir
            ( prod(PIECE_SIZE([1 3])) * ...
            ...
            ( spatial_heat_flow(old(i,j,k),old(i,j-1,k),PIECE_SIZE(2)) ...
            + spatial_heat_flow(old(i,j,k),old(i,j+1,k),PIECE_SIZE(2)))...
            )...
            ...
            + ... % Z dir (THIS ONE HAS RADIATION)
            ( prod(PIECE_SIZE(1:2)) * ...
            ...
            ( spatial_heat_flow(old(i,j,k),old(i,j,k-1),PIECE_SIZE(3)) ...
            + spatial_heat_flow(old(i,j,k),old(i,j,k+1),PIECE_SIZE(3)))...
            )...
            );  
    end
end

change_in_temp(i,2:end-1,3:end-1) = ...
    change_in_temp(i,2:end-1,3:end-1) ./ ( SPECHEAT * PIECE_MASS );
%(L)EFT (J fixed at 1, -J exposed)
j = 1;
parfor i=2:nx-1 % y
    for k=3:(size(old,3)-1) % z
        change_in_temp(i,j,k) = ...
            ( ... % overall heat added for segment
            ...
            ...
            ... % X direction
            ( prod(PIECE_SIZE(2:3)) * ...
            ...
            ( spatial_heat_flow(old(i,j,k),old(i-1,j,k),PIECE_SIZE(1)) ...
            + spatial_heat_flow(old(i,j,k),old(i+1,j,k),PIECE_SIZE(1)))...
            )...
            ...
            + ... % Y dir
            ( prod(PIECE_SIZE([1 3])) * ...
            ...
            ( spatial_heat_flow(old(i,j,k),old(i,j+1,k),PIECE_SIZE(2)) ...
            + radiate(old(i,j,k)))...
            )...
            ...
            + ... % Z dir (THIS ONE HAS RADIATION)
            ( prod(PIECE_SIZE(1:2)) * ...
            ...
            ( spatial_heat_flow(old(i,j,k),old(i,j,k-1),PIECE_SIZE(3)) ...
            + spatial_heat_flow(old(i,j,k),old(i,j,k+1),PIECE_SIZE(3)))...
            )...
            );  
    end
end

change_in_temp(2:end-1,j,3:end-1) = ...
    change_in_temp(2:end-1,j,3:end-1) ./ ( SPECHEAT * PIECE_MASS );
%     (R)IGHT (J fixed at end, +J exposed)



% EDGES
%     TR
%     TL
%     TA
%     TF
%     BR %SPECIAL
%     BL %SPECIAL
%     BA %SPECIAL
%     BF %SPECIAL
% #DIV BY HEAT CAP AND MASS
% CORNERS
%     ALT
%     ART
%     FLT
%     FRT
%     ALB %SPECIAL
%     ARB %SPECIAL
%     FLB %SPECIAL
%     FRB %SPECIAL

% #DIV BY HEAT CAP AND MASS


% # MULTIPLY ALL BY TIMESTEP
    
if (~warned_yet) && abs(max(change_in_temp)) > Tgrill
        warning('VALUES ARE BLOWING UP');
        warned_yet = true;
end

temp_array = temp_array + change_in_temp;

    
end % time

%# TRIM GRILL TEMP OFF



end % main hamburger function


function [qpera] = radiate(Tint)
%
% ARGUMENTS
% 
%     Tint     The internal temperature of the radiating body,
%                  in degrees Celcius.
%                 
% RETURNS
% 
%     qpera    The rate of heat flow per square meter (J/m^2)


STFB = 5.60E-8; % stefan-boltzmann
CTOKEL = 273.15; %C to kelvin by adding. Kelvin to C by subtracting.
TEXT = 20 + CTOKEL;
Tint = Tint + CTOKEL;

qpera = (TEXT^4 - Tint^4) * STFB;

end % radiate function


function [qpera] = spatial_heat_flow(Tcurr, Tcomp, distance)
% 
% ARGUMENTS
% 
%     Tcurr     The temperature of the "current" piece in either degrees
%                     Celcius or in Kelvin.
%     Tcomp     The temperature of the piece to compare against the
%                   "current" piece. In degrees Celcius or in Kelvin.
%                   
% (Since we are using the difference between these two temperatures, either
% degrees C or Kelvin works.)
% 
% RETURNS
% 
%     qpera     The rate of heat flow per square meter (J/m^2)

global DENS;
global SPECHEAT;
global DIFFUS;

q = - DIFFUS * DENS * SPECHEAT * ( (Tcurr-Tcomp) / distance ) ;

end % spatial heat flow function