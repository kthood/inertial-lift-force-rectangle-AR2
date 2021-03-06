README - INERTIAL MIGRATION RECTANGULAR CHANNEL AR=2 CODE
Kaitlyn Hood 9/1/2015
DOI 10.5281/zenodo.832327

For more details about the derivation of this code, see the preprint "Direct
Measurement of Particle Inertial Migration in Rectangular Channels" by K. 
Hood, S. Kahkeshani, D. Di Carlo, and M. Roper at: 
http://arxiv.org/abs/1509.01643

Cite as:
Kaitlyn Hood. (2017, July 19). kthood/inertial-lift-force-rectangle-AR2: First release of inertial lift force rectangle AR2. Zenodo. http://doi.org/10.5281/zenodo.832328

%——————————————————————————————————————————————————————————————————————————%
getLiftForce_AR2.m

REQUIRES: HoLeal_channel_AR2_Re1_mesh8_08-25-15.mat

INPUT: 		x 		- x coordinate of center of particle
		y 		- y coordinate of center of particle
		rsph 		- radius of particle
		U		- maximum velocity of background flow
		rho		- fluid density
		L		- side length of square channel
		
OUTPUT:		forcex	- lift force in the x-direction including order rsph^5
		forcey 	- lift force in the y-direction including order rsph^5
		
NOTES:	Results may be inaccurate for large particle radius (rsph > .2L)
		

Code computes the fourth and fifth rider terms of the (inertial) lift force 
of a particle in a square channel.  
Side length of square channel is L, and the coordinates are chosen:
	-.5L < x < .5L, 	-.5L < y < .5L
The result is in dimensional coordinates.

%——————————————————————————————————————————————————————————————————————————%
getLiftForceO4_AR2.m

REQUIRES: HoLeal_channel_AR2_Re1_mesh8_08-25-15.mat

INPUT: 		x 		- x coordinate of center of particle
		y 		- y coordinate of center of particle
		rsph 		- radius of particle
		U		- maximum velocity of background flow
		rho		- fluid density
		L		- side length of square channel
		
OUTPUT:		forcex	- lift force in the x-direction including order rsph^4
		forcey 	- lift force in the y-direction including order rsph^4
		
NOTES:	Results may be inaccurate for moderate to large particle radius
		

Code computes the fourth order term of the (inertial) lift force of a particle 
in a square channel.  
Side length of square channel is L, and the coordinates are chosen:
	-.5L < x < .5L, 	-.5L < y < .5L
The result is in dimensional coordinates.
%——————————————————————————————————————————————————————————————————————————%
getLiftImage_AR2.m

REQUIRES: HoLeal_channel_AR2_Re1_mesh8_08-25-15.mat

INPUT: 		x 		- x coordinate of center of particle
		y 		- y coordinate of center of particle
		rsph 		- radius of particle
		L		- side length of square channel

OUTPUT
  v1 - image velocity due to the stokeslet in the x-direction (2x1 vector)
  v2 - image velocity due to the stokeslet in the x-direction (2x1 vector)

NOTES:	Results may be inaccurate for moderate to large particle radius
		

Code computes the image velocities of the (inertial) lift force of a particle 
in a square channel.  
Side length of square channel is L, and the coordinates are chosen:
	-.5L < x < .5L, 	-.5L < y < .5L
The result is in non-dimensional coordinates. 

%——————————————————————————————————————————————————————————————————————————%
liftforce_AR2_channel_plot.m

REQUIRES: 	HoLeal_channel_AR2_Re1_mesh8_08-25-15.mat
		getLiftForce_AR2.m

INPUT:

OUTPUT:		plot of the lift force field (.eps)

%——————————————————————————————————————————————————————————————————————————%
animation_focusing_AR2.m

REQUIRES: 	getLiftForce_AR2.m
		inertial_constants_AR2_n201_09012015.mat

INPUT:

OUTPUT:		animation of particles moving according to the lift forces 
		in a square channel (.gif)

%——————————————————————————————————————————————————————————————————————————%

