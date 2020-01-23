# This code was generated to test the feasability of using convex optimization 
# for the tasking of Earth Observation constellations. It uses a combination of
# ADMM, iLQR, and cost function heuristics. It is meant to be compared to 
# the current standards in cooperative tasking like MILP (from optimization
# research), as well as random assignment and genetic algorithms. 

# First, we'll establish the parameters of a satellite constellation

# TODO: 
# - normalize all state vectors
# - find better way to transfer global variables

using LinearAlgebra
using SatelliteToolbox
using Plots
using UnicodePlots

include("simulation.jl")
include("dynamics.jl")
include("satellite.jl")

R_E = 6371

num_sats = 4
num_planes = 2
sats_per_plane = Int(num_sats/num_planes)
alt = [800,1200] # km
ecc = [0,0] # unitless
inc = [57.5,96.6] # degrees
Ω = [0,0] # degrees
ω = [0,0] # degrees
ν = evenly_distribute(num_sats,num_planes,sats_per_plane) # degrees

c = generate_constellation(num_sats,num_planes,sats_per_plane,alt,ecc,inc,Ω,ω,ν)

MJD = 7321+51544.5 #days since Jan 1,2000, 12 h converted to MJD

r_eci,v_eci = keplerian_to_ECI(c,MJD)

# define a satellite
sat = define_satellite(10,[.15 .015 .015;.015 .15 .015;.015 .015 .15])

## For satellites:

t = MJD:sec_to_days(.1):MJD+1/24/60

X = fill(fill(fill(0.0,17),length(t)),num_sats)

@time for i in 1:num_sats
	# initial attitude
	ω_0 = [.1,.1,.1] # rad/s
	q_0 = [1,0,0,0] # unitless
	ρ_0 = [0,0,0] # kgm^2/s

	x_0 = [r_eci[i]/R_E;v_eci[i];ω_0;q_0;ρ_0;MJD]

	# now simulate

	U = fill([ones(3);zeros(3)],length(t))

	X[i] = simulate_trajectory(x_0,U,t)
end

# find the magnetic field for the first satellite 

B_ECI_temp = get_mag_field(X[1],t)
Plots.plot(collect(t),B_ECI_temp)

