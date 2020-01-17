# This code was generated to test the feasability of using convex optimization 
# for the tasking of Earth Observation constellations. It uses a combination of
# ADMM, iLQR, and cost function heuristics. It is meant to be compared to 
# the current standards in cooperative tasking like MILP (from optimization
# research), as well as random assignment and genetic algorithms. 

# First, we'll establish the parameters of a satellite constellation

using LinearAlgebra

include("simulation.jl")

num_sats = 100
num_planes = 2
sats_per_plane = Int(num_sats/num_planes)
alt = [800,1200] # km
ecc = [0,0] # unitless
inc = [57.5,96.6] # degrees
Ω = [0,0] # degrees
ω = [0,0] # degrees
ν = evenly_distribute(num_sats,num_planes,sats_per_plane) # degrees

c = generate_constellation(num_sats,num_planes,sats_per_plane,alt,ecc,inc,Ω,ω,ν)

t = 2000 #days since Jan 1,2000, 12 h

r_eci,v_eci = keplerian_to_ECI(c,t)

# the state for each satellite is [r,v]
