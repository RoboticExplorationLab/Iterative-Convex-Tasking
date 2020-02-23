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
using ForwardDiff

include("simulation.jl")
include("dynamics.jl")
include("satellite.jl")
include("sensors.jl")
include("least_squares.jl")

R_E = 6378.137

num_sats = 3
num_planes = 1
sats_per_plane = Int(num_sats/num_planes)
alt = [800] # km
ecc = [.1] # unitless
inc = [57.5] # degrees
Ω = [10] # degrees
ω = [20] # degrees
M = [0,10,20]
# M = evenly_distribute(num_sats,num_planes,sats_per_plane) # degrees

c = generate_constellation(num_sats,num_planes,sats_per_plane,alt,ecc,inc,Ω,ω,M)

MJD = 7321+51544.5 #days since Jan 1,2000, 12 h converted to MJD

r_eci,v_eci = keplerian_to_ECI(c,MJD)

# define a satellite and its ranging sensor
sat = define_satellite(10,[.15 .015 .015;.015 .15 .015;.015 .015 .15])
sensor = sensor_struct(0,0)

## For satellites:

t = MJD:sec_to_days(.1):MJD+sec_to_days(60*60*1.5)

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

## get relative ranging measurements

X_rel = zeros(num_sats,length(t),num_sats)
for i in 1:num_sats
	for j in 1:length(t)
		for k in 1:num_sats
			X_rel[i,j,k] = rel_range(X[i][j][1:3],X[k][j][1:3],sensor)
		end
	end
end

three_sat_X_rel = zeros(num_sats*length(t))
for j in 1:length(t)
	three_sat_X_rel[3*(j-1)+1:3*(j-1)+3] = [X_rel[1,j,2];
		X_rel[2,j,3];
		X_rel[1,j,3]]
end

## levenberg marquardt section
# initial guess
a = 6378.137+800
e = .1
inc = 57.5
Ω = 10
ω = 20
M_1 = 0
M_2 = 5
M_3 = 20

iters = 1

β_vec = fill(fill(0.0,8),iters)
β_vec[1] = [a,e,inc,Ω,ω,M_1,M_2,M_3]
λ = 1E-4
β = β_vec[1]

@time for k in 1:iters-1
	n = 8 # number of scaling parameters
	m = 3 # number of outputs from function
	meas = length(t) # number of measurements

	# get jacobian
	J = zeros(m*meas,n)
	g = x -> ForwardDiff.jacobian(three_sat_rel_ranges,x);
	for k in 1:length(t)
		x = [β;(t[k]-t[1])*24*60*60]
		J[m*(k-1)+1:m*(k-1)+m,:] = g(x)[:,1:n]
	end

	# find all measurements
	f_meas = zeros(m*meas)
	for k in 1:length(t)
		x = [β;(t[k]-t[1])*24*60*60]
		f_meas[m*(k-1)+1:m*(k-1)+m] = three_sat_rel_ranges(x)
	end

	δ = inv(J'*J+λ*diagm(diag(J'*J)))*J'*(three_sat_X_rel-f_meas)
	
	f_meas_new = zeros(m*meas)
	for k in 1:length(t)
		x = [β+δ;(t[k]-t[1])*24*60*60]
		f_meas_new[m*(k-1)+1:m*(k-1)+m] = three_sat_rel_ranges(x)
	end
 
 	if ((three_sat_X_rel-f_meas_new)'*(three_sat_X_rel-f_meas_new)
 		 < (three_sat_X_rel-f_meas)'*(three_sat_X_rel-f_meas))
		β += δ
		λ/=1.5
	else
		λ*=5
	end

end

