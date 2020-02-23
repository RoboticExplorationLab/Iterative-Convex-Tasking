# This code was generated to test the feasability of using satellite relative ranging
# measurements for better orbit determination.

# First, we'll establish the parameters of a satellite constellation

# TODO:
# - normalize all state vectors
# - find better way to transfer global variables

using LinearAlgebra
using SatelliteToolbox
using Plots
using UnicodePlots
using ForwardDiff
using SparseArrays

include("simulation.jl")
include("dynamics.jl")
include("satellite.jl")
include("sensors.jl")
include("least_squares.jl")

R_E = 6378.137

# true constellation (just one sat cluster for now)
num_sats = 1
num_planes = 1
sats_per_plane = Int(num_sats/num_planes)
alt = [800] # km
ecc = [.01] # unitless
inc = [57.5] # degrees
Ω = [10] # degrees
ω = [20] # degrees
M = [0]
# M = evenly_distribute(num_sats,num_planes,sats_per_plane) # degrees

c = generate_constellation(num_sats,num_planes,sats_per_plane,alt,ecc,inc,Ω,ω,M)

MJD = 7321+51544.5 #days since Jan 1,2000, 12 h converted to MJD

r_eci,v_eci = keplerian_to_ECI(c,MJD)

# This converts desired Keplerian parameters to an initial ECI
# estimated constellation (just assume one for now)
num_sats = 1
num_planes = 1
sats_per_plane = Int(num_sats/num_planes)
alt = [810] # km
ecc = [.01] # unitless
inc = [57.5] # degrees
Ω = [10] # degrees
ω = [20] # degrees
M = [0,10,20]
# M = evenly_distribute(num_sats,num_planes,sats_per_plane) # degrees

c_est = generate_constellation(num_sats,num_planes,sats_per_plane,alt,ecc,inc,Ω,ω,M)

r_eci_est,v_eci_est = keplerian_to_ECI(c_est,MJD)

# define a satellite and its ranging sensor
sat = define_satellite(10,[.15 .015 .015;.015 .15 .015;.015 .015 .15],2.2*10/.1)
sensor = sensor_struct(0,0)
perfect_sensor = sensor_struct(0,0)

## For satellites:

dt = 60

t = MJD:sec_to_days(dt):MJD+sec_to_days(600)

# X = fill(fill(fill(0.0,17),length(t)),num_sats)

# # X_simple = fill(fill(fill(0.0,17),length(t)),num_sats)

# @time for i in 1:num_sats
# 	# initial attitude
# 	ω_0 = [.1,.1,.1] # rad/s
# 	q_0 = [1,0,0,0] # unitless
# 	ρ_0 = [0,0,0] # kgm^2/s

# 	# add random dV to the satellite on the order of 1 m/s

# 	x_0 = [r_eci[i]/R_E;v_eci[i]+rand(3)/1000;ω_0;q_0;ρ_0;MJD]
# 	x_0_est = [r_eci_est[i]/R_E;v_eci_est[i];ω_0;q_0;ρ_0;MJD]

# 	# now simulate

# 	U = fill([ones(3);zeros(3)],length(t))

# 	X[i] = simulate_trajectory(x_0,U,t)

# 	# now do the same thing with J2 off (to generate error)

# 	# X_simple[i] = simulate_trajectory_simple(x_0_est,U,t)

# end

## get relative ranging measurements TODO: ORGANIZE THIS BETTER

# X_rel = zeros(num_sats,length(t),num_sats)
# for i in 1:num_sats
# 	for j in 1:length(t)
# 		for k in 1:num_sats
# 			X_rel[i,j,k] = rel_range(X_simple[i][j][1:3],X_simple[k][j][1:3],perfect_sensor)
# 		end
# 	end
# end

# g_estimate = zeros(num_sats*length(t))
# for j in 1:length(t)
# 	g_estimate[3*(j-1)+1:3*(j-1)+3] =
# 		[rel_range(X_simple[1][j][1:3],X_simple[2][j][1:3],perfect_sensor);
# 		rel_range(X_simple[2][j][1:3],X_simple[3][j][1:3],perfect_sensor);
# 		rel_range(X_simple[1][j][1:3],X_simple[3][j][1:3],perfect_sensor)]
# end

# Now state goes:
# [x1,1 x2,1 x3,1 x1,2 x2,2 ...]
# where x is just the position and velocity

# Build the gradient and Jacobian:
# f(x) passes in 18 state vectors and gets out the next 18
# g(x) passes in 18 state vectors and gets out 3 relative ranges

function f(x)
	# Needs state in form:
	# [x1;x2;x3] where x1 = [pos;vel]

	dt = 6

	return [rk4(x[1:6],0,posvel_dynamics_dot,dt)[1:6];
		rk4(x[7:12],0,posvel_dynamics_dot,dt)[1:6];
		rk4(x[13:18],0,posvel_dynamics_dot,dt)[1:6]]

end

function g(x)
	# Needs state in form:
	# [x1;x2;x3] where x1 = [pos;vel]

	return [rel_range(x[1:3],x[7:9],perfect_sensor);
		rel_range(x[7:9],x[13:15],perfect_sensor);
		rel_range(x[1:3],x[13:15],perfect_sensor)]

end

function unpackx(X)
	# this function just unpacks X into a vector
	N = length(X)
	m = size(X[1],1)
	x_unpacked = zeros(N*m)

	for k = 1:N
		x_unpacked[m*(k-1)+1:m*k] = X[k]
	end

	return x_unpacked
end

function packx(x_unpacked,m)
	# this function repacks X into the right form
	N = Int(length(x_unpacked)/m)

	X = fill(fill(0.0,m),N)

	for k = 1:N
		X[k] = x_unpacked[m*(k-1)+1:m*k]
	end

	return X

end

X_correct = fill(fill(0.0,18),length(t))
X_est = fill(fill(0.0,18),length(t))

# preturb by initial deltaV of 1 m/s
X_correct[1] = [r_eci[1]/R_E;v_eci[1]+rand(3)/1000;r_eci[1]/R_E;v_eci[1]+rand(3)/1000;r_eci[1]/R_E;v_eci[1]+rand(3)/1000]

for j = 2:length(t)
	X_correct[j] = f(X_correct[j-1])
end

for j in 1:length(t)
	X_est[j] = X_correct[j] + rand(18)/100
end

X_est-X_correct

Y = fill(fill(0.0,3),length(t))
for j in 1:length(t)
	Y[j] =
		[rel_range(X_correct[j][1:3],X_correct[j][7:9],perfect_sensor);
		rel_range(X_correct[j][7:9],X_correct[j][13:15],perfect_sensor);
		rel_range(X_correct[j][1:3],X_correct[j][13:15],perfect_sensor)]
end

R = 1.0*Matrix(I,3,3) # weighting matrices
Q = 1.0*Matrix(I,18,18)

function traj_cost(X)
	# given an input trajectory and a set of measurements
	# determine the cost that you want to minimize

	# INPUTS:
	# X: toal vector of trajectory (n states)
	# GLOBALS:
	# R: cost of measurement residual
	# Q: cost of dynamics residual
	# Y: total vector of measurements
	# f: dynamics function
	# g: cost function
	# N: number of time steps per input
	# m: number of states per time step

	# OUTPUTS:
	# T: total cost

	m = 3*6
	N = length(t)
	T = 0
	for k = 1:N-1
		T += 1/2*(Y[k]-g(X[m*(k-1)+1:m*k]))'*inv(R)*(Y[k]-g(X[m*(k-1)+1:m*k])) +
			(X[m*(k)+1:m*(k+1)]-f(X[m*(k-1)+1:m*k]))'*inv(Q)*(X[m*(k)+1:m*(k+1)]-f(X[m*(k-1)+1:m*k]))
	end
	T += 1/2*(Y[N]-g(X[m*(N-1)+1:m*N]))'*inv(R)*(Y[N]-g(X[m*(N-1)+1:m*N]))

	return T

end

# cfg_grad = ForwardDiff.GradientConfig(traj_cost, unpackx(X_est))
# cfg_jac = ForwardDiff.HessianConfig(traj_cost, unpackx(X_est))

grad = x -> ForwardDiff.gradient(traj_cost,x)
jac = x -> ForwardDiff.hessian(traj_cost,x)
iters = 3
err = zeros(iters+1)
err[1] = traj_cost(unpackx(X_est))
# λ = zeros(iters+1)
# λ[1] = 1E8 # set to 0 for GN
#
# Newton step parameters
# ρ_αu = 1.2
# ρ_αd = 0.5
# ρ_ls = 0.01
# ρ_λd = 0.5
# ρ_λu = 1.0
#
# @time let X_est = X_est
# 	α = 1
# 	for i in 1:iters
# 		λ_temp = λ[i]
#
# 		∇f = grad(unpackx(X_est))
#
# 		∇f2 = jac(unpackx(X_est))
#
# 		δx = -inv(∇f2'*∇f2+λ_temp*Matrix(I,size(∇f2,1),size(∇f2,1)))*∇f2'*∇f
#
# 		# line search
# 		while traj_cost(unpackx(X_est) + α*δx) > traj_cost(unpackx(X_est)) + ρ_ls*∇f'*(α*δx)
# 			α *= ρ_αd
#
# 			#reassign λ and recompute δx
# 			λ_temp *= ρ_λu
# 			δx = -inv(∇f2'*∇f2+λ_temp*Matrix(I,size(∇f2,1),size(∇f2,1)))*∇f2'*∇f
# 		end
#
# 		X_est += packx(α*δx,18)
#
# 		# reassign step size
# 		α = min(ρ_αu*α,1)
#
# 		# reconfigure λ
# 		λ[i+1] = λ_temp*ρ_λd
#
# 		err[i+1] = traj_cost(unpackx(X_est))
#
# 		println("Iteration: ",i)
# 	end
# end

#Standard Newton Method
@time let X_est = X_est
	α = 1
	for i in 1:iters

		current_cost = traj_cost(unpackx(X_est))

		∇f = grad(unpackx(X_est))
		∇f2 = jac(unpackx(X_est))

		δx = -∇f2\∇f

		# dumbest possible line search
		α = 1.0
		new_X = unpackx(X_est) + α*δx
		new_cost = traj_cost(new_X)
		while  new_cost > current_cost
			α = 0.5*α
			new_X = unpackx(X_est) + α*δx
			new_cost = traj_cost(new_X)
		end
		println("alpha= ",α)
		X_est += packx(new_X,18)

		err[i+1] = new_cost

		println("Iteration: ",i)
	end
end

plot([err,λ], layout = 2)
