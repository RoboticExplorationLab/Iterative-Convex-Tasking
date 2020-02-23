## This file contains all the functions for least squares curve fittings

## first let's do a 3 satellite model 
## assume the input parameters are as follows:
	# a: semi-major (km)
	# e: eccentricity (unitless)
	# i: inclination (deg)
	# Ω: RAAN (deg)
	# ω: arg of perigee (deg)
		# M_1: first sat mean anomaly 
		# M_2: second sat mean anomaly
		# M_3: thrid sat mean anomaly (deg)

function levenberg_marquardt_update(y,f,β,t,λ)
	# this function develops the δ for a levenberg_marquardt update

	# INPUTS:
	# y: the vector of measurements over a time period 
	# f: function (f(x,β))
	# β: the scaling parameters (in this case)
		# a: semi-major (km)
		# e: eccentricity (unitless)
		# i: inclination (deg)
		# Ω: RAAN (deg)
		# ω: arg of perigee (deg)
		# M_1: first sat mean anomaly 
		# M_2: second sat mean anomaly
		# M_3: thrid sat mean anomaly (deg)
	# t: the time vector in MJD
	# λ: the scaling factor

	# OUTPUTS:
	# δ: the change in the scaling parameters

	# assume that y is [y1,y2,...]
	# find the size of the measurements

	n = 8 # number of scaling parameters
	m = 3 # number of outputs from function
	meas = length(t) # number of measurements

	# get jacobian
	J = zeros(m*meas,n)
	g = x -> ForwardDiff.jacobian(f,x);
	for k in 1:length(t)
		x = [β;(t[k]-t[1])*24*60*60]
		J[m*(k-1)+1:m*(k-1)+m,:] = g(x)[:,1:n]
	end

	# find all measurements
	f_meas = zeros(m*meas)
	for k in 1:length(t)
		x = [β;(t[k]-t[1])*24*60*60]
		f_meas[m*(k-1)+1:m*(k-1)+m] = f(x)
	end

	δ = inv(J'*J+λ*diagm(diag(J'*J)))*J'*(y-f_meas)

	return δ

end

function three_sat_rel_ranges(x::Vector)
	# this function generates the relative ranges between
	# the three satellites

	# INPUT:
	# β:
		# a: semi-major (km)
		# e: eccentricity (unitless)
		# i: inclination (deg)
		# Ω: RAAN (deg)
		# ω: arg of perigee (deg)
		# M_1: first sat mean anomaly 
		# M_2: second sat mean anomaly
		# M_3: thrid sat mean anomaly (deg)
		# t: time since start of orbit

	# OUTPUT:
	# r_rel: relative distnaces (km) between sats in the form
	# 	[r_1->2; r_2->3; r_1->3]

	μ = 3.986E5 #km3/s2 

	#unpack
	a = x[1]
	e = x[2]
	i = x[3]
	Ω = x[4]
	ω = x[5]
	t = x[9]

	n = sqrt(μ/a^3)
	M1 = x[6]+n*t
	M2 = x[7]+n*t
	M3 = x[8]+n*t

	rel_1 = norm(kep_to_eci(a,e,i,Ω,ω,M1,t) - kep_to_eci(a,e,i,Ω,ω,M2,t))
	rel_2 = norm(kep_to_eci(a,e,i,Ω,ω,M2,t) - kep_to_eci(a,e,i,Ω,ω,M3,t))
	rel_3 = norm(kep_to_eci(a,e,i,Ω,ω,M1,t) - kep_to_eci(a,e,i,Ω,ω,M3,t))

	return [rel_1,rel_2,rel_3]

end

function kep_to_eci(a,e,inc,Ω,ω,M,t)
	# this function simply converts to eci 
	# INPUTS:
	# keplecrian elements
	# t: time since start of orbit

	# OUTPUTS:
	# r_eci (km)

	mu = 3.986E5 #km3/s2
	R_E = 6378.137 #km
	T = 2*pi*sqrt.(a^3/mu) # seconds
	n = 2*pi ./T #radians/sec
	J2 = 1.0826E-3

	# Gauss variational equations 
	δΩ = -1.5*n*J2*(R_E/a/(1-e^2))^2*cos(inc)*180/pi*t
	δω = .75*n*J2*(R_E/a/(1-e^2))^2*(5*cos(inc)-1)*180/pi*t

	# convert to radians 
	M = M/180*pi

	# kepler's equation
	E = pi
	for i in 1:10
		E = e*sin(E)+M
	end

	r_eci = R_z(-Ω+δΩ)*R_x(-inc)*R_z(-ω+δω)*[a*(cos(E)-e);a*sqrt(1-e^2)*sin(E);0]

	# r_eci = [(-a*sin(E)*sqrt(1-e^2)*(cos(Ω)*sin(ω)+sin(Ω)*cos(i)*cos(ω)) 
	# 	- a*(e-cos(E))*(cos(Ω)*cos(ω)-sin(Ω)*cos(i)*sin(ω)));
	# 	(-a*sin(E)*sqrt(1-e^2)*(sin(Ω)*sin(ω)-cos(Ω)*cos(i)*cos(ω)) 
	# 	- a*(e-cos(E))*(sin(Ω)*cos(ω)+cos(Ω)*cos(i)*sin(ω)));
	# 	-a*sin(i)*sin(ω)*(e-cos(E))+a*sin(E)*cos(ω)*sin(i)*sqrt(1-e^2)]

	return r_eci

end

function R_x(a)
	# in degrees
	return [1 0 0;0 cosd(a) sind(a);0 -sind(a) cosd(a)]
end

function R_z(a)
	# in degrees
	return [cosd(a) sind(a) 0;-sind(a) cosd(a) 0;0 0 1]
end

function levenberg_marquardt(iters,X_est,Y,f,g,R,Q)
	J = zeros((3*length(t))+(18*length(t)-18),18*length(t))
	residuals = zeros(size(Y,1)+size(X_est,1)-18)

	RMSE = zeros(iters)
	RMSE[1] = sqrt(sum((X_est-X_correct).^2)/length(X_est))
	for i in 2:iters
		for k in 1:length(t)
			# get jacobians of functions 
			G_temp = G(X_est[18(k-1)+1:18*k])
			F_temp = F(X_est[18(k-1)+1:18*k])

			# fill out residuals 
			residuals[3*(k-1)+1:3*(k-1)+3] = (
				Y[3*(k-1)+1:3*(k-1)+3] - g(X_est[18(k-1)+1:18*k])
				)

			J[3(k-1)+1:3*k,18(k-1)+1:18*k] = (
				-sqrt(inv(R))*G_temp
				)

			if k < length(t)
				J[18(k-1)+1+3*length(t):18(k-1)+18+3*length(t),18(k-1)+1:18*k] = (
					-sqrt(inv(Q))*F_temp
					)
				J[18(k-1)+1+length(t):18(k-1)+18+length(t),18k+1:18(k+1)] = (
					sqrt(inv(Q))
					)
				residuals[18*(k-1)+1+3*length(t):18*(k-1)+18+3*length(t)] = (
					X_est[18(k)+1:18(k+1)] - f(X_est[18(k-1)+1:18k])
					)

			end

		end

		# Levenberg Marquardt
		n = 8 # number of scaling parameters
		m = 3 # number of outputs from function
		meas = length(t) # number of measurements
		λ = 1000

		# make sparse
		Js = sparse(J)


		δ = \(Js'*Js+λ*sparse(diagm(diag(Js'*Js))),Js'*(residuals))

		X_est -= δ

		RMSE[i] = sqrt(sum((X_est-X_correct).^2)/length(X_est))

		# if RMSE[i] > RMSE[i-1]
		# 	λ *= 5
		# else 
		# 	λ /= 1.2
		# end
	end
	return RMSE
end


