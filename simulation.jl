## This file contains the basis functions for simulating the satellite constellation

# Marcos:
macro name(arg)
    x = string(arg)
    quote
        $x
    end
end

## Let's define a custom structure for the constelaltion with all the important data
struct constellation
	num_sats::Int64
	num_planes::Int64
	sats_per_plane::Int64
	alt::AbstractArray{Float64,1}
	ecc::AbstractArray{Float64,1}
	inc::AbstractArray{Float64,1}
	Ω::AbstractArray{Float64,1}
	ω::AbstractArray{Float64,1}
	ν::AbstractArray{Float64,1}
end

function generate_constellation(num_sats,num_planes,sats_per_plane,alt,ecc,inc,Ω,ω,ν)
	# this function takes generates a constellation struct
	# 
	# INPUTS:
	# num_sats: Int 
	# num_planes: Int 
	# alt (altitudes): array of size 1 x num_planes (km)
	# ecc (eccentricies): array of size 1 x num_planes (unitless)
	# inc (inclinations): array of size 1 x num_planes (degrees)
	# Ω (RAAN): array of size 1 x num_planes (degrees)
	# ω (arg of periapsis): array of size 1 x num_planes (degrees)
	# ν (true anomaly): array of size 1 x num_planes (degrees)
	# 
	# OUTPUTS:
	# constellation (if possible)

	size_flag = true
	for x = [alt,ecc,inc,Ω,ω]
		if size(x) != (num_planes,)
			print("The size of the input, ", x, " is wrong.\n")
			size_flag = false
		end
	end
	if size(ν) != (num_sats,)
		print("the size of the ν vector is wrong.\n")
	end
	if size_flag
		return constellation(num_sats,num_planes,sats_per_plane,alt,ecc,inc,Ω,ω,ν)
	end
end

function evenly_distribute(num_sats,num_planes,sats_per_plane)
	# this function generates the ν values for an evently distributed constellation

	# INPUTS:
	# num_sats: Int
	# num_planes: Int
	# sats_per_plane: Int

	# OUTPUS:
	# ν: array of size 1xnum_sats that has all the true anomalies for the satellites

	ν = zeros(num_sats)
	for i in 1:num_planes
		ν[(i-1)*sats_per_plane+1:i*sats_per_plane] = collect(0:360/sats_per_plane:360-360/sats_per_plane)
	end

	return ν
end

## Now list conversions between orbital parameters to ECI for dynamics stuff

function R_x(a)
	# in degrees
	return [1 0 0;0 cosd(a) sind(a);0 -sind(a) cosd(a)]
end

function R_z(a)
	# in degrees
	return [cosd(a) sind(a) 0;-sind(a) cosd(a) 0;0 0 1]
end

function true_anomaly_to_ecc_anomaly(ν,e)
	# this function takes in a true anomaly and outputs an eccentric anomaly

	# INPUTS:
	# ν: true anomaly (degrees)
	# e: eccentricity (unitless)

	# OUTPUTS:
	# E: eccentric anomaly (degrees)

	E = 2*atan(sqrt(1-e)*tand(ν/2),sqrt(1-e))*180/pi
	return E
end

function keplerian_to_pqr(c::constellation)
	# this function converts a constellation and a time into a series of r_pqr

	# INPUTS:
	# c: constellation struct 
	# t: time since Jan 1, 2000, 12 h in solar days 

	# OUTPUTS:
	# r_pqr: array of r_pqr positions for each satellite
	# v_pqr: array of velocities for each satellite in pqr

	mu = 3.986E5 #km3/s2
	R_E = 6371 #km
	T = 2*pi*sqrt.((c.alt.+R_E).^3/mu) # seconds
	n = 2*pi ./T #radians/sec

	r_pqr = fill(fill(0.0,3),c.num_sats)
	v_pqr = fill(fill(0.0,3),c.num_sats)
	for i in 1:c.num_sats
		plane_ind = Int(ceil(i/c.sats_per_plane))
		E = true_anomaly_to_ecc_anomaly(c.ν[i],c.ecc[plane_ind])
		r_pqr[i] = [(c.alt[plane_ind]+R_E)*(cosd(E)-c.ecc[plane_ind])
			(c.alt[plane_ind]+R_E)*sqrt(1-c.ecc[plane_ind]^2)*sind(E)
			0]
	
		Ė = n[plane_ind]/(1-c.ecc[plane_ind]*cosd(E)) # rad/s
		v_pqr[i] = [-(c.alt[plane_ind]+R_E)*Ė*sind(E)
			(c.alt[plane_ind]+R_E)*Ė*sqrt(1-c.ecc[plane_ind]^2)*cosd(E)
			0]

	end

	return r_pqr,v_pqr
end

function pqr_to_ijk(c::constellation,r_pqr,v_pqr)
	# this function converts from pqr to ijk coordinates for all satellites

	# INPUTS:
	# c: the constellation struct
	# r_pqr: array of the the pqr coordinates (km)
	# v_pqr: array of the velocities in pqr (km/s)

	# OUTPUTS:
	# r_ijk: array of the ijk coordinates (km)
	# v_ijk: array of the ijk velocities (km)

	r_ijk = fill(fill(0.0,3),c.num_sats)
	v_ijk = fill(fill(0.0,3),c.num_sats)

	conv = []
	for i in 1:num_planes
		push!(conv,R_z(-c.Ω[i])*R_x(-c.inc[i])*R_z(-c.ω[i]))
	end

	for i in 1:c.num_sats
		plane_ind = Int(ceil(i/c.sats_per_plane))
		r_ijk[i] = conv[plane_ind]*r_pqr[i]
		v_ijk[i] = conv[plane_ind]*v_pqr[i]
	end

	return r_ijk,v_ijk
end

function ijk_to_eci(c::constellation,r_ijk,v_ijk,MJD)
	# this function converts from pqr to ijk coordinates for all satellites

	# INPUTS:
	# c: the constellation struct
	# r_ijk: array of the the ijk coordinates (km)
	# v_ijk: array of the velocities in ijk (km/s)
	# t: time in MJD

	# OUTPUTS:
	# r_eci: array of the eci coordinates (km)
	# v_eci: array of the eci velocities (km)

	t = MJD - 51544.5
	conv = R_z(280.4606 + 360.9856473 * t)

	r_eci = fill(fill(0.0,3),c.num_sats)
	v_eci = fill(fill(0.0,3),c.num_sats)

	for i in 1:c.num_sats
		plane_ind = Int(ceil(i/c.sats_per_plane))
		r_eci[i] = conv*r_ijk[i]
		v_eci[i] = conv*v_ijk[i]
	end

	return r_eci,v_eci
end

function keplerian_to_ECI(c::constellation, MJD)
	# this function takes a constellation struct and a time and converts it to 
	# a series of ECI positions and velocities for each satellite

	# INPUTS:
	# constellation: struct defined earlier
	# t: time in MJD

	# OUTPUTS:
	# positions: array of 3x1 arrays that lists the positions of each satellite
	# velocities: array of 3x1 arrays that lists the velocities for each satellite

	r_pqr,v_pqr = keplerian_to_pqr(c)
	r_ijk,v_ijk = pqr_to_ijk(c,r_pqr,v_pqr)
	r_eci,v_eci = ijk_to_eci(c,r_ijk,v_ijk,MJD)

	return r_eci,v_eci
end



function get_mag_field(X,t)

	# INPUTS:
	# X: lists of states
	# [r_eci/R_E;
	# v_eci;
	# ω; r
	# q;
	# ρ;
	# MJD]
	# U: control over entire trajectory
	# t: time vector in MJD

	# OUTPUTS:
	# B_ECI: 

	B_ECI = zeros(length(t),3)

	for i in 1:length(t)

		r = X[i][1:3]*R_E #km
		v = X[i][4:6] #km/s

		MJD = X[i][17] - 51544.5
		t = MJD - 51544.5

		GMST = (280.4606 + 360.9856473*(t))./180*pi #calculates Greenwich Mean Standard Time
		ROT = [cos(GMST) sin(GMST) 0;-sin(GMST) cos(GMST) 0; 0 0 1]
		pos_ecef = ROT*r
		lat= pos_ecef[3]/norm(pos_ecef)
		long= atan(pos_ecef[2]/pos_ecef[1])
		NED_to_ENU = [0 1 0;1 0 0;0 0 -1]
		R_ENU_to_XYZ =
	            [-sin(long) -sin(lat)*cos(long) cos(lat)*cos(long);
	            cos(long) -sin(lat)*sin(long) cos(lat)*sin(long);
	            0 cos(lat) sin(lat)]
		B_ECI[i,:] = Rz(GMST)'* R_ENU_to_XYZ *NED_to_ENU * collect(igrf12syn(0, 2019, 2,(norm(r)), lat*180/pi, ((long+2π) % 2π)*180/pi)[1:3])*1E-9
	end

	return B_ECI

end


function simulate_trajectory(x_0,U,t)
	# this function simulates a satellite trajectory over a given time vector

	# INPUTS:
	# x_0: initial state
	# [r_eci/R_E;
	# v_eci;
	# ω; r
	# q;
	# ρ;
	# MJD]
	# U: control over entire trajectory
	# t: time vector in MJD

	# OUTPUTS:
	# X: state over all time points

	X = fill(fill(0.0,17),length(t))

	X[1] = x_0

	B_ECI_log = fill(fill(0.0,3),length(t))

	for i in 1:length(t)-1

		X[i+1]= rk4(X[i],U[i],dynamics_dot,days_to_sec(t[i+1]-t[i]))

	end

	return X

end

function days_to_sec(days)
	# simple function that converts between days and seconds

	return days*24*60*60

end

function sec_to_days(seconds)
	# simple function that converts between seconds and days

	return seconds/60/60/24

end

