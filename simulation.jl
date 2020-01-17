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

	Θ = (280.4606 + 360.9856473 * t)/180*pi # gives radians

	mu = 3.986E5 #km3/s2
	R = 6371 #km
	T = 2*pi*sqrt.((c.alt.+R).^3/mu) # seconds
	n = 2*pi ./T #radians/sec

	r_pqr = fill(fill(0.0,3),c.num_sats)
	v_pqr = fill(fill(0.0,3),c.num_sats)
	for i in 1:c.num_sats
		plane_ind = Int(ceil(i/c.sats_per_plane))
		E = true_anomaly_to_ecc_anomaly(c.ν[i],c.ecc[plane_ind])
		r_pqr[i] = [(c.alt[plane_ind]+R)*(cosd(E)-c.ecc[plane_ind])
			(c.alt[plane_ind]+R)*sqrt(1-c.ecc[plane_ind]^2)*sind(E)
			0]
	
		Ė = n[plane_ind]/(1-c.ecc[plane_ind]*cosd(E)) # rad/s
		v_pqr[i] = [-(c.alt[plane_ind]+R)*Ė*sind(E)
			(c.alt[plane_ind]+R)*Ė*sqrt(1-c.ecc[plane_ind]^2)*cosd(E)
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

function ijk_to_eci(c::constellation,r_ijk,v_ijk,t)
	# this function converts from pqr to ijk coordinates for all satellites

	# INPUTS:
	# c: the constellation struct
	# r_ijk: array of the the ijk coordinates (km)
	# v_ijk: array of the velocities in ijk (km/s)
	# t: solar days since since Jan 1, 2000, 12 h

	# OUTPUTS:
	# r_eci: array of the eci coordinates (km)
	# v_eci: array of the eci velocities (km)

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

function keplerian_to_ECI(c::constellation, t)
	# this function takes a constellation struct and a time and converts it to 
	# a series of ECI positions and velocities for each satellite

	# INPUTS:
	# constellation: struct defined earlier
	# t: solar time in days since January 1, 2000 at 12h

	# OUTPUTS:
	# positions: array of 3x1 arrays that lists the positions of each satellite
	# velocities: array of 3x1 arrays that lists the velocities for each satellite

	r_pqr,v_pqr = keplerian_to_pqr(c)
	r_ijk,v_ijk = pqr_to_ijk(c,r_pqr,v_pqr)
	r_eci,v_eci = ijk_to_eci(c,r_ijk,v_ijk,t)

	return r_eci,v_eci
end

## Now let's generate some dynamics

function OrbitPlotter(x,p,t)
#this function plots an orbit of a satellite to get the position over a time period


r = x[1:3] #km
v = x[4:6] #km/s

#force balance
#gravity
#2-body gravity
f_grav=GM/(norm(r)^2)*-r/(norm(r)) #km/s^2

# atmospheric drag (assuming that time is invariant)
# date is new years on 2019
altitude=norm(r)-R_E #meters

#find the latitude and longitude in a geometric frame
GMST = (280.4606 + 360.9856473*(t/24/60/60 .+ (MJD_0.+t./(24*60*60)) .- 51544.5))./180*pi #calculates Greenwich Mean Standard Time
ROT = [cos(GMST) sin(GMST) 0;-sin(GMST) cos(GMST) 0; 0 0 1]
pos_ecef = ROT*r
lat= pos_ecef[3]/norm(pos_ecef)
long= atan(pos_ecef[2]/pos_ecef[1])

# out=SatelliteToolbox.nrlmsise00(DatetoJD(2019, 1, 1, 00, 00, 00), # Julian Day
    #              altitude,# Altitude [m]
    #              -lat*pi/180,# Latitude [rad]
    #              -long*pi/180,# Longitude [rad]
    #              83.7,# 81 day average F10.7 flux
    #              102.5,# Daily F10.7 for previous day
    #              15;# Magnetic index (daily)
    #              output_si = false# Output in cm^-3 and g/cm^-3
    # );
# ω_earth = [0;0;7.2921150E-5] #rad/s
# v_ecef = v-cross(ω_earth,r)
# f_drag=.5*out.den_Total*norm(v_ecef)^2*mass/BC*-v_ecef/norm(v_ecef)

#J2 oblateness term
J2=0.0010826359;
#J2=0;
f_J2=[J2*r[1]/norm(r)^7*(6*r[3]-1.5*(r[1]^2+r[2]^2))
     J2*r[2]/norm(r)^7*(6*r[3]-1.5*(r[1]^2+r[2]^2))
     J2*r[3]/norm(r)^7*(3*r[3]-4.5*(r[1]^2+r[2]^2))];


#acceleration
a=(f_grav+f_J2); #km/s^2

return [v;a]

end

