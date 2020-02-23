using LinearAlgebra
using SatelliteToolbox

## FIRST, VERY IMPORTANT:
# STATE:
# [r_eci/R_E (unitless);
# v_eci (km/s);
# ω; (rad/s)
# q; (unitless)
# ρ; (kgm^2/s^2)
# MJD]

## Now let's generate some dynamics

function dynamics_dot_simple(x,u)
	# INPUTS:
	# STATE:
	# [r_eci/R_E;
	# v_eci;
	# ω; r
	# q;
	# ρ;
	# MJD]

	# IMPORTANT:
	# assume you have J, mass defined

	# u: 
	# [m (Am^2);
	# ρ̇ (kgm^2/s^2)]

	# OUTPUTS:
	# [v;a]

	R_E = 6378.137 #km

	r = x[1:3]*R_E #km
	v = x[4:6] #km/s
	ω = x[7:9] # rad/s
	q = x[10:13] # unitless 
	ρ = x[14:16]
	MJD = x[17]

	m = u[1:3] 
	ρ̇ = u[4:6]

	#force balance
	#gravity
	#2-body gravity
	mu = 3.986E5 #km3/s2
	f_grav=mu/(norm(r)^2)*-r/(norm(r)) #km/s^2

	# atmospheric drag (assuming that time is invariant)
	# date is new years on 2019

	#find the latitude and longitude in a geometric frame
	t = MJD - 51544.5
	GMST = (280.4606 + 360.9856473*(t))./180*pi #calculates Greenwich Mean Standard Time
	ROT = [cos(GMST) sin(GMST) 0;-sin(GMST) cos(GMST) 0; 0 0 1]
	pos_ecef = ROT*r
	lat= asin(pos_ecef[3]/norm(pos_ecef))
	long= atan(pos_ecef[2],pos_ecef[1])

	#acceleration
	a=(f_grav) #km/s^2

	# magnetic field
	NED_to_ENU = [0 1 0;1 0 0;0 0 -1]
	R_ENU_to_XYZ =
            [-sin(long) -sin(lat)*cos(long) cos(lat)*cos(long);
            cos(long) -sin(lat)*sin(long) cos(lat)*sin(long);
            0 cos(lat) sin(lat)]
	B_ECI = Rz(GMST)'* R_ENU_to_XYZ * NED_to_ENU * igrf12(MJD_to_year(MJD),norm(r)*1E3, lat, long)*1E-9
	
	#rotation (Euler)
	ω̇ = inv(sat.J)*(cross(m,B_ECI) - ρ̇ - cross(ω,sat.J*ω+ρ))

	#quaternion 
	q̇ = 1/2*(Quaternion(q)*Quaternion([0;ω]))[:]

	return [v/R_E;a;ω̇;q̇;ρ̇;0]

end


function dynamics_dot(x,u)
	# INPUTS:
	# STATE:
	# [r_eci/R_E;
	# v_eci;
	# ω; r
	# q;
	# ρ;
	# MJD]

	# IMPORTANT:
	# assume you have J, mass defined

	# u: 
	# [m (Am^2);
	# ρ̇ (kgm^2/s^2)]

	# OUTPUTS:
	# [v;a]

	R_E = 6378.137 #km

	r = x[1:3]*R_E #km
	v = x[4:6] #km/s
	ω = x[7:9] # rad/s
	q = x[10:13] # unitless 
	ρ = x[14:16]
	MJD = x[17]

	m = u[1:3] 
	ρ̇ = u[4:6]

	#force balance
	#gravity
	#2-body gravity
	mu = 3.986E5 #km3/s2
	f_grav=mu/(norm(r)^2)*-r/(norm(r)) #km/s^2

	# atmospheric drag (assuming that time is invariant)
	# date is new years on 2019

	#find the latitude and longitude in a geometric frame
	t =  58865.5  - 51544.5 # fix for now
	GMST = (280.4606 + 360.9856473*(t))./180*pi #calculates Greenwich Mean Standard Time
	ROT = [cos(GMST) sin(GMST) 0;-sin(GMST) cos(GMST) 0; 0 0 1]
	pos_ecef = ROT*r
	lat= asin(pos_ecef[3]/norm(pos_ecef))
	long= atan(pos_ecef[2],pos_ecef[1])

	# exponential model of atmosphere
	den = SatelliteToolbox.expatmosphere((norm(r)-R_E)*1000)

	# #nrmlsise00
	# out = SatelliteToolbox.nrlmsise00(DatetoJD(2019, 1, 1, 00, 00, 00), # Julian Day
	#                  (norm(r)-R_E)*1000,# Altitude [m]
	#                  -lat,# Latitude [rad]
	#                  -long,# Longitude [rad]
	#                  83.7,# 81 day average F10.7 flux
	#                  102.5,# Daily F10.7 for previous day
	#                  15;# Magnetic index (daily)
	#                  output_si = true# Output in cm^-3 and g/cm^-3
	#     );
	# den = out.den_Total
	ω_earth = [0;0;7.2921150E-5] #rad/s
	v_ecef = v-cross(ω_earth,r)
	f_drag_ecef = .5*den*norm(v_ecef)^2*sat.mass/sat.BC*-v_ecef/norm(v_ecef)
	f_drag = ROT'*f_drag_ecef

	#J2 oblateness term
	J2=0.0010826359
	# J2=0
	f_J2=[J2*r[1]/norm(r)^7*(6*r[3]-1.5*(r[1]^2+r[2]^2))
	     J2*r[2]/norm(r)^7*(6*r[3]-1.5*(r[1]^2+r[2]^2))
	     J2*r[3]/norm(r)^7*(3*r[3]-4.5*(r[1]^2+r[2]^2))]


	#acceleration
	a=(f_grav+f_J2+f_drag) #km/s^2

	# magnetic field
	NED_to_ENU = [0 1 0;1 0 0;0 0 -1]
	R_ENU_to_XYZ =
            [-sin(long) -sin(lat)*cos(long) cos(lat)*cos(long);
            cos(long) -sin(lat)*sin(long) cos(lat)*sin(long);
            0 cos(lat) sin(lat)]
	B_ECI = Rz(GMST)'* R_ENU_to_XYZ * NED_to_ENU * igrf12(MJD_to_year(MJD),norm(r)*1E3, lat, long)*1E-9
	
	#rotation (Euler)
	ω̇ = inv(sat.J)*(cross(m,B_ECI) - ρ̇ - cross(ω,sat.J*ω+ρ))

	#quaternion 
	q̇ = 1/2*(Quaternion(q)*Quaternion([0;ω]))[:]

	return [v/R_E;a;ω̇;q̇;ρ̇;0]

end

function posvel_dynamics_dot(x,u)
	# INPUTS:
	# STATE:
	# [r_eci/R_E;
	# v_eci]
	# CONTROLS:
	# dummy

	# OUTPUTS:
	# [v;a]

	R_E = 6378.137 #km

	r = x[1:3]*R_E #km
	v = x[4:6] #km/s

	#force balance
	#gravity
	#2-body gravity
	mu = 3.986E5 #km3/s2
	f_grav=mu/(norm(r)^2)*-r/(norm(r)) #km/s^2

	# atmospheric drag (assuming that time is invariant)
	# date is new years on 2019

	#find the latitude and longitude in a geometric frame
	t = 58865.5 - 51544.5 # fix MJD for now
	GMST = (280.4606 + 360.9856473*(t))./180*pi #calculates Greenwich Mean Standard Time
	ROT = [cos(GMST) sin(GMST) 0;-sin(GMST) cos(GMST) 0; 0 0 1]
	pos_ecef = ROT*r
	lat= asin(pos_ecef[3]/norm(pos_ecef))
	long= atan(pos_ecef[2],pos_ecef[1])

	# exponential model of atmosphere
	den = SatelliteToolbox.expatmosphere((norm(r)-R_E)*1000)

	# #nrmlsise00
	# out = SatelliteToolbox.nrlmsise00(DatetoJD(2019, 1, 1, 00, 00, 00), # Julian Day
	#                  (norm(r)-R_E)*1000,# Altitude [m]
	#                  -lat,# Latitude [rad]
	#                  -long,# Longitude [rad]
	#                  83.7,# 81 day average F10.7 flux
	#                  102.5,# Daily F10.7 for previous day
	#                  15;# Magnetic index (daily)
	#                  output_si = true# Output in cm^-3 and g/cm^-3
	#     );
	# den = out.den_Total
	ω_earth = [0;0;7.2921150E-5] #rad/s
	v_ecef = v-cross(ω_earth,r)
	f_drag_ecef = .5*den*norm(v_ecef)^2*sat.mass/sat.BC*-v_ecef/norm(v_ecef)
	f_drag = ROT'*f_drag_ecef

	#J2 oblateness term
	J2=0.0010826359
	# J2=0
	f_J2=[J2*r[1]/norm(r)^7*(6*r[3]-1.5*(r[1]^2+r[2]^2))
	     J2*r[2]/norm(r)^7*(6*r[3]-1.5*(r[1]^2+r[2]^2))
	     J2*r[3]/norm(r)^7*(3*r[3]-4.5*(r[1]^2+r[2]^2))]

	#acceleration
	a=(f_grav+f_J2+f_drag) #km/s^2

	return [v/R_E;a]
end

function rk4(x, u, xdot, dt)
	# INPUTS:
	# x: state 
	# u: controls
	# xdot: state dot FUNCTION
	# dt: time step 

	# OUTPUTS:
	# x_new: new state

	k1 = dt*xdot(x,u);
	k2 = dt*xdot(x+k1/2,u);
	k3 = dt*xdot(x+k2/2,u);
	k4 = dt*xdot(x+k3,u);

	return x + 1/6*(k1+2*k2+2*k3+k4)

end

function euler(x,u,xdot,dt)
	# INPUTS:
	# x: state 
	# u: controls
	# xdot: state dot FUNCTION
	# dt: time step 

	# OUTPUTS:
	# x_new: new state

	return x + xdot(x,u)*dt

end


function MJD_to_year(MJD)
	# simple function to convert MJD to what year it is 

	return round(MJD / 365.25 + 1858.8767)

end

function Rx(theta)

    return [1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]

end

function Rz(theta)

    return [cos(theta) sin(theta) 0 ; -sin(theta) cos(theta) 0;0 0 1]

end

