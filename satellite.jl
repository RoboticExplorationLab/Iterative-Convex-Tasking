## this function is to define satellite characteristics

struct satellite
	mass::Float64
	J::AbstractArray{Float64,2}
	BC::Float64 
end

function define_satellite(mass, J, BC)
	# simple function that creates a satellite struct
	# INPUTS:
	# mass (kg)
	# Moment of inertia (kgm2)
	# Ballistic Coefficient (kg/m2)

	return satellite(mass,J,BC)

end