## this function is to define satellite characteristics

struct satellite
	mass::Float64
	J::AbstractArray{Float64,2}
end

function define_satellite(mass, J)
	# simple function that creates a satellite struct

	return satellite(mass,J)

end