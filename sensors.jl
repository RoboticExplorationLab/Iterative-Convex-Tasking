## This file is for sensor modeling

using Distributions
using LinearAlgebra

mutable struct sensor_struct
	μ::Float64 # mean
	Ω # covariance
end

function rel_range(r1::Array,r2::Array,sensor)
	# this function takes in two state positions in ECI (in R_E units)

	if (length(r1) != 3 || length(r2) != 3)
		print("Relative ranges are incorrect length")
	end

	R_E = 6378 # km
	# NOTE: OUTPUTS IN R_E UNITS

	if r1 == r2
		range = 0
	else
		range = norm(r1-r2) + rand(Normal(sensor.μ,sensor.Ω))
	end

	return range

end