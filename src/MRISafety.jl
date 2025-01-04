
module MRISafety
	
	using FFTW

	"""
		Not sure what effect δf has

	"""
	function gradient_power(g::AbstractVector{<: Real}, δt::Real, δf::Real, forbidden::NTuple{N, NTuple{2, Real}}) where N
		max_time = floor(Int, 1 / (δf * δt))
		max_time += mod(max_time, 2)
		@assert 0 < max_time ≤ length(g)
		g = @view g[1:max_time]
		power_spectrum = abs2.(rfft(g))
		frequency = rfftfreq(max_time, 1/δt)
		forbidden_power = Vector{Float64}(undef, N)
		total_power = sum(power_spectrum)
		for (i, (f, Δf)) in enumerate(forbidden)
			idx = (f - 0.5Δf .< frequency .< f + 0.5Δf)
			forbidden_power[i] = sum(power_spectrum[idx]) / total_power
		end
		return forbidden_power, power_spectrum, frequency
	end
end

