using Distributions

# Degrees of freedom for the t-distribution
ν = 100

# Desired mean and standard deviation
μ = -2.0
σ = 3.0

# Create a standard Student's t-distribution
t_dist = TDist(ν)

# To modify the distribution: multiply by σ (for std dev) and add μ (for mean)
# Shifted and scaled t-distribution: X = μ + σ * T
shifted_scaled_cdf = x -> cdf(t_dist, (x - μ) / σ)

# Calculate the CDF at some point x
x = 1.5
cdf_value = shifted_scaled_cdf(x)
print("tail: ",1-cdf_value)
