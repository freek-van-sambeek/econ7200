# Load packages
using Plots
using ProgressMeter


# Set up the plotting canvas
gr()

# Define file logistics
PATH_OUTPUT = "section_3/"

# Initialize parameter values
alpha = 0.3
beta = 0.8
# rho = 1 / 0.8 - 1 (This is wrong, right?)
rho = (1 - beta) / beta
delta = 0.2
B = (rho + delta) / alpha
max_iters = 1000
kss = (1 / (alpha * beta * B) + (delta - 1) / (alpha * B))^(1 / (alpha - 1)) # K in the steady state
K_cardinality = 501 # Cardinality of the set K
tolerance = 10^-6 # Tolerance of discrepancy nth value function and true value function

#= Fill the grid of values of starting capital containing
the steady state level of capital =#
K = zeros(Float64, (1, K_cardinality))
for i in 1:K_cardinality
   K[i] = (1 + 0.25 / floor(K_cardinality / 2) * (i - ceil(K_cardinality / 2))) * kss
end

# Initial Value Function
vfun = zeros(Float64, (max_iters, K_cardinality))
gfun = copy(vfun)

# Set up a progressbar for the optimization algorithm
pbar = Progress(max_iters)

# Run the optimization algorithm
i = 1
diff = 10^10

while i <= max_iters && diff > tolerance
   global i = i + 1
   for col_k in 1:K_cardinality
      vfun[i, col_k] = -10^10
      for col_k_prime in 1:K_cardinality
         # Define a helper value to set the value function to if consumption is negative
         vhelp = -10^10

         # Compute consumption as the production minus the investment
         consumption = B * K[col_k]^alpha + (1 - delta) * K[col_k] - K[col_k_prime]
         if consumption <= 0.0
            vhelp = -10^12
         else
            # Compute the new guess for the value function
            vhelp = log(consumption) + beta * vfun[i-1, col_k_prime]
         end

         #= Check if the new guess for the value function is a valid one
         (i.e.) not one with negative consumption=#
         if vhelp >= vfun[i, col_k]
            # If so, log it as the new value function guess and policy function
            vfun[i, col_k] = vhelp
            gfun[i, col_k] = K[col_k_prime]
         end
      end
   end

   # Compute the supremum over all k in K of the absolute difference with the previous value function
   global diff = maximum(abs.(vfun[i, :] - vfun[i-1, :]))
   next!(
      pbar; showvalues=[(string("Sup abs difference with previous value function (with tolerance of ", tolerance, "): "), diff)]
   )
end

# Convert K to an x-axis series
K = vec(K)

# Plot various iterations of the value function over the grid of K
plot(K, permutedims(vfun[[1, 2, 3, 11, i], :]), label=["v0" "v1" "v2" "v10" "Converged v"])
xlabel!("Capital Stock k Today")
ylabel!("Value Function")
title!("Value Function: True and Approximated")
xlims!(K[1], K[end])
ylims!(0, 2)
savefig(string(PATH_OUTPUT, "value_function_iteration.pdf"))

# Plot various iterations of the policy function over the grid of K
plot(K, permutedims(gfun[[2, 3, 11, i], :]), label=["g2" "g3" "g11" "Converged g"])
xlabel!("Capital Stock k Today")
ylabel!("Policy Function")
title!("Policy Function: True and Approximated")
savefig(string(PATH_OUTPUT, "policy_function.pdf"))

#= Plot capital over time assuming we start at capital below the steady state,
and using the converged policy function =#
k0 = 0.9 * kss
mink = K[1]
maxk = K[end]
consumption = zeros(Float64, (101, 1)) # Consumption
capital = zeros(Float64, (102, 1)) # Capital
capital[1] = k0
for j in 1:101
   capital[j + 1] = gfun[i, convert(Int64, floor((capital[j] - mink) / (0.25 / floor(K_cardinality / 2) * kss) + 1))]
   consumption[j] = B * capital[j]^alpha + (1 - delta) * capital[j] - capital[j + 1]
end

t = range(start=0, stop=100, length=101)

plot(t, consumption)
xlabel!("Time")
ylabel!("c(t)")
title!("c(t) over time")
savefig(string(PATH_OUTPUT, "evolution_consumption.pdf"))

plot(t, capital[2:102])
xlabel!("Time")
ylabel!("k(t)")
title!("k(t) over time")
savefig(string(PATH_OUTPUT, "evolution_capital.pdf"))
