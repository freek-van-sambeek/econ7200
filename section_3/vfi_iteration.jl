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
rho = 1 / 0.8 - 1
# rho = (1 - beta) / beta # (Shouldn't rho be this?)
delta = 0.2
B = (rho + delta) / alpha
max_iters = 1000
kss = (1 / (alpha * beta * B) + (delta - 1) / (alpha * B))^(1 / (alpha - 1)) # K in the steady state
K_cardinality = 501 # Cardinality of the set K
tolerance = 10^-6 # Tolerance of discrepancy nth value function and true value function

#= Fill the grid of values of starting capital containing
the steady state level of capital =#
K = [(1 + 0.25 / floor(K_cardinality / 2) * (i - ceil(K_cardinality / 2))) * kss for i in 1:K_cardinality]


function run_vfi()
   # Initial Value Function
   vfun = zeros(Float64, (max_iters + 1, K_cardinality))
   gfun = copy(vfun)

   # Set up a progressbar for the optimization algorithm
   pbar = Progress(max_iters)

   # Run the optimization algorithm
   i = 1
   diff = 10^10

   while i <= max_iters && diff > tolerance
      i = i + 1
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
      diff = maximum(abs.(vfun[i, :] - vfun[i-1, :]))
      next!(
         pbar; showvalues=[(string("Sup abs difference with previous value function (with tolerance of ", tolerance, "): "), diff)]
      )
   end

   # Set the progress bar to complete
   finish!(pbar)

   return vfun[1:i, :], gfun[1:i, :]
end


function plot_results(; vfun::Matrix{Float64}, gfun::Matrix{Float64})
   # Convert K to an x-axis series
   global K = vec(K)

   # Plot various iterations of the value function over the grid of K
   plot(K, permutedims(vfun[[1, 2, 3, 11, end], :]), label=["v0" "v1" "v2" "v10" "Converged v"])
   xlabel!("Capital Stock k Today")
   ylabel!("Value Function")
   title!("Value Function: True and Approximated")
   xlims!(K[1], K[end])
   ylims!(0, 2)
   savefig(string(PATH_OUTPUT, "value_function_iteration.pdf"))

   # Plot various iterations of the policy function over the grid of K
   plot(K, permutedims(gfun[[2, 3, 11, end], :]), label=["g2" "g3" "g11" "Converged g"])
   xlabel!("Capital Stock k Today")
   ylabel!("Policy Function")
   title!("Policy Function: True and Approximated")
   savefig(string(PATH_OUTPUT, "policy_function.pdf"))

   #= Plot capital over time assuming we start at capital below the steady state,
   and using the converged policy function =#
   k0 = 0.9 * kss
   mink = K[1]
   consumption = zeros(Float64, (101, 1)) # Consumption
   capital = zeros(Float64, (102, 1)) # Capital
   capital[1] = k0
   for i in 1:101
      capital[i+1] = gfun[end, convert(Int64, floor((capital[i] - mink) / (0.25 / floor(K_cardinality / 2) * kss) + 1))]
      consumption[i] = B * capital[i]^alpha + (1 - delta) * capital[i] - capital[i+1]
   end

   # Plot consumption and capital over 101 time periods (t=0 to t=100)
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
end


function main()
   vfun, gfun = run_vfi()
   plot_results(vfun=vfun, gfun=gfun)
end


# Only call main when it is not imported
if abspath(PROGRAM_FILE) == @__FILE__
   main()
end
