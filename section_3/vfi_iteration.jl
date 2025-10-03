# Load packages
using DataFrames
using Plots
using ProgressMeter


# Set up the plotting canvas
gr()

# Define file logistics
PATH_OUTPUT = "section_3/"

# Define global variables
max_iters = 1000
tolerance = 10^-6 # Tolerance of discrepancy nth value function and true value function


function return_grid(; kss::Float64, K_cardinality::Int64)
   return [(1 / 2 + (i -  0.3) / K_cardinality) * kss for i in 1:K_cardinality]
end


"""
   run_vfi() -> Tuple{Matrix{Float64}, Matrix{Float64}}

Estimate the value and policy functions and returns the grid-evaluated functions
at the different iterations.
"""
function run_vfi(; alpha::Float64, beta::Float64, delta::Float64, A::Float64, K::Array, u::Function)
   # Initial Value Function
   K_cardinality = length(K)
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
            consumption = A * K[col_k]^alpha + (1 - delta) * K[col_k] - K[col_k_prime]
            if consumption <= 0.0
               vhelp = -10^12
            else
               # Compute the new guess for the value function
               vhelp = u(consumption) + beta * vfun[i-1, col_k_prime]
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


"""
   plot_results(; vfun::Matrix{Float64}, gfun::Matrix{Float64}) -> Nothing

Plot iterations of the value and policy functions.

# Arguments
- `vfun::Matrix{Float64}`: The matrix with in the iterations indexing the rows,
   and the capital grid index indexing the columns of the value function.

# Returns
- `gfun::Matrix{Float64}`: The matrix with in the iterations indexing the rows,
   and the capital grid index indexing the columns of the policy function.
"""
function plot_results(; alpha::Float64, beta::Float64, A::Float64, delta::Float64, kss::Float64, K::Vector{Float64}, vfun::Matrix{Float64}, gfun::Matrix{Float64}, model::String)
   # # Convert K to an x-axis series
   # K = vec(K)

   # Obtain the cardinality of K again
   K_cardinality = length(K)

   # Compute the analytical solution
   a_1 = alpha / (1 - alpha * beta)
   a_0 = ((1 + beta * a_1) * (log(A) - log(1 + beta * a_1)) + beta * a_1 * log(beta * a_1)) / (1 - beta)
   # a0 = (log(B * (1 - alpha * beta)) + alpha * beta / (1 - alpha * beta) * log(alpha * beta * B)) / (1 - beta)
   v_analytical = [(a_0 + a_1 * log(K[i])) for i in 1:K_cardinality]
   g_analytical = [(alpha * beta * A * K[i]^alpha) for i in 1:K_cardinality]

   # Plot various iterations of the value function over the grid of K
   plot(K, permutedims(vfun[[1, 2, 3, 11, end], :]), label=["v0" "v1" "v2" "v10" "Converged v"])
   if !occursin("sigma", model)
      plot!(K, v_analytical, label="Analytical v", ls=:dash, color=:black)
   end
   xlabel!("Capital Stock k Today")
   ylabel!("Value Function")
   title!("Value Function: True and Approximated")
   xlims!(K[1], K[end])
   # ylims!(0, 2)
   savefig(string(PATH_OUTPUT, "value_function_iteration_", model, ".pdf"))
   
   # Plot various iterations of the policy function over the grid of K
   plot(K, permutedims(gfun[[2, 3, 11, end], :]), label=["g2" "g3" "g11" "Converged g"])
   if !occursin("sigma", model)
      plot!(K, g_analytical, label="Analytical g", ls=:dash, color=:black)
   end
   xlabel!("Capital Stock k Today")
   ylabel!("Policy Function")
   title!("Policy Function: True and Approximated")
   savefig(string(PATH_OUTPUT, "policy_function_", model, ".pdf"))
   
   #= Plot capital over time assuming we start at capital below the steady state,
   and using the converged policy function =#
   k0 = 0.8 * kss
   min_k = K[1] / kss
   max_k = K[end] / kss
   consumption = zeros(Float64, 101) # Consumption
   capital = zeros(Float64, 102) # Capital
   capital[1] = k0
   # println(K)
   # println(gfun[end, :])
   for i in 1:101
      # println(convert(Int64, floor((capital[i] / kss - min_k) / (max_k - min_k) * K_cardinality) + 1))
      capital[i+1] = gfun[end, convert(Int64, floor((capital[i] / kss - min_k) / (max_k - min_k) * K_cardinality) + 1)]
      consumption[i] = A * capital[i]^alpha + (1 - delta) * capital[i] - capital[i+1]
   end

   # Print a table with the capital and consumption paths
   df = DataFrame(Capital = capital[2:end], Consumption = consumption)
   println(df)

   # Plot consumption and capital over 101 time periods (t=0 to t=100)
   t = range(start=0, stop=100, length=101)

   plot(t, consumption)
   xlabel!("Time")
   ylabel!("c(t)")
   title!("c(t) over time")
   savefig(string(PATH_OUTPUT, "evolution_consumption_", model, ".pdf"))

   plot(t, capital[2:102])
   xlabel!("Time")
   ylabel!("k(t)")
   title!("k(t) over time")
   savefig(string(PATH_OUTPUT, "evolution_capital_",  model, ".pdf"))
end


function main()
   ### Q3 - Q4 ### - Run VFI
   # Initialize parameter values
   alpha = 0.3
   beta = 0.6
   rho = 1 / beta - 1
   delta = 1.0
   A = (rho + delta) / alpha

   # Compute the steady state capital
   kss = (rho / (1 - rho))^(1 / (1 - alpha))

   #= Fill the grid of values of starting capital containing
   the steady state level of capital =#
   K = return_grid(kss=kss, K_cardinality=5)
   
   # Run the VFI
   vfun, gfun = run_vfi(alpha=alpha, beta=beta, delta=delta, A=A, K=K, u=x->log(x))
   
   # Plot the results
   plot_results(alpha=alpha, beta=beta, A=A, delta=delta, kss=kss, K=K, vfun=vfun, gfun=gfun, model="coarse")
   
   ### Q5 - Q6 ### - Run the VFI with a different, more dense grid of values
   K = return_grid(kss=kss, K_cardinality=200)

   # Run the VFI
   vfun, gfun = run_vfi(alpha=alpha, beta=beta, delta=delta, A=A, K=K, u=x->log(x))
   
   # Plot the results
   plot_results(alpha=alpha, beta=beta, A=A, delta=delta, kss=kss, K=K, vfun=vfun, gfun=gfun, model="fine")
   
   ### Q7 ###
   # Adjust the starting parameter values
   beta = 0.8
   rho = 1 / beta - 1
   A = (rho + 0.2) / alpha

   # Run the VFI
   vfun, gfun = run_vfi(alpha=alpha, beta=beta, delta=delta, A=A, K=K, u=x->log(x))
   
   # Plot the results
   plot_results(alpha=alpha, beta=beta, A=A, delta=delta, kss=kss, K=K, vfun=vfun, gfun=gfun, model="fine_alternative")
   
   ### Q8 ### - Change the utility function
   for sigma in [0.5, 2]
      # Run the VFI
      vfun, gfun = run_vfi(alpha=alpha, beta=beta, delta=delta, A=A, K=K, u=x->(x^(1 - sigma) - 1) / (1 - sigma))

      # Plot the results
      plot_results(alpha=alpha, beta=beta, A=A, delta=delta, kss=kss, K=K, vfun=vfun, gfun=gfun, model=string("fine_sigma_", sigma))
   end
end


# Only call main when it is not imported
if abspath(PROGRAM_FILE) == @__FILE__
   main()
end
