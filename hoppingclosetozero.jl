#= 
just hopping, close to limiting radius, along one  axis.
=#

using JLD
include("discs_in_box.jl")

w = 1.5
h = 1.0
r = 0.24


C = 3π/2   # E = 1/2 => abs(v^2) = 1
rhmax = h / 4 #maximal radius in which a horizontal hopp is possible
rvmax= w/4 #maximal radius in which a vertical hopp is possible (any hopp is not posible in our geometry after this)
rmax=  (h+w-√(2*h*w))/2  
#precompile
times, positions, velocities, collision_types = dynamics(w, h, r, 10^2);

"""
Calculate the times at which horizontal hops occur
"""
function horizontal_hopping_times(times, positions, velocities)
    Δxs = [x[3] - x[1] for x in positions]  # x_2 - x_1
    Δus = [v[3] - v[1] for v in velocities]  # u_2 - u_1;
    
    # indices where there is a hop between collisions i and i+1: 
    horiz_hop_indices = findall(i->sign(Δxs[i]) != sign(Δxs[i+1]),1:(length(positions)-1) );  
    
    # x + t*u = 0   so   t = -x/u
    horiz_hopping_times = times[horiz_hop_indices] - (Δxs[horiz_hop_indices] ./ Δus[horiz_hop_indices])
    
    return horiz_hopping_times
end




function simulation_horizontal_hopping(w, h, rs, num_collisions=10^7)
    
    num_times = Float64[]
   
    
    σ = Float64[]
    stderror = Float64[]
  
    rs = rs
 
    for r in rs
        print(r, " <- vamos a hacer ese radio")
    
        a = w/2 - r
        b = h/2 - r
    
        hor_hopping_times=[]
        nhits=0
    
    #    while(nhits<200)
            times, positions, velocities, collision_types = dynamics(w, h, r, num_collisions);
            hor_hopping_times = diff(horizontal_hopping_times(times, positions, velocities))
            nhits=length(hor_hopping_times)
#        end
             
        push!(num_times, mean(hor_hopping_times))   # diff gives inter-hop times
        
        aux=std(hor_hopping_times)
        push!(σ, aux)
        push!(stderror, aux/sqrt(nhits))

        
    end

    
    return rs, num_times, σ, stderror

end

logdeltas=-10:(0.25):-3
diffs=exp.(logdeltas)
## esta version del código no da resultados analiticos

## estos son los radios
rs=rhmax.-diffs


#precompile
nevents = 10^1
(rs, num_hor_hop_data,
 sigma, stderror) = simulation_horizontal_hopping(w, h, rs,nevents);


#en serio
nevents=10^6
(rs, num_hor_hop_data,
 sigma, stderror) = simulation_horizontal_hopping(w, h, rs,nevents);


save("hoppinglongclosetozero_2D.jld",
     "radiolimiteparahopping", rhmax,
    "rs", rs, 
    "num_hor_hop_data", num_hor_hop_data,
     "nevents", nevents,
     "sigma", sigma,
     "stderror", stderror)
    


