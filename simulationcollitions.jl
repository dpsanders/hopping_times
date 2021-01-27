function simulation_disc_collisions(h1, h2, h3, num_collisions=10^5)
    
    num_times = Float64[]
    σ=Float64[]
    stderror=Float64[]

    n_min=num_collisions/500 # minumum number of choosen events to get out of dynamics  
    rs = 0.0051:0.00256:rmax
    
    rs = 0.0051:0.005:rmax
    
    for r in rs
        print(r, " ")
        
        a = h1/2 - r
        b = h2/2 - r
        c = h3/2 - r
        
        times, positions, velocities, collision_types = dynamics(h1, h2, h3, r, num_collisions);
        
        collision_times = diff(times[collision_types .== 7])
        nhits=length(collision_times)
        
        push!(num_times, mean(collision_times))   # diff gives inter-hop times
        push!(exact_times, disc_collision_analytical(a, b, r))
        aux=std(num_times)
        push!(σ, aux)
        push!(stderror, aux/sqrt(nhits))
        
        
    end
    
    return rs, num_times, exact_times, σ, stderror
end
