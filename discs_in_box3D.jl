using Distributions, LinearAlgebra


normsq(v) = sum(abs2, v)

"""
Generate initial condition for two spheres with radii r1 and r2,
in a box of size 2a*2b*2c.
Returns positions first disc, velocity first disc, positions second disc
velocities second disc
"""

function initial_condition(h1, h2, h3, r)

    a = h1/2 - r
    b = h2/2 - r
    c = h3/2 - r

   # print(" ea r ini ", r, " ")
    # TODO: Add check for max number of tries to place discs

    # position:
    x1, x2 = rand(Uniform(-a, a), 2)
    y1, y2 = rand(Uniform(-b, b), 2)
    z1, z2 = rand(Uniform(-c, c), 2)
    distsq = normsq([x1-x2, y1-y2, z1-z2])

    n=0
    # repeat until right, more common with more dimensions.
    while distsq <= (2r)^2

        x1, x2 = rand(Uniform(-a, a), 2)
        y1, y2 = rand(Uniform(-b, b), 2)
        z1, z2 = rand(Uniform(-c, c), 2)
        distsq = normsq([x1-x2, y1-y2, z1-z2])
        n+=1
    #    if mod(n, 10^6)==0
    #        print( "van $n intentos de iniciar")
    #        flush(stdout)
    #    end
    end
 
   # generate velocities whose sum squared is 1 by rejection method
    v = ones(6)
    avs1= normsq(v)
    k=0
    println(k)
    while avs1 > 1
        v = rand(Uniform(-1.0, 1.0), 6)
        avs1=normsq(v)
        #=
        k=k+1
        if mod(k,1000)==0
            println( "haaay $k basuras y $avs1")
        end
        =#
    end

    ## normalize is now on LinearAlgebra
    normalize!(v)  # sum squared is 1

    return [x1, y1, z1], v[1:3], [x2, y2, z2], v[4:6]
end

#initial_condition(a, b, r) = initial_condition(a, b, r, r)
#println("hola pipo")
initial_condition(1, 1.2, 1.5, 0.1)

"""
Simulate the collision of two discs of radii r1 and r2
in the box [-w/2, w/2] × [-h/2, h/2]
"""
function dynamics(h1, h2, h3, r, num_collisions)

    t = 0.0

    times_data = Float64[]
    position_data = Vector{Float64}[]
    velocity_data = Vector{Float64}[]
    collision_type_data = Int[]

#    print(" ea r dynamics yeah  ", r, " ")
    
    x1, v1, x2, v2 = initial_condition(h1, h2, h3, r)
    x = vcat(x1, x2)  # all position coords of both discs
    v = vcat(v1, v2)  # all velocity coords

    # Code for wall collisions (in variable which):
    # 1, 2, 3, correspond to first disc colliding with positive h1,h2, h3 direction.
    # 4, 5, 6, correspond to second disc colliding with positive h1,h2, h3 direction.
    # negative correspond to negative direction walls
    # 7 corresponds to two discs colliding

    # set up arrays so that all wall collisions achieved
    # point particle at centre of disc collides with wall moved by radius

    

    a = h1/2 - r
    b = h2/2 - r
    c = h3/2 - r
    
    high_walls = [ a,  b,  c,  a,  b, c]
    low_walls =  [-a, -b, -c,  -a, -b, -c]

    collision_times = zeros(7)

    which = 17  # which obstacle whichly hit

	for n in 1:num_collisions

        # find wall collision times:
        # for each component of each disc, check for collisions with both walls:

        # (c - x) / v  is the time to collide with wall at c in that direction

        for i in 1:6

	    t1 = (low_walls[i] - x[i]) / v[i]  # wall at positive direction
	    t2 = (high_walls[i] - x[i]) / v[i]  # wall at opposite position
	    collision_times[i] = max(t1, t2) # just one time is positive

        end


		# find disc collision time:
        # TODO: Use StaticArrays?

        Δx = x1 - x2
	Δv = v1 - v2

        # chicharronera, David please do not repeat variables for completely other
        # things
	ach = Δv⋅Δv
	bch = 2 * Δx⋅Δv
	cch = normsq(Δx) - (2r)^2
           
        discriminant = bch^2 - 4*ach*cch

        collision_times[7] = Inf  # assume no collision

        # the condition `which != 5` excludes the sphere from consideration if was hit at previous step
     
            
	if discriminant >= 0 && which != 7
	    d = √discriminant

            t1 = (-bch + d) / (2*ach)
            t2 = (-bch - d) / (2*ach)
          
            
	    if t1 > 0
		if t2 > 0
		    collision_times[7] = min(t1, t2)
		else
		    collision_times[7] = t1
                end
            end
            
        end
            

	min_time, which = findmin(collision_times)

        t += min_time
	x += v * min_time


	# implement collision:
	if which == 7  # collision of two discs           
            # interchange the velocitie components along Delta x
            # that is, the vector joining their centers.
            
            vp1 = ((v1 ⋅ Δx) / normsq(Δx)) * Δx
	    vn1 = v1 - vp1

            vp2 = ((v2 ⋅ Δx) / normsq(Δx)) * Δx
	    vn2 = v2 - vp2

            
	    v[1:3] = vn1 + vp2
	    v[4:6] = vn2 + vp1

        else
	    v[which] = -v[which]  # reflect off hard wall

            if x[which] < 0
                which = -which  # indicate which wall collided with using negative sign for low walls
            end

        end

        x1=x[1:3]
        x2=x[4:6]
        v1=v[1:3]
        v2=v[4:6]    
            
        push!(position_data, copy(x))
        push!(velocity_data, copy(v))
        push!(times_data, t)
        push!(collision_type_data, which)

    end # for n collitions

    return times_data, position_data, velocity_data, collision_type_data
end
