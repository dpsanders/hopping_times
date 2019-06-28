using Distributions, LinearAlgebra


normsq(v) = sum(abs2, v)

"""
Generate initial condition for two spheres with radii r1 and r2, in a box of size 2a*2b
"""
function initial_condition(w, h, r)

    a = w/2 - r
    b = h/2 - r

    # TODO: Add check for max number of tries to place discs

    # position:
    x1, x2 = rand(Uniform(-a, a), 2)
    y1, y2 = rand(Uniform(-b, b), 2)

    while (x1 - x2)^2 + (y1 - y2)^2 <= (2r)^2

        x1, x2 = rand(Uniform(-a, a), 2)
        y1, y2 = rand(Uniform(-b, b), 2)
    end
 
   # generate velocities whose sum squared is 1 by rejection method
    v = ones(4)

    while normsq(v) > 1
        v = rand(Uniform(-1.0, 1.0), 4)
    end

    ## normalize is now on LinearAlgebra
    normalize!(v)  # sum squared is 1

    return [x1, y1], v[1:2], [x2, y2], v[3:4]
end

#initial_condition(a, b, r) = initial_condition(a, b, r, r)

initial_condition(1, 2, 0.1)

"""
Simulate the collision of two discs of radii r1 and r2
in the box [-w/2, w/2] × [-h/2, h/2]
"""
function dynamics(w, h, r, num_collisions)

    t = 0.0

    times_data = Float64[]
    position_data = Vector{Float64}[]
    velocity_data = Vector{Float64}[]
    collision_type_data = Int[]

    x1, v1, x2, v2 = initial_condition(w, h, r)
    x = vcat(x1, x2)  # all position coords of both discs
    v = vcat(v1, v2)  # all velocity coords

    # Code for wall collisions (in variable which):
    # 1, 2, 3, 4 correspond to x1, y1, x2, y2 colliding with wall
    # positive corresponds to right or top wall; negative corresponds to left or bottom wall
    # 5 corresponds to two discs colliding

    # set up arrays so that all wall collisions achieved
    # point particle at centre of disc collides with wall moved by radius

    # high_walls = [a - r1, b - r2, a - r1, b - r2]

    a = w/2 - r
    b = h/2 - r

    high_walls = [ a,  b,  a,  b]
    low_walls =  [-a, -b, -a, -b]

    collision_times = zeros(5)

    which = 10  # which obstacle whichly hit

	for n in 1:num_collisions

        # find wall collision times:
        # for each component of each disc, check for collisions with both walls:

        # (c - x) / v  is the time to collide with wall at c in that direction

        for i in 1:4

			t1 = (low_walls[i] - x[i]) / v[i]  # wall at lowest coordinate in that direction (on left or bottom)
			t2 = (high_walls[i] - x[i]) / v[i]  # wall at opposite position
			collision_times[i] = max(t1, t2)
        end


		# find disc collision time:
        # TODO: Use StaticArrays?

        Δx = [ x[1] - x[3], x[2] - x[4] ]
		Δv = [ v[1] - v[3], v[2] - v[4] ]

		a = Δv⋅Δv
		b = 2 * Δx⋅Δv
		c = normsq(Δx) - (2r)^2

        discriminant = b^2 - 4a*c

        collision_times[5] = Inf  # assume no collision

        # the condition `which != 5` excludes the sphere from consideration if was hit at previous step

		if discriminant >= 0 && which != 5
			d = √discriminant

            t1 = (-b + d) / (2a)
		    t2 = (-b - d) / (2a)

			if t1 > 0
				if t2 > 0
					collision_times[5] = min(t1, t2)
				else
					collision_times[5] = t1
                end
            end
        end

		min_time, which = findmin(collision_times)

        t += min_time
		x += v * min_time


		# implement collision:
		if which == 5  # collision of two discs

			Δx = [ x[1] - x[3], x[2] - x[4] ]

            v1 = v[1:2]
			v2 = v[3:4]

			vp1 = ((v1 ⋅ Δx) / normsq(Δx)) * Δx
			vn1 = v1 - vp1

            vp2 = ((v2 ⋅ Δx) / normsq(Δx)) * Δx
			vn2 = v2 - vp2

			v[1:2] = vn1 + vp2
			v[3:4] = vn2 + vp1

        else
			v[which] = -v[which]  # reflect off hard wall

            if x[which] < 0
                which = -which  # indicate which wall collided with using negative sign for low walls
            end

        end


        push!(position_data, copy(x))
        push!(velocity_data, copy(v))
        push!(times_data, t)
        push!(collision_type_data, which)

    end

    return times_data, position_data, velocity_data, collision_type_data
end
