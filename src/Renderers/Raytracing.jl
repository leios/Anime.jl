function Base.isapprox(n1::Nothing, n2::Nothing)
    return true
end

struct Rays
    # Velocity direction vector
    l::T where T <: Union{Array{Float64,3}, CuArray{Float64,3}}

    # Position vector
    p::T where T <: Union{Array{Float64,3}, CuArray{Float64,3}}

    # Color
    c::T where T <: Union{Array{Float64,3}, CuArray{Float64,3}}

end

struct Surface

    # Reflectivity
    r::Float64

    # Transmission
    t::Float64

    # Color
    c::RGBA

    # index of refraction
    ior::Float64

    function Surface(in_r, in_t, in_c, in_ior)
        if !isapprox(in_r+in_t+in_c.alpha, 1)
            error("invalid surface definition, RTC < 1")
        end
        new(in_r,in_t,in_c, in_ior)
    end

    Surface(in_r, in_t, in_c::Float64, in_ior) =
         new(in_r, in_t, RGBA(0,0,0,0), in_ior)
end

abstract type Object end

struct Sphere <: Object
    # Lens position
    p::Vector{Float64}

    # Lens radius
    r::Float64

    s::Surface
end

function Lens(p, r, ior)
    return Sphere(p, r, Surface(0,1,RGBA(0,0,0,0),ior))
end

function ReflectingSphere(p, r)
    return Sphere(p,r,Surface(1,0,RGBA(0,0,0,0),0))
end

function ColoredSphere(p, r, c::RGB)
    return Sphere(p, r, Surface(0,0,RGBA(c), 0))
end

struct SkyBox <: Object
    # Skybox position
    p::Vector{Float64}

    # Skybox radius
    r::Float64
end

# NOTE: rename without "sphere" so it can be used for other objects
function sphere_normal_at(ray_pos, sphere)
    n = normalize(ray_pos .- sphere.p)

    return n
end

function inside_of(pos, sphere)
    x = sphere.p[1] - pos[1]
    y = sphere.p[2] - pos[2]
    if (x^2 + y^2 <= sphere.r^2)
        return true
    else
        return false
    end
end

# NOTE: light moves as a particular speed with respect to the medium it is
#       moving through, so...
#       n_2*v = n_1*l + (n_1*cos(theta_1) - n_2*cos(theta_2))*n
#       Other approximations: ior = n_1/n_2, c = -n*l
function refract!(ray_pos, ray_dir, lens::Sphere, ior)
    n = sphere_normal_at(ray_pos, lens)

    if dot(n, ray_dir) > 0
        n .*= -1
    end
    c = dot(-n, ray_dir);
    d = 1.0 - ior^2 * (1.0 - c^2);

    if (d < 0.0)
        reflect!(ray_dir, ray_pos, n)
    end

    ray_dir = ior * raydir + (ior * c - sqrt(d)) * n;
end

abstract type Wall <: Object end;

struct Mirror <: Wall
    # Normal vector
    n::Vector{Float64}

    # Position of mirror
    p::Vector{Float64}

    # Mirror size
    scale::Float64

    Mirror(in_n, in_p) = new(in_n, in_p, 2.5)
end

# NOTE: Probably remove?
function is_behind(ray, mirror)
    if dot(ray.p.-mirror.p, mirror.n) >= 0
        return true
    else
        return false
    end
end

# note: for reflection, l_x -> l_x, but l_y -> -l_y
#       In this case, the y component of l = cos(theta)*n
#       so new vector: v = l + 2cos(theta)*n
function reflect!(ray_pos, ray_dir, n)
    ray_dir = ray_dir .- 2*dot(ray_dir, n).*n
    ray_pos = ray_pos .+ 0.001*ray_dir
end

function intersection(ray_pos, ray_dir, sphere::S;
                      threshold = 0.01) where
                      {S <: Union{Sphere, SkyBox}}
    relative_dist = ray_pos-sphere.p
    a = dot(ray_dir, ray_dir)
    b = 2.0 * dot(relative_dist, ray_dir)
    c = dot(relative_dist, relative_dist) - sphere.r*sphere.r
    discriminant = b*b - 4*a*c

    if discriminant < 0
        return nothing
    elseif discriminant > 0
        roots = [(-b + sqrt(discriminant)) / (2*a),
                 (-b - sqrt(discriminant)) / (2*a)]
        min = minimum(roots)
        max = maximum(roots)

        if min > threshold
            return (min)*ray_dir
        elseif max > threshold
            return (max)*ray_dir
        else
            return nothing
        end
    else
        # Returns nothing if tangential
        return nothing
        #return (-b/(2*a))*ray_dir
    end 
end


function intersection(ray_pos, ray_dir, wall::W) where {W <: Wall}
    intersection_pt = -dot((ray_pos .- wall.p),wall.n)/dot(ray_dir, wall.n)

    if isfinite(intersection_pt) && intersection_pt > 0 &&
       intersection_pt != NaN
        return intersection_pt*ray_dir
    else
        return nothing
    end
end

function propagate!(rays::Rays, objects::Vector{O},
                    num_intersections) where {O <: Object}
    @floop ThreadedEx() for i = 1:size(rays.p)[1]
        for j = 1:size(rays.p)[2]
            output = propagate!(rays.p[i,j,:], rays.l[i,j,:], rays.c[i,j,:],
                                objects, num_intersections)
            rays.p[i,j,:] .= output[1][:]
            rays.l[i,j,:] .= output[2][:]
            rays.c[i,j,:] .= output[3][:]
        end
    end

    return rays
end

function propagate!(ray_pos, ray_dir, ray_clr, objects::Vector{O},
                   num_intersections) where {O <: Object}

    for i = 1:num_intersections
        if ray_dir != zeros(length(ray_dir))
            intersect_final = [Inf, Inf]
            intersected_object = nothing
            for object in objects
                intersect = intersection(ray_pos, ray_dir, object)
                if intersect != nothing &&
                   sum(intersect[:].^2) < sum(intersect_final[:].^2)
                    intersect_final = intersect
                    intersected_object = object
                end
            end

            if intersect_final != nothing
                #ray = Ray(ray.l, ray.p .+ intersect_final, ray.c)
                if typeof(intersected_object) == Sphere
                    #reflected_ray = ray
                    #refracted_ray = ray
                    #colored_ray = ray
                    if !isapprox(intersected_object.s.t, 0)
                        ior = 1/intersected_object.s.ior
                        if dot(ray_dir,
                               sphere_normal_at(ray_pos, ray_dir,
                                                intersected_object)) > 0
                            ior = intersected_object.s.ior
                        end

                        refract!(ray_pos, ray_direction,
                                 intersected_object, ior)
                        propagate!(ray_pos, ray_dir, ray_clr, objects,
                                   num_intersections-i)
                    end

                    if !isapprox(intersected_object.s.r, 0)
                        n = sphere_normal_at(ray_pos, ray_dir,
                                             intersected_object)
                        reflect(ray_pos, ray_dir, n)
                        propagate!(ray_pos, ray_dir, ray_clr, objects,
                                   num_intersections-i)
                    end

                    if !isapprox(intersected_object.s.c.alpha, 0)
                        ray_color = RGB(intersected_object.s.c)
                        ray_dir .= zeros(length(ray_dir))
                        ray_clr[1] = ray_color.r
                        ray_clr[2] = ray_color.g
                        ray_clr[3] = ray_color.b
                    end

#=
                    ray_color = intersected_object.s.t*refracted_ray.c +
                                intersected_object.s.r*reflected_ray.c +
                                intersected_object.s.c.alpha*colored_ray.c

                    ray = Ray(zeros(length(ray.l)), ray.p, ray_color)
=#

                elseif typeof(intersected_object) == Mirror
                    reflect!(ray_pos, ray_dir, intersected_object.n)
                elseif typeof(intersected_object) == SkyBox
                    temp_clr = pixel_color(ray_pos, 100)
                    ray_clr[1] = temp_clr.r
                    ray_clr[2] = temp_clr.g
                    ray_clr[3] = temp_clr.b
                    ray_dir .= zeros(length(ray_dir))
                end
            else
                println("hit nothing")
            end
        end
    end

    return [ray_pos, ray_dir, ray_clr]

end

function pixel_color(position, extents)
    c = RGB(0)
    if position[1] < extents && position[1] > -extents
        c += RGB((position[1]+extents)/(2.0*extents), 0, 0)
    else
        println(position)
    end

    if position[2] < extents && position[2] > -extents
        c += RGB(0,0,(position[2]+extents)/(2.0*extents))
    else
        println(position)
    end

    if position[3] < extents && position[3] > -extents
        c += RGB(0,(position[3]+extents)/(2.0*extents), 0)
    else
        println(position)
    end

    return c
end

# NOTE: Speed this up!
function convert_to_img(ray_colors::Array{Float64}, cam::Camera, filename;
                        AT = Array)
    for i = 1:size(cam.pixels)[1]
        for j = 1:size(cam.pixels)[2]
            cam.pixels[i,j] = RGB(ray_colors[i,j,1],
                                  ray_colors[i,j,2],
                                  ray_colors[i,j,3])
        end
    end

    save(filename, transpose(cam.pixels))
end

# NOTE: Extra allocations here when creating the appropriate camera arrays
#       to pass to  GPU kernel
function init_rays!(rays::Rays, cam::Camera; numcores = 4, numthreads = 256)
    AT = Array
    if isa(rays.p, Array)
        kernel! = init_rays_kernel!(CPU(), numcores)
    else
        kernel! = init_rays_kernel!(CUDADevice(), numthreads)
        AT = CuArray
    end

    kernel!(rays.p, rays.l, rays.c,
            AT(cam.p), AT([size(cam.pixels)[1], size(cam.pixels)[2]]),
            AT([cam.size[1], cam.size[2]]),
            cam.focal_length, ndrange=size(rays.p))
end

@kernel function init_rays_kernel!(ray_positions, ray_directions, ray_colors,
                                   cam_p, res, dim, focal_length)
    i,j = @index(Global, NTuple)
    pixel_width_x = dim[1] / res[1]
    pixel_width_y = dim[2] / res[2]

    # create a set of rays that go through every pixel in our grid.
    ray_positions[i,j,1] = cam_p[1] + 0.5*dim[1] - i*dim[1]/res[1] + 
                          0.5*pixel_width_x
    ray_positions[i,j,2] = cam_p[2] + 0.5*dim[2] - j*dim[2]/res[2] +
                          0.5*pixel_width_y
    ray_positions[i,j,3] = cam_p[3]+focal_length

    # This is for normalization
    temp_sum = 0
    for k = 1:3
        ray_colors[i,j,k] = 0

        ray_directions[i,j,k] = ray_positions[i,j,k] - cam_p[k]
        temp_sum += ray_directions[i,j,k]^2
    end

    temp_sum = sqrt(temp_sum)

    for k = 1:3
        ray_directions[i,j,k] /= temp_sum
    end

end

function ray_trace(objects::Vector{O}, cam::Camera, rays;
                   filename="check.png",
                   num_intersections = 10, AT = Array) where {O <: Object}

    CUDA.@time wait(init_rays!(rays, cam))

    @time rays = propagate!(rays, objects, num_intersections)

    @time convert_to_img(Array(rays.c), cam, filename; AT = AT)

end
