function Base.isapprox(n1::Nothing, n2::Nothing)
    return true
end

abstract type Object end;

struct Sphere <: Object
    sdf::Function
    texture::Function
end

function initialize(t::Type{Sphere}; radius = 1, position = (0,0,0))

    function texture(x, y, z) 
        return (1,1,1)
    end

    return initialize(t, texture; radius = radius, position = position)
end

function initialize(t::Type{Sphere}, texture; radius = 1, position = (0,0,0))
    sdf(x, y, z) = sqrt((x-position[1])^2 +
                        (y-position[2])^2 +
                        (z-position[3])^2) - radius

    return Sphere(sdf, texture)
end

struct Rays
    # Velocity direction vector
    l::T where T <: Union{Array{Float64,3}, CuArray{Float64,3}}

    # Position vector
    p::T where T <: Union{Array{Float64,3}, CuArray{Float64,3}}

    # Color
    c::T where T <: Union{Array{Float64,3}, CuArray{Float64,3}}
end


function propagate!(rays::Rays, objects::Vector{O};
                    numcores = 4, numthreads = 256,
                    num_intersections = 1) where {O <: Object}

    sdf_tuple = Tuple(i.sdf for i in objects)
    texture_tuple = Tuple(i.texture for i in objects)

    println(sdf_tuple)
    println(texture_tuple)

    if isa(rays.p, Array)
        kernel! = propagate_kernel!(CPU(), numcores)
    else
        kernel! = propagate_kernel!(CUDADevice(), numthreads)
    end

    kernel!(rays.p, rays.l, rays.c, num_intersections,
            sdf_tuple, texture_tuple, 
            ndrange=(size(rays.p)[1], size(rays.p)[2]))


end

@kernel function propagate_kernel!(ray_pos, ray_dir, ray_clr, 
                                   sdf_tuple, texture_tuple,
                                   num_intersections)
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
            cam.focal_length, ndrange=(size(rays.p)[1], size(rays.p)[2]))
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

function ray_march(objects::Vector{O}, cam::Camera, rays;
                  filename="check.png",
                  num_intersections = 10, AT = Array) where {O <: Object}

    CUDA.@time wait(init_rays!(rays, cam))

    CUDA.@time wait(propagate!(rays, objects;
                    num_intersections = num_intersections))

    @time convert_to_img(Array(rays.c), cam, filename; AT = AT)

end
