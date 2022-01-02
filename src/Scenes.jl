function sphere_rand(num_spheres, max_extents, max_radius)
    spheres = [Sphere([0],0,Surface(0,0,0.,0)) for i = 1:num_spheres]
    for i = 1:num_spheres
        r = rand()
        t = rand()
        ca = rand()
        all = r+t+ca
        r = r / all
        t = t / all
        ca = ca / all
        ior = 1 + rand()
        sphere_surface = Surface(r, t, RGBA(rand(), rand(), rand(), ca), ior)
        spheres[i] = Sphere([rand()*(max_extents[1, 2]-max_extents[1, 1])+
                             max_extents[1, 1],
                             rand()*(max_extents[2, 2]-max_extents[2, 1])+
                             max_extents[2, 1],
                             rand()*(max_extents[3, 2]-max_extents[3, 1])+
                             max_extents[3, 1]],
                            rand()*max_radius, sphere_surface)

    end
    return spheres
end

function main(;AT=Array)
    sky = [SkyBox([0.0, 0.0, 0.0], 1000)]
    #spheres = [ColoredSphere([-50,0,-25], 20, RGB(0, 0, 1))]
    #spheres = [Lens([0,0,-25], 20, 1.5), ReflectingSphere([0,50,-100],20),
    #           ColoredSphere([-50,0,-25], 20, RGB(0, 0, 1))]
#=
    spheres = [Lens([0,0,-25], 20, 1.5), ReflectingSphere([0,50,-100],20),
               ColoredSphere([-50,0,-25], 20, RGB(0, 0, 1)),
               Sphere([30, 25, -60], 20,
                      Surface(0.0, 0.75, RGBA(1,0,0,0.25), 1.5)),
               Sphere([50, 0, -25], 20,
                      Surface(0.5, 0.0, RGBA(0,1,0,0.5), 1.5)),
               Sphere([-30, 25, -60], 20,
                      Surface(0.5, 0.5, RGBA(1,1,1,0), 1.5))]
=#
    #bg_spheres = sphere_rand(10, [-400 400; -400 400; -350 -400], 50)


    blank_img = Array{RGB}(undef, 1920, 1080)
    #blank_img = Array{RGB}(undef, 20, 10)
    blank_img[:] .= RGB(0)

    cam = Camera(blank_img, [160,90], -100, [0,0,100])
    #rays = [Ray([0.0,0.0,0.0],[0.0,0.0,0.0],0) for i = 1:length(blank_img)]
    rays = Rays(AT((zeros(size(blank_img)...,3))),
                AT((zeros(size(blank_img)...,3))),
                AT((zeros(size(blank_img)...,3))))

    last_frame = 1

    for i = 1:last_frame
        angle = 2*pi*(i-1)/last_frame
        pos = [[sin(angle)*30,cos(angle)*30,-10],
               [sin(angle+(2*pi/3))*30,cos(angle+(2*pi/3))*30,-10],
               [sin(angle+(4*pi/3))*30,cos(angle+(4*pi/3))*30,-10]]
        spheres = [ColoredSphere(pos[1], 15, RGB(1, 0.25, 0.25)),
                   ColoredSphere(pos[2], 15, RGB(0.25, 1, 0.25)),
                   ColoredSphere(pos[3], 15, RGB(0.25, 0.25, 1))]


        objects = vcat(sky, spheres)
        #objects = vcat(sky, spheres, bg_spheres)

        ray_trace(objects, cam, rays; num_intersections=10,
                  filename="check"*lpad(i-1,5,"0")*".png")
    end

end
