# For now, all cameras are aligned on the z axis
# NOTE: think about camera dimensionality, so create a separate 2/3d camera?
#       maybe make Camera an abstract type?
struct Camera
    # Set of all pixels, counts as scene resolution
    # NOTE: Think about having a separate array for each color, RGB (or LAB)
    pixels

    # physical size of aperture
    size::Vector{Float64}

    # camera's distance from screen
    focal_length::Float64

    # camera's position
    # NOTE: add rotation (quaternions?)
    p::Vector{Float64}
end
