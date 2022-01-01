module Anime

using LinearAlgebra
using Test
using FLoops
using Images

#include("Objects.jl")
#include("Encoders.jl")

include("Camera.jl")
include("Renderers/Raytracing.jl")
include("Scenes.jl")

end # module
