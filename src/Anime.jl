module Anime

using LinearAlgebra
using Test
using FLoops
using Images

using KernelAbstractions
using CUDA

if has_cuda_gpu()
    using CUDAKernels
end

#include("Objects.jl")
#include("Encoders.jl")

include("Camera.jl")
include("Renderers/Raytracing.jl")
include("Scenes.jl")

end # module
