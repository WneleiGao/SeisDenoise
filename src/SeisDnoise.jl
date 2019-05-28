module SeisDenoise

    # the dependency of this module
    using LinearAlgebra,
          Printf,           # formated print
          DSP,              # smoothing and convolution
          FFTW              # fourier transform

    include("synthetic/synthetic.jl")          # read and write segy data, internally defined regular sampled data (borrow from rsf)
    # include("tensor/tensor.jl")  # finite-difference method for acoustic wave equation
    # include("fxy/fxy.jl")
    # include("mssa/mssa.jl")
    # include("patch/patch.jl")

end # module
