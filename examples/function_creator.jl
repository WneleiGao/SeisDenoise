"""
   return two functions
"""
function hankel_multiplication_creator(d::Array{Tv}) where {Tv <: Union{Float32, Float64, Complex{Float32}, Complex{Float64}}}

    # dimensions of array
    dims = size(d)
    N    = length(dims)

    # Fourier transform of d
    d_hat = fft(d)

    # compute L, K
    L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
    L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
    L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1
    L4 = floor(Int64,dims[4]/2)+1; K4 = dims[4]-L4+1

    # create forward operator
    op_forward = function (x::Array{Tv})

                     # dimension of x
                     M = ndims(x)

                     # check the size of input
                     K4*K3*K2*K1 == length(x) || throw(DimensionMismatch("check length of input"))

                     # reshape x to a multi-dimensonal array, reshape doesn't change the input
                     x = reshape(v, K1, K2, K3, K4)

                     # reverse x and padding zeros, new memory is allocated
                     x_hat = reverse_order(x; n1=dims[1], n2=dims[2], n3=dims[3], n4=dims[4])
                     fft!(x_hat)

                     # element-wise multiplication
                     x_hat .= x_hat .* d_hat
                     ifft!(x_hat)

                     # make the output have same data format
                     if Tv <: AbstractFloat

                        # make the output have same dimension
                        if M == 1
                           return vec(real(x_hat[K1:dims[1],K2:dims[2],K3:dims[3],K4:dims[4]]))
                        else
                           return real(x_hat[K1:dims[1],K2:dims[2],K3:dims[3],K4:dims[4]])
                        end

                     else

                        # make the output have same dimension
                        if M == 1
                           return vec(x_hat[K1:dims[1],K2:dims[2],K3:dims[3],K4:dims[4]])
                        else
                           return x_hat[K1:dims[1],K2:dims[2],K3:dims[3],K4:dims[4]]
                        end
                     end
                 end

    # create adjoint operator
    op_adjoint = function (x::Array{Tv})

                     # dimension of x
                     M = ndims(x)

                     # check the size of input
                     L4*L3*L2*L1 == length(x) || throw(DimensionMismatch("check length of input"))

                     # reshape x to a multi-dimensonal array, reshape doesn't change the input
                     x = reshape(v, L1, L2, L3, L4)

                     # reverse v and padding zeros
                     x_hat = reverse_order(x; n1=dims[1], n2=dims[2], n3=dims[3], n4=dims[4])
                     fft!(x_hat)

                     # conjugate property of Fourier transform
                     if Tv <: Complex
                        d_hat = conj(d)
                        fft!(d_hat)
                     end

                     # element-wise multiplication
                     x_hat .= x_hat .* d_hat
                     ifft!(x_hat)

                     # make the output have same data format
                     if Tv <: AbstractFloat

                        # make the output have same dimension
                        if M == 1
                           return vec(real(x_hat[L1:dims[1],L2:dims[2],L3:dims[3],L4:dims[4]]))
                        else
                           return real(x_hat[L1:dims[1],L2:dims[2],L3:dims[3],L4:dims[4]])
                        end

                     else

                        # make the output have same dimension
                        if M == 1
                           return vec(x_hat[L1:dims[1],L2:dims[2],L3:dims[3],L4:dims[4]])
                        else
                           return x_hat[L1:dims[1],L2:dims[2],L3:dims[3],L4:dims[4]]
                        end
                     end
                 end

    return op_forward, op_adjoint
end
