"""
   build a multi-level circulant matrix up to five dimension
"""
function build_circulant_matrix(d::Array{Tv}) where {Tv<:Number}

    # determine the dimensions of an array
    dims = size(d)
    N    = length(dims)

    # level 1
    if N == 1

       # allocate memory for Circulant matrix
       C = Array{Tv,2}(undef, dims[1], dims[1])

       for j1 = 1 : dims[1]
           for i1 = 1 : dims[1]
               n1 = i1 - j1 + 1
               n1 = n1 >= 1 ? n1 : n1+dims[1]

               C[i1,j1] = d[n1]
           end
       end


    elseif N == 2

       # allocate memory for Circulant matrix
       C = Array{Tv,2}(undef, prod(dims[1:2]), prod(dims[1:2]))

       # build circulant
       for j2 = 1 : dims[2]
           c2 = (j2-1)*dims[1]
           for i2 = 1 : dims[2]
               r2 = (i2-1)*dims[1]
               n2 = i2 - j2 + 1
               n2 = n2 >= 1 ? n2 : n2+dims[2]

               for j1 = 1 : dims[1]
                   c1 = c2+j1
                   for i1 = 1 : dims[1]
                       r1 = r2+i1
                       n1 = i1 - j1 + 1
                       n1 = n1 >= 1 ? n1 : n1+dims[1]

                       C[r1,c1] = d[n1,n2]
                   end
               end
           end
       end


    elseif N == 3

       # allocate memory for Circulant matrix
       C = Array{Tv,2}(undef, prod(dims[1:3]), prod(dims[1:3]))

       # build circulant
       for j3 = 1 : dims[3]
           c3 = (j3-1)*dims[2]*dims[1]
           for i3 = 1 : dims[3]
               r3 = (i3-1)*dims[2]*dims[1]
               n3 = i3 - j3 + 1
               n3 = n3 >= 1 ? n3 : n3+dims[3]

               for j2 = 1 : dims[2]
                   c2 = (j2-1)*dims[1]
                   for i2 = 1 : dims[2]
                       r2 = (i2-1)*dims[1]
                       n2 = i2 - j2 + 1
                       n2 = n2 >= 1 ? n2 : n2+dims[2]

                       for j1 = 1 : dims[1]
                           c1 = c3+c2+j1
                           for i1 = 1 : dims[1]
                               r1 = r3+r2+i1
                               n1 = i1 - j1 + 1
                               n1 = n1 >= 1 ? n1 : n1+dims[1]

                               C[r1,c1] = d[n1,n2,n3]
                           end
                       end
                   end
               end
           end
       end


    elseif N == 4

       # allocate memory for Circulant matrix
       C = Array{Tv,2}(undef, prod(dims[1:4]), prod(dims[1:4]))

       # build circulant
       for j4 = 1 : dims[4]
           c4 = (j4-1)*dims[3]*dims[2]*dims[1]
           for i4 = 1 : dims[4]
               r4 = (i4-1)*dims[3]*dims[2]*dims[1]
               n4 = i4 - j4 + 1
               n4 = n4 >= 1 ? n4 : n4+dims[4]

               for j3 = 1 : dims[3]
                   c3 = (j3-1)*dims[2]*dims[1]
                   for i3 = 1 : dims[3]
                       r3 = (i3-1)*dims[2]*dims[1]
                       n3 = i3 - j3 + 1
                       n3 = n3 >= 1 ? n3 : n3+dims[3]

                       for j2 = 1 : dims[2]
                           c2 = (j2-1)*dims[1]
                           for i2 = 1 : dims[2]
                               r2 = (i2-1)*dims[1]
                               n2 = i2 - j2 + 1
                               n2 = n2 >= 1 ? n2 : n2+dims[2]

                               for j1 = 1 : dims[1]
                                   c1 = c4+c3+c2+j1
                                   for i1 = 1 : dims[1]
                                       r1 = r4+r3+r2+i1
                                       n1 = i1 - j1 + 1
                                       n1 = n1 >= 1 ? n1 : n1+dims[1]

                                       C[r1,c1] = d[n1,n2,n3,n4]
                                   end
                               end
                           end
                       end
                   end
               end
           end
       end


    elseif N == 5

       # allocate memory for Circulant matrix
       C = Array{Tv,2}(undef, prod(dims[1:5]), prod(dims[1:5]))

       # build circulant
       for j5 = 1 : dims[5]
           c5 = (j5-1)*dims[4]*dims[3]*dims[2]*dims[1]
           for i5 = 1 : dims[5]
               r5 = (i5-1)*dims[4]*dims[3]*dims[2]*dims[1]
               n5 = i5 - j5 + 1
               n5 = n5 >= 1 ? n5 : n5+dims[5]

               for j4 = 1 : dims[4]
                   c4 = (j4-1)*dims[3]*dims[2]*dims[1]
                   for i4 = 1 : dims[4]
                       r4 = (i4-1)*dims[3]*dims[2]*dims[1]
                       n4 = i4 - j4 + 1
                       n4 = n4 >= 1 ? n4 : n4+dims[4]

                       for j3 = 1 : dims[3]
                           c3 = (j3-1)*dims[2]*dims[1]
                           for i3 = 1 : dims[3]
                               r3 = (i3-1)*dims[2]*dims[1]
                               n3 = i3 - j3 + 1
                               n3 = n3 >= 1 ? n3 : n3+dims[3]

                               for j2 = 1 : dims[2]
                                   c2 = (j2-1)*dims[1]
                                   for i2 = 1 : dims[2]
                                       r2 = (i2-1)*dims[1]
                                       n2 = i2 - j2 + 1
                                       n2 = n2 >= 1 ? n2 : n2+dims[2]

                                       for j1 = 1 : dims[1]
                                           c1 = c5+c4+c3+c2+j1
                                           for i1 = 1 : dims[1]
                                               r1 = r5+r4+r3+r2+i1
                                               n1 = i1 - j1 + 1
                                               n1 = n1 >= 1 ? n1 : n1+dims[1]

                                               C[r1,c1] = d[n1,n2,n3,n4,n5]
                                           end
                                       end
                                   end
                               end
                           end
                       end
                   end
               end
           end
       end


    else
       error("only surpport up to five dimension")
    end

    return C
end

"""
   build a multi-level Toeplitze matrix from a multi-dimensional array
"""
function build_toeplitz_matrix(d::Array{Tv}) where {Tv <: Number}

    # determine the dimensions of an array
    dims = size(d)
    N    = length(dims)

    # level 1
    if N == 1
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1

       # allocate memory for Hankel matrix
       T  = Array{Tv,2}(undef, L1, K1)

       # build hankel
       for j1 = 1 : K1
           for i1 = 1 : L1
               n1 = K1 + i1 - j1

               T[i1,j1] = d[n1]
           end
       end


    # level 2
    elseif N == 2
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1

       # allocate memory for Hankel matrix
       T  = Array{Tv,2}(undef, L2*L1, K2*K1)

       # build Hankel
       for j2 = 1 : K2
           c2 = (j2-1)*K1
           for i2 = 1 : L2
               r2 = (i2-1)*L1
               n2 = K2 + i2 - j2

               # first layer
               for j1 = 1 : K1
                   c1 = c2 + j1
                   for i1 = 1 : L1
                       r1 = r2 + i1
                       n1 = K1 + i1 - j1

                       T[r1, c1] = d[n1,n2]
                   end
               end #first
           end
       end #second


    # level 3
    elseif N == 3
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
       L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1

       # allocate memory for Hankel matrix
       T  = Array{Tv,2}(undef, L3*L2*L1, K3*K2*K1)

       # build Hankel
       for j3 = 1 : K3
           c3 = (j3-1)*K2*K1
           for i3 = 1 : L3
               r3 = (i3-1)*L2*L1
               n3 = K3 + i3 - j3

               # second layer
               for j2 = 1 : K2
                   c2 = (j2-1)*K1
                   for i2 = 1 : L2
                       r2 = (i2-1)*L1
                       n2 = K2 + i2 - j2

                       # first layer
                       for j1 = 1 : K1
                           c1 = c3 + c2 + j1
                           for i1 = 1 : L1
                               r1 = r3 + r2 + i1
                               n1 = K1 + i1 - j1

                               T[r1, c1] = d[n1,n2,n3]
                           end
                       end #first
                   end
               end #second
           end
       end #third


    # level 4
    elseif N == 4
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
       L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1
       L4 = floor(Int64,dims[4]/2)+1; K4 = dims[4]-L4+1

       # allocate memory for Hankel matrix
       T  = Array{Tv,2}(undef, L4*L3*L2*L1, K4*K3*K2*K1)

       # build Hankel
       for j4 = 1 : K4
           c4 = (j4-1)*K3*K2*K1
           for i4 = 1 : L4
               r4 = (i4-1)*L3*L2*L1
               n4 = K4 + i4 - j4

               # third layer
               for j3 = 1 : K3
                   c3 = (j3-1)*K2*K1
                   for i3 = 1 : L3
                       r3 = (i3-1)*L2*L1
                       n3 = K3 + i3 - j3

                       # second layer
                       for j2 = 1 : K2
                           c2 = (j2-1)*K1
                           for i2 = 1 : L2
                               r2 = (i2-1)*L1
                               n2 = K2 + i2 - j2

                               # first layer
                               for j1 = 1 : K1
                                   c1 = c4 + c3 + c2 + j1
                                   for i1 = 1 : L1
                                       r1 = r4 + r3 + r2 + i1
                                       n1 = K1 + i1 - j1

                                       T[r1, c1] = d[n1,n2,n3,n4]
                                   end
                               end #first
                           end
                       end #second
                   end
               end #third
           end
       end #fourth


    # level 5
    elseif N == 5
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
       L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1
       L4 = floor(Int64,dims[4]/2)+1; K4 = dims[4]-L4+1
       L5 = floor(Int64,dims[5]/2)+1; K5 = dims[5]-L5+1

       # allocate memory for Hankel matrix
       T  = Array{Tv,2}(undef, L5*L4*L3*L2*L1, K5*K4*K3*K2*K1)

       # build Hankel
       for j5 = 1 : K5
           c5 = (j5-1)*K4*K3*K2*K1
           for i5 = 1 : L5
               r5 = (i5-1)*L4*L3*L2*L1
               n5 = K5 + i5 - j5

               for j4 = 1 : K4
                   c4 = (j4-1)*K3*K2*K1
                   for i4 = 1 : L4
                       r4 = (i4-1)*L3*L2*L1
                       n4 = K4 + i4 - j4

                       # third layer
                       for j3 = 1 : K3
                           c3 = (j3-1)*K2*K1
                           for i3 = 1 : L3
                               r3 = (i3-1)*L2*L1
                               n3 = K3 + i3 - j3

                               # second layer
                               for j2 = 1 : K2
                                   c2 = (j2-1)*K1
                                   for i2 = 1 : L2
                                       r2 = (i2-1)*L1
                                       n2 = K2 + i2 - j2

                                       # first layer
                                       for j1 = 1 : K1
                                           c1 = c5 + c4 + c3 + c2 + j1
                                           for i1 = 1 : L1
                                               r1 = r5 + r4 + r3 + r2 + i1
                                               n1 = K1 + i1 - j1

                                               T[r1, c1] = d[n1,n2,n3,n4,n5]
                                           end
                                       end #first
                                   end
                               end #second
                           end
                       end #third
                   end
               end #fourth
           end
       end # fiveth


    else
       error("only support up to 5 dimensions")
    end

    return T
end

"""
   efficient way to compute the product between toeplitze matrix and multi-dimensional array
"""
function toeplitz_multiplication(d::Array{Tv}, v::Array{Tv}) where {Tv <: Number}

    # dimensions of array
    dims = size(d)
    N    = length(dims)

    # Fourier transform of d
    d_hat = fft(d)

    if N == 1
       # compute L, K
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1

       if flag == "forward"

          # check the size of v
          size(v, 1) == K1 || error(DimensionMismatch("check the size of v"))

          # padding zeros
          if Tv <: AbstractFloat
             v_hat = zeros(Complex{Tv}, dims[1])
          else
             v_hat = zeros(Tv, dims[1])
          end
          v_hat[1:K1] .= v

          # fourier transform
          fft!(v_hat)

          # element-wise multiplication
          r = ifft(d_hat .* v_hat)

          # return the last L1 element
          return r[K1:dims[1]]

       elseif flag == "adjoint"

          # check the size of v
          size(v, 1) == L1 || error(DimensionMismatch("check the size of v"))

          # padding zeros
          if Tv <: AbstractFloat
             v_hat = zeros(Complex{Tv}, dims[1])
          else
             v_hat = zeros(Tv, dims[1])
          end
          v_hat[1:L1] .= v

          # fourier transform
          fft!(v_hat)

          # conjugate property
          d_tilde = copy(d_hat)
          for i1 = 1 : dims[1]
              j1 = i1 == 1 ? 1 : dims[1]-i1+2

              d_tilde[i1] = conj(d_hat[j1])
          end

          # element-wise multiplication
          r = ifft(d_tilde .* v_hat)

          # return the last K1 element
          return r[L1:dims[1]]

       else
          error("non-surpported operation")
       end


    elseif N == 2
       # compute L, K
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1

       if flag == "forward"

          # check the size of v
          size(v, 1) == K1 || error(DimensionMismatch("check the size of v"))
          size(v, 2) == K2 || error(DimensionMismatch("check the size of v"))

          # padding zeros
          if Tv <: AbstractFloat
             v_hat = zeros(Complex{Tv}, dims[1], dims[2])
          else
             v_hat = zeros(Tv, dims[1], dims[2])
          end
          v_hat[1:L1,1:L2] .= v

          # fourier transform
          fft!(v_hat)

          # return the last L1 element
          return r[K1:dims[1],K2:dims[2]]

       elseif flag == "adjoint"

          # check the size of v
          size(v, 1) == L1 || error(DimensionMismatch("check the size of v"))
          size(v, 2) == L2 || error(DimensionMismatch("check the size of v"))

          # reverse v and padding zeros
          v_hat = reverse_order(v; n1=dims[1], n2=dims[2])
          v_hat = fft(v_hat)

          # conjugate property
          d_tilde = copy(d_hat)
          for i2 = 1 : dims[2]
              j2 = i2 == 1 ? 1 : dims[2]-i2+2

              for i1 = 1 : dims[1]
                  j1 = i1 == 1 ? 1 : dims[1]-i1+2

                  d_tilde[i1,i2] = conj(d_hat[j1,j2])
              end
          end

          # element-wise multiplication
          r = ifft(d_tilde .* v_hat)

          # return the last K1 element
          return r[L1:dims[1], L2:dims[2]]

       else
          error("non-surpported operation")
       end


    elseif N == 3
       # compute L, K
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
       L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1

       if flag == "forward"

          # check the size of v
          size(v, 1) == K1 || error(DimensionMismatch("check the size of v"))
          size(v, 2) == K2 || error(DimensionMismatch("check the size of v"))
          size(v, 3) == K3 || error(DimensionMismatch("check the size of v"))

          # reverse v and padding zeros
          v_hat = reverse_order(v; n1=dims[1], n2=dims[2], n3=dims[3])
          v_hat = fft(v_hat)

          # element-wise multiplication
          r = ifft(d_hat .* v_hat)

          # return the last L1 element
          return r[K1:dims[1],K2:dims[2],K3:dims[3]]

       elseif flag == "adjoint"

          # check the size of v
          size(v, 1) == L1 || error(DimensionMismatch("check the size of v"))
          size(v, 2) == L2 || error(DimensionMismatch("check the size of v"))
          size(v, 3) == L3 || error(DimensionMismatch("check the size of v"))

          # reverse v and padding zeros
          v_hat = reverse_order(v; n1=dims[1], n2=dims[2], n3=dims[3])
          v_hat = fft(v_hat)

          # conjugate property
          d_tilde = copy(d_hat)
          for i3 = 1 : dims[3]
              j3 = i3 == 1 ? 1 : dims[3]-i3+2

              for i2 = 1 : dims[2]
                  j2 = i2 == 1 ? 1 : dims[2]-i2+2

                  for i1 = 1 : dims[1]
                      j1 = i1 == 1 ? 1 : dims[1]-i1+2

                      d_tilde[i1,i2,i3] = conj(d_hat[j1,j2,j3])
                  end
              end
          end

          # element-wise multiplication
          r = ifft(d_tilde .* v_hat)

          # return the last K1 element
          return r[L1:dims[1],L2:dims[2],L3:dims[3]]

       else
          error("non-surpported operation")
       end


    elseif N == 4
       # compute L, K
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
       L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1
       L4 = floor(Int64,dims[4]/2)+1; K4 = dims[4]-L4+1

       if flag == "forward"

          # check the size of v
          size(v, 1) == K1 || error(DimensionMismatch("check the size of v"))
          size(v, 2) == K2 || error(DimensionMismatch("check the size of v"))
          size(v, 3) == K3 || error(DimensionMismatch("check the size of v"))
          size(v, 4) == K4 || error(DimensionMismatch("check the size of v"))

          # reverse v and padding zeros
          v_hat = reverse_order(v; n1=dims[1], n2=dims[2], n3=dims[3], n4=dims[4])
          v_hat = fft(v_hat)

          # element-wise multiplication
          r = ifft(d_hat .* v_hat)

          # return the last L1 element
          return r[K1:dims[1],K2:dims[2],K3:dims[3],K4:dims[4]]

       elseif flag == "adjoint"

          # check the size of v
          size(v, 1) == L1 || error(DimensionMismatch("check the size of v"))
          size(v, 2) == L2 || error(DimensionMismatch("check the size of v"))
          size(v, 3) == L3 || error(DimensionMismatch("check the size of v"))
          size(v, 4) == L4 || error(DimensionMismatch("check the size of v"))

          # reverse v and padding zeros
          v_hat = reverse_order(v; n1=dims[1], n2=dims[2], n3=dims[3], n4=dims[4])
          v_hat = fft(v_hat)

          # conjugate property
          d_tilde = copy(d_hat)
          for i4 = 1 : dims[4]
              j4 = i4 == 1 ? 1 : dims[4]-i4+2

              for i3 = 1 : dims[3]
                  j3 = i3 == 1 ? 1 : dims[3]-i3+2

                  for i2 = 1 : dims[2]
                      j2 = i2 == 1 ? 1 : dims[2]-i2+2

                      for i1 = 1 : dims[1]
                          j1 = i1 == 1 ? 1 : dims[1]-i1+2

                          d_tilde[i1,i2,i3,i4] = conj(d_hat[j1,j2,j3,j4])
                      end
                  end
              end
          end

          # element-wise multiplication
          r = ifft(d_tilde .* v_hat)

          # return the last K1 element
          return r[L1:dims[1],L2:dims[2],L3:dims[3],L4:dims[4]]

       else
          error("non-surpported operation")
       end


    elseif N == 5
       # compute L, K
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
       L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1
       L4 = floor(Int64,dims[4]/2)+1; K4 = dims[4]-L4+1
       L5 = floor(Int64,dims[5]/2)+1; K5 = dims[5]-L5+1

       if flag == "forward"

          # check the size of v
          size(v, 1) == K1 || error(DimensionMismatch("check the size of v"))
          size(v, 2) == K2 || error(DimensionMismatch("check the size of v"))
          size(v, 3) == K3 || error(DimensionMismatch("check the size of v"))
          size(v, 4) == K4 || error(DimensionMismatch("check the size of v"))
          size(v, 5) == K5 || error(DimensionMismatch("check the size of v"))

          # reverse v and padding zeros
          v_hat = reverse_order(v; n1=dims[1], n2=dims[2], n3=dims[3], n4=dims[4], n5=dims[5])
          v_hat = fft(v_hat)

          # element-wise multiplication
          r = ifft(d_hat .* v_hat)

          # return the last L1 element
          return r[K1:dims[1],K2:dims[2],K3:dims[3],K4:dims[4],K5:dims[5]]

       elseif flag == "adjoint"

          # check the size of v
          size(v, 1) == L1 || error(DimensionMismatch("check the size of v"))
          size(v, 2) == L2 || error(DimensionMismatch("check the size of v"))
          size(v, 3) == L3 || error(DimensionMismatch("check the size of v"))
          size(v, 4) == L4 || error(DimensionMismatch("check the size of v"))
          size(v, 5) == L5 || error(DimensionMismatch("check the size of v"))

          # reverse v and padding zeros
          v_hat = reverse_order(v; n1=dims[1], n2=dims[2], n3=dims[3], n4=dims[4], n5=dims[5])
          v_hat = fft(v_hat)

          # conjugate property
          d_tilde = copy(d_hat)
          for i5 = 1 : dims[5]
              j5 = i5 == 1 ? 1 : dims[5]-i5+2

              for i4 = 1 : dims[4]
                  j4 = i4 == 1 ? 1 : dims[4]-i4+2

                  for i3 = 1 : dims[3]
                      j3 = i3 == 1 ? 1 : dims[3]-i3+2

                      for i2 = 1 : dims[2]
                          j2 = i2 == 1 ? 1 : dims[2]-i2+2

                          for i1 = 1 : dims[1]
                              j1 = i1 == 1 ? 1 : dims[1]-i1+2

                              d_tilde[i1,i2,i3,i4,i5] = conj(d_hat[j1,j2,j3,j4,j5])
                          end
                      end
                  end
              end
          end

          # element-wise multiplication
          r = ifft(d_tilde .* v_hat)

          # return the last K1 element
          return r[L1:dims[1],L2:dims[2],L3:dims[3],L4:dims[4],L5:dims[5]]

      else
         error("only supported up to five dimension")
      end # forward or adjoint

   end # dimension case


end

"""
   build a multi-level block Hankel matrix from a multi-dimensional array
"""
function build_hankel_matrix(d::Array{Tv}) where {Tv <: Number}

    # determine the dimensions of an array
    dims = size(d)
    N    = length(dims)

    # level 1
    if N == 1
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1

       # allocate memory for Hankel matrix
       H  = Array{Tv,2}(undef, L1, K1)

       # build hankel
       for j1 = 1 : K1
           for i1 = 1 : L1
               n1 = i1+j1-1
               H[i1,j1] = d[n1]
           end
       end


    # level 2
    elseif N == 2
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1

       # allocate memory for Hankel matrix
       H  = Array{Tv,2}(undef, L2*L1, K2*K1)

       # build Hankel
       for j2 = 1 : K2
           c2 = (j2-1)*K1
           for i2 = 1 : L2
               r2 = (i2-1)*L1
               n2 = i2 + j2 - 1

               # first layer
               for j1 = 1 : K1
                   c1 = c2 + j1
                   for i1 = 1 : L1
                       r1 = r2 + i1
                       n1 = i1 + j1 - 1
                       H[r1, c1] = d[n1,n2]
                   end
               end #first
           end
       end #second


    # level 3
    elseif N == 3
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
       L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1

       # allocate memory for Hankel matrix
       H  = Array{Tv,2}(undef, L3*L2*L1, K3*K2*K1)

       # build Hankel
       for j3 = 1 : K3
           c3 = (j3-1)*K2*K1
           for i3 = 1 : L3
               r3 = (i3-1)*L2*L1
               n3 = i3 + j3 - 1

               # second layer
               for j2 = 1 : K2
                   c2 = (j2-1)*K1
                   for i2 = 1 : L2
                       r2 = (i2-1)*L1
                       n2 = i2 + j2 - 1

                       # first layer
                       for j1 = 1 : K1
                           c1 = c3 + c2 + j1
                           for i1 = 1 : L1
                               r1 = r3 + r2 + i1
                               n1 = i1 + j1 - 1
                               H[r1, c1] = d[n1,n2,n3]
                           end
                       end #first
                   end
               end #second
           end
       end #third


    # level 4
    elseif N == 4
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
       L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1
       L4 = floor(Int64,dims[4]/2)+1; K4 = dims[4]-L4+1

       # allocate memory for Hankel matrix
       H  = Array{Tv,2}(undef, L4*L3*L2*L1, K4*K3*K2*K1)

       # build Hankel
       for j4 = 1 : K4
           c4 = (j4-1)*K3*K2*K1
           for i4 = 1 : L4
               r4 = (i4-1)*L3*L2*L1
               n4 = i4 + j4 - 1

               # third layer
               for j3 = 1 : K3
                   c3 = (j3-1)*K2*K1
                   for i3 = 1 : L3
                       r3 = (i3-1)*L2*L1
                       n3 = i3 + j3 - 1

                       # second layer
                       for j2 = 1 : K2
                           c2 = (j2-1)*K1
                           for i2 = 1 : L2
                               r2 = (i2-1)*L1
                               n2 = i2 + j2 - 1

                               # first layer
                               for j1 = 1 : K1
                                   c1 = c4 + c3 + c2 + j1
                                   for i1 = 1 : L1
                                       r1 = r4 + r3 + r2 + i1
                                       n1 = i1 + j1 - 1
                                       H[r1, c1] = d[n1,n2,n3,n4]
                                   end
                               end #first
                           end
                       end #second
                   end
               end #third
           end
       end #fourth


    # level 5
    elseif N == 5
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
       L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1
       L4 = floor(Int64,dims[4]/2)+1; K4 = dims[4]-L4+1
       L5 = floor(Int64,dims[5]/2)+1; K5 = dims[5]-L5+1

       # allocate memory for Hankel matrix
       H  = Array{Tv,2}(undef, L5*L4*L3*L2*L1, K5*K4*K3*K2*K1)

       # build Hankel
       for j5 = 1 : K5
           c5 = (j5-1)*K4*K3*K2*K1
           for i5 = 1 : L5
               r5 = (i5-1)*L4*L3*L2*L1
               n5 = i5 + j5 - 1

               for j4 = 1 : K4
                   c4 = (j4-1)*K3*K2*K1
                   for i4 = 1 : L4
                       r4 = (i4-1)*L3*L2*L1
                       n4 = i4 + j4 - 1

                       # third layer
                       for j3 = 1 : K3
                           c3 = (j3-1)*K2*K1
                           for i3 = 1 : L3
                               r3 = (i3-1)*L2*L1
                               n3 = i3 + j3 - 1

                               # second layer
                               for j2 = 1 : K2
                                   c2 = (j2-1)*K1
                                   for i2 = 1 : L2
                                       r2 = (i2-1)*L1
                                       n2 = i2 + j2 - 1

                                       # first layer
                                       for j1 = 1 : K1
                                           c1 = c5 + c4 + c3 + c2 + j1
                                           for i1 = 1 : L1
                                               r1 = r5 + r4 + r3 + r2 + i1
                                               n1 = i1 + j1 - 1
                                               H[r1, c1] = d[n1,n2,n3,n4,n5]
                                           end
                                       end #first
                                   end
                               end #second
                           end
                       end #third
                   end
               end #fourth
           end
       end # fiveth


    else
       error("only support up to 5 dimensions")
    end

    return H
end

"""
   count the number of times of the elements of a multi-dimensional array get repeated
when building a hankel matrix
"""
function count_hankel(dims::Union{Ti,Vector{Ti}}) where {Ti<:Int64}

    # convert to a vector when the input is a scalar
    if typeof(dims) == Int64
       dims = [dims]
    end

    # determine the dimensions of an array
    N = length(dims)

    # level 1
    if N == 1

       # allocate memory to record the number
       count_num = zeros(Int64, dims[1])

       # size of hankel matrix
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1

       # assign value
       for i1 = 1 : dims[1]
           if i1 <= K1
              n1 = i1
           end
           if L1 > K1 && i1 == L1
              n1 = L1-1
           end
           if i1 > L1
              n1 = dims[1]-i1+1
           end

           count_num[i1] = n1
       end

    # level 2
    elseif N == 2

       # allocate memory to record the number
       count_num = zeros(Int64, dims[1], dims[2])

       # size of hankel matrix
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1

       # second dimension
       for i2 = 1 : dims[2]
           if i2 <= K2
              n2 = i2
           end
           if L2 > K2 && i2 == L2
              n2 = L2-1
           end
           if i2 > L2
              n2 = dims[2]-i2+1
           end

           # first dimension
           for i1 = 1 : dims[1]
               if i1 <= K1
                  n1 = i1
               end
               if L1 > K1 && i1 == L1
                  n1 = L1-1
               end
               if i1 > L1
                  n1 = dims[1]-i1+1
               end

               count_num[i1,i2] = n2 * n1
           end
       end


    # level 3
    elseif N == 3

       # allocate memory to record the number
       count_num = zeros(Int64, dims[1], dims[2], dims[3])

       # size of hankel matrix
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
       L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1

       # third dimension
       for i3 = 1 : dims[3]
           if i3 <= K3
              n3 = i3
           end
           if L3 > K3 && i3 == L3
              n3 = L3-1
           end
           if i3 > L3
              n3 = dims[3]-i3+1
           end

           # second dimension
           for i2 = 1 : dims[2]
               if i2 <= K2
                  n2 = i2
               end
               if L2 > K2 && i2 == L2
                  n2 = L2-1
               end
               if i2 > L2
                  n2 = dims[2]-i2+1
               end

               # first dimension
               for i1 = 1 : dims[1]
                   if i1 <= K1
                      n1 = i1
                   end
                   if L1 > K1 && i1 == L1
                      n1 = L1-1
                   end
                   if i1 > L1
                      n1 = dims[1]-i1+1
                   end

                   count_num[i1,i2,i3] = n3 * n2 * n1
               end
           end
       end


    # level 4
    elseif N == 4

       # allocate memory to record the number
       count_num = zeros(Int64, dims[1], dims[2], dims[3], dims[4])

       # size of hankel matrix
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
       L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1
       L4 = floor(Int64,dims[4]/2)+1; K4 = dims[4]-L4+1

       # fourth dimension
       for i4 = 1 : dims[4]
           if i4 <= K4
              n4 = i4
           end
           if L4 > K4 && i4 == L4
              n4 = L4-1
           end
           if i4 > L4
              n4 = dims[4]-i4+1
           end

           # third dimension
           for i3 = 1 : dims[3]
               if i3 <= K3
                  n3 = i3
               end
               if L3 > K3 && i3 == L3
                  n3 = L3-1
               end
               if i3 > L3
                  n3 = dims[3]-i3+1
               end

               # second dimension
               for i2 = 1 : dims[2]
                   if i2 <= K2
                      n2 = i2
                   end
                   if L2 > K2 && i2 == L2
                      n2 = L2-1
                   end
                   if i2 > L2
                      n2 = dims[2]-i2+1
                   end

                   # first dimension
                   for i1 = 1 : dims[1]
                       if i1 <= K1
                          n1 = i1
                       end
                       if L1 > K1 && i1 == L1
                          n1 = L1-1
                       end
                       if i1 > L1
                          n1 = dims[1]-i1+1
                       end

                       count_num[i1,i2,i3,i4] = n4 * n3 * n2 * n1
                   end
               end
           end
       end


    # level 5
    elseif N == 5

       # allocate memory to record the number
       count_num = zeros(Int64, dims[1], dims[2], dims[3], dims[4], dims[5])

       # size of hankel matrix
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
       L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1
       L4 = floor(Int64,dims[4]/2)+1; K4 = dims[4]-L4+1
       L5 = floor(Int64,dims[5]/2)+1; K5 = dims[5]-L5+1

       # assign value
       for i5 = 1 : dims[5]
           if i5 <= K5
              n5 = i5
           end
           if L5 > K5 && i5 == L5
              n5 = L5-1
           end
           if i5 > L5
              n5 = dims[5]-i5+1
           end

           # fourth dimension
           for i4 = 1 : dims[4]
               if i4 <= K4
                  n4 = i4
               end
               if L4 > K4 && i4 == L4
                  n4 = L4-1
               end
               if i4 > L4
                  n4 = dims[4]-i4+1
               end

               # third dimension
               for i3 = 1 : dims[3]
                   if i3 <= K3
                      n3 = i3
                   end
                   if L3 > K3 && i3 == L3
                      n3 = L3-1
                   end
                   if i3 > L3
                      n3 = dims[3]-i3+1
                   end

                   # second dimension
                   for i2 = 1 : dims[2]
                       if i2 <= K2
                          n2 = i2
                       end
                       if L2 > K2 && i2 == L2
                          n2 = L2-1
                       end
                       if i2 > L2
                          n2 = dims[2]-i2+1
                       end

                       # first dimension
                       for i1 = 1 : dims[1]
                           if i1 <= K1
                              n1 = i1
                           end
                           if L1 > K1 && i1 == L1
                              n1 = L1-1
                           end
                           if i1 > L1
                              n1 = dims[1]-i1+1
                           end

                           count_num[i1,i2,i3,i4,i5] = n5 * n4 * n3 * n2 * n1
                       end
                   end
               end
           end
       end


    else
       error("only support up to 5 dimensions")
    end

    return count_num
end

"""
   reverse the order of the element of a multi-dimensional array and padding zeros
to the resultant multi-dimensional array
"""
function reverse_order(d::Array{Tv};
                       n1=0, n2=0, n3=0, n4=0, n5=0) where {Tv<:Number}

    # get the dimensions of input array
    dims = size(d)
    N    = length(dims)

    if N == 1

       # padding zeros if needed
       n1 = n1 < dims[1] ? dims[1] : n1

       # allocate memory for the new array
       r = zeros(Tv, n1)

       # first dimension
       for i1 = 1 : dims[1]
           j1 = dims[1]-i1+1

           # assign value
           r[i1] = d[j1]
       end


    elseif N == 2

      # padding zeros if needed
      n1 = n1 < dims[1] ? dims[1] : n1
      n2 = n2 < dims[2] ? dims[2] : n2

      # allocate memory for the new array
      r = zeros(Tv, n1, n2)

      # second dimension
      for i2 = 1 : dims[2]
          j2 = dims[2]-i2+1

          # first dimension
          for i1 = 1 : dims[1]
              j1 = dims[1]-i1+1

              # assign value
              r[i1,i2] = d[j1,j2]
          end
      end


    elseif N == 3

       # padding zeros if needed
       n1 = n1 < dims[1] ? dims[1] : n1
       n2 = n2 < dims[2] ? dims[2] : n2
       n3 = n3 < dims[3] ? dims[3] : n3

       # allocate memory for the new array
       r = zeros(Tv, n1, n2, n3)

       # third dimension
       for i3 = 1 : dims[3]
           j3 = dims[3]-i3+1

           # second dimension
           for i2 = 1 : dims[2]
               j2 = dims[2]-i2+1

               # first dimension
               for i1 = 1 : dims[1]
                   j1 = dims[1]-i1+1

                   # assign value
                   r[i1,i2,i3] = d[j1,j2,j3]
               end
           end
       end


    elseif N == 4

       # padding zeros if needed
       n1 = n1 < dims[1] ? dims[1] : n1
       n2 = n2 < dims[2] ? dims[2] : n2
       n3 = n3 < dims[3] ? dims[3] : n3
       n4 = n4 < dims[4] ? dims[4] : n4

       # allocate memory for the new array
       r = zeros(Tv, n1, n2, n3, n4)

       # fourth dimension
       for i4 = 1 : dims[4]
           j4 = dims[4]-i4+1

           # third dimension
           for i3 = 1 : dims[3]
               j3 = dims[3]-i3+1

               # second dimension
               for i2 = 1 : dims[2]
                   j2 = dims[2]-i2+1

                   # first dimension
                   for i1 = 1 : dims[1]
                       j1 = dims[1]-i1+1

                       # assign value
                       r[i1,i2,i3,i4] = d[j1,j2,j3,j4]
                   end
               end
           end
       end


    elseif N == 5

       # padding zeros if needed
       n1 = n1 < dims[1] ? dims[1] : n1
       n2 = n2 < dims[2] ? dims[2] : n2
       n3 = n3 < dims[3] ? dims[3] : n3
       n4 = n4 < dims[4] ? dims[4] : n4
       n5 = n5 < dims[5] ? dims[5] : n5

       # allocate memory for the new array
       r = zeros(Tv, n1, n2, n3, n4, n5)

       # fiveth dimension
       for i5 = 1 : dims[5]
           j5 = dims[5]-i5+1

           # fourth dimension
           for i4 = 1 : dims[4]
               j4 = dims[4]-i4+1

               # third dimension
               for i3 = 1 : dims[3]
                   j3 = dims[3]-i3+1

                   # second dimension
                   for i2 = 1 : dims[2]
                       j2 = dims[2]-i2+1

                       # first dimension
                       for i1 = 1 : dims[1]
                           j1 = dims[1]-i1+1

                           # assign value
                           r[i1,i2,i3,i4,i5] = d[j1,j2,j3,j4,j5]
                       end
                   end
               end
           end
       end


    else
       error("only support up to 5 dimensions")
    end

    return r
end


# compute toeplitz times a vector
function toeplitz_times_vector(c::Vector{Tv}, v::Vector{Tv}; embed_flag=1) where {Tv<:Number}

    # length of vector
    N = length(c)

    # compute L, K
    L = floor(Int64, N/2) + 1
    K = N + 1 - L

    # check the length of v
    length(v) == K || DimensionMismatch()

    # padding zeros to v
    v_pad = vcat(v, zeros(Tv,L-1))

    # fast computation
    if embed_flag == 1

       r_pad = real(ifft( fft(c) .* fft(v_pad) ) )

       # return the last L element N-L+1 : N and K = N-L+1
       return r_pad[K:N]

    elseif embed_flag == 2

       c_hat = Vector{Tv}(undef, N)
       # the second way to compute
       c_hat[1:L]   .= c[K:N]

       # move the first K-l element to the end
       c_hat[L+1:N] .= c[1:K-1]

       # fast computation
       r_pad = real(ifft( fft(c_hat) .* fft(v_pad) ) )

       # return the first L element
       return r_pad[1:L]
    end
end


# Hankel matrix or its Hermitian transpose times a vector
function hankel_multiplication(d::Array{Tv}, v::Array{Tv}; flag="forward") where {Tv<:Number}

    # dimensions of array
    dims = size(d)
    N    = length(dims)

    # Fourier transform of d
    d_hat = fft(d)

    if N == 1
       # compute L, K
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1

       if flag == "forward"

          # check the size of v
          size(v, 1) == K1 || error(DimensionMismatch("check the size of v"))

          # reverse v and padding zeros
          v_hat = reverse_order(v; n1=dims[1])
          v_hat = fft(v_hat)

          # element-wise multiplication
          r = ifft(d_hat .* v_hat)

          # return the last L1 element
          return r[K1:dims[1]]

       elseif flag == "adjoint"

          # check the size of v
          size(v, 1) == L1 || error(DimensionMismatch("check the size of v"))

          # reverse v and padding zeros
          v_hat = reverse_order(v; n1=dims[1])
          v_hat = fft(v_hat)

          # conjugate property
          d_tilde = copy(d_hat)
          for i1 = 1 : dims[1]
              j1 = i1 == 1 ? 1 : dims[1]-i1+2

              d_tilde[i1] = conj(d_hat[j1])
          end

          # element-wise multiplication
          r = ifft(d_tilde .* v_hat)

          # return the last K1 element
          return r[L1:dims[1]]

       else
          error("non-surpported operation")
       end


    elseif N == 2
       # compute L, K
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1

       if flag == "forward"

          # check the size of v
          size(v, 1) == K1 || error(DimensionMismatch("check the size of v"))
          size(v, 2) == K2 || error(DimensionMismatch("check the size of v"))

          # reverse v and padding zeros
          v_hat = reverse_order(v; n1=dims[1], n2=dims[2])
          v_hat = fft(v_hat)

          # element-wise multiplication
          r = ifft(d_hat .* v_hat)

          # return the last L1 element
          return r[K1:dims[1],K2:dims[2]]

       elseif flag == "adjoint"

          # check the size of v
          size(v, 1) == L1 || error(DimensionMismatch("check the size of v"))
          size(v, 2) == L2 || error(DimensionMismatch("check the size of v"))

          # reverse v and padding zeros
          v_hat = reverse_order(v; n1=dims[1], n2=dims[2])
          v_hat = fft(v_hat)

          # conjugate property
          d_tilde = copy(d_hat)
          for i2 = 1 : dims[2]
              j2 = i2 == 1 ? 1 : dims[2]-i2+2

              for i1 = 1 : dims[1]
                  j1 = i1 == 1 ? 1 : dims[1]-i1+2

                  d_tilde[i1,i2] = conj(d_hat[j1,j2])
              end
          end

          # element-wise multiplication
          r = ifft(d_tilde .* v_hat)

          # return the last K1 element
          return r[L1:dims[1], L2:dims[2]]

       else
          error("non-surpported operation")
       end


    elseif N == 3
       # compute L, K
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
       L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1

       if flag == "forward"

          # check the size of v
          size(v, 1) == K1 || error(DimensionMismatch("check the size of v"))
          size(v, 2) == K2 || error(DimensionMismatch("check the size of v"))
          size(v, 3) == K3 || error(DimensionMismatch("check the size of v"))

          # reverse v and padding zeros
          v_hat = reverse_order(v; n1=dims[1], n2=dims[2], n3=dims[3])
          v_hat = fft(v_hat)

          # element-wise multiplication
          r = ifft(d_hat .* v_hat)

          # return the last L1 element
          return r[K1:dims[1],K2:dims[2],K3:dims[3]]

       elseif flag == "adjoint"

          # check the size of v
          size(v, 1) == L1 || error(DimensionMismatch("check the size of v"))
          size(v, 2) == L2 || error(DimensionMismatch("check the size of v"))
          size(v, 3) == L3 || error(DimensionMismatch("check the size of v"))

          # reverse v and padding zeros
          v_hat = reverse_order(v; n1=dims[1], n2=dims[2], n3=dims[3])
          v_hat = fft(v_hat)

          # conjugate property
          d_tilde = copy(d_hat)
          for i3 = 1 : dims[3]
              j3 = i3 == 1 ? 1 : dims[3]-i3+2

              for i2 = 1 : dims[2]
                  j2 = i2 == 1 ? 1 : dims[2]-i2+2

                  for i1 = 1 : dims[1]
                      j1 = i1 == 1 ? 1 : dims[1]-i1+2

                      d_tilde[i1,i2,i3] = conj(d_hat[j1,j2,j3])
                  end
              end
          end

          # element-wise multiplication
          r = ifft(d_tilde .* v_hat)

          # return the last K1 element
          return r[L1:dims[1],L2:dims[2],L3:dims[3]]

       else
          error("non-surpported operation")
       end


    elseif N == 4
       # compute L, K
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
       L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1
       L4 = floor(Int64,dims[4]/2)+1; K4 = dims[4]-L4+1

       if flag == "forward"

          # check the size of v
          size(v, 1) == K1 || error(DimensionMismatch("check the size of v"))
          size(v, 2) == K2 || error(DimensionMismatch("check the size of v"))
          size(v, 3) == K3 || error(DimensionMismatch("check the size of v"))
          size(v, 4) == K4 || error(DimensionMismatch("check the size of v"))

          # reverse v and padding zeros
          v_hat = reverse_order(v; n1=dims[1], n2=dims[2], n3=dims[3], n4=dims[4])
          v_hat = fft(v_hat)

          # element-wise multiplication
          r = ifft(d_hat .* v_hat)

          # return the last L1 element
          return r[K1:dims[1],K2:dims[2],K3:dims[3],K4:dims[4]]

       elseif flag == "adjoint"

          # check the size of v
          size(v, 1) == L1 || error(DimensionMismatch("check the size of v"))
          size(v, 2) == L2 || error(DimensionMismatch("check the size of v"))
          size(v, 3) == L3 || error(DimensionMismatch("check the size of v"))
          size(v, 4) == L4 || error(DimensionMismatch("check the size of v"))

          # reverse v and padding zeros
          v_hat = reverse_order(v; n1=dims[1], n2=dims[2], n3=dims[3], n4=dims[4])
          v_hat = fft(v_hat)

          # conjugate property
          d_tilde = copy(d_hat)
          for i4 = 1 : dims[4]
              j4 = i4 == 1 ? 1 : dims[4]-i4+2

              for i3 = 1 : dims[3]
                  j3 = i3 == 1 ? 1 : dims[3]-i3+2

                  for i2 = 1 : dims[2]
                      j2 = i2 == 1 ? 1 : dims[2]-i2+2

                      for i1 = 1 : dims[1]
                          j1 = i1 == 1 ? 1 : dims[1]-i1+2

                          d_tilde[i1,i2,i3,i4] = conj(d_hat[j1,j2,j3,j4])
                      end
                  end
              end
          end

          # element-wise multiplication
          r = ifft(d_tilde .* v_hat)

          # return the last K1 element
          return r[L1:dims[1],L2:dims[2],L3:dims[3],L4:dims[4]]

       else
          error("non-surpported operation")
       end


    elseif N == 5
       # compute L, K
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
       L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1
       L4 = floor(Int64,dims[4]/2)+1; K4 = dims[4]-L4+1
       L5 = floor(Int64,dims[5]/2)+1; K5 = dims[5]-L5+1

       if flag == "forward"

          # check the size of v
          size(v, 1) == K1 || error(DimensionMismatch("check the size of v"))
          size(v, 2) == K2 || error(DimensionMismatch("check the size of v"))
          size(v, 3) == K3 || error(DimensionMismatch("check the size of v"))
          size(v, 4) == K4 || error(DimensionMismatch("check the size of v"))
          size(v, 5) == K5 || error(DimensionMismatch("check the size of v"))

          # reverse v and padding zeros
          v_hat = reverse_order(v; n1=dims[1], n2=dims[2], n3=dims[3], n4=dims[4], n5=dims[5])
          v_hat = fft(v_hat)

          # element-wise multiplication
          r = ifft(d_hat .* v_hat)

          # return the last L1 element
          return r[K1:dims[1],K2:dims[2],K3:dims[3],K4:dims[4],K5:dims[5]]

       elseif flag == "adjoint"

          # check the size of v
          size(v, 1) == L1 || error(DimensionMismatch("check the size of v"))
          size(v, 2) == L2 || error(DimensionMismatch("check the size of v"))
          size(v, 3) == L3 || error(DimensionMismatch("check the size of v"))
          size(v, 4) == L4 || error(DimensionMismatch("check the size of v"))
          size(v, 5) == L5 || error(DimensionMismatch("check the size of v"))

          # reverse v and padding zeros
          v_hat = reverse_order(v; n1=dims[1], n2=dims[2], n3=dims[3], n4=dims[4], n5=dims[5])
          v_hat = fft(v_hat)

          # conjugate property
          d_tilde = copy(d_hat)
          for i5 = 1 : dims[5]
              j5 = i5 == 1 ? 1 : dims[5]-i5+2

              for i4 = 1 : dims[4]
                  j4 = i4 == 1 ? 1 : dims[4]-i4+2

                  for i3 = 1 : dims[3]
                      j3 = i3 == 1 ? 1 : dims[3]-i3+2

                      for i2 = 1 : dims[2]
                          j2 = i2 == 1 ? 1 : dims[2]-i2+2

                          for i1 = 1 : dims[1]
                              j1 = i1 == 1 ? 1 : dims[1]-i1+2

                              d_tilde[i1,i2,i3,i4,i5] = conj(d_hat[j1,j2,j3,j4,j5])
                          end
                      end
                  end
              end
          end

          # element-wise multiplication
          r = ifft(d_tilde .* v_hat)

          # return the last K1 element
          return r[L1:dims[1],L2:dims[2],L3:dims[3],L4:dims[4],L5:dims[5]]

      else
         error("only supported up to five dimension")
      end # forward or adjoint

   end # dimension case

end # end function

"""
   Fast anti-diagonal averaging for multi-level Hankel matrix
"""
function anti_diagonal_summation1(u::Vector{Tv}, v::Vector{Tv},
                                  L::Union{Ti,Vector{Ti}},
                                  K::Union{Ti,Vector{Ti}}) where {Tv<:Number, Ti<:Int64}

    # order of hankel matrix
    order =  length(L)
    order == length(K) || error(DimensionMismatch("length of L mismatch length of K"))

    if order == 1

       # check the size of input
       L[1] == length(u) || error(DimensionMismatch("check the length of u"))
       K[1] == length(v) || error(DimensionMismatch("check the length of v"))

       # the size of array
       N1 = L[1]+K[1]-1

       # padding zeros
       if Tv <: Complex
          u_hat = zeros(Tv, N1)
          v_hat = zeros(Tv, N1)
       else
          u_hat = zeros(Complex{Tv}, N1)
          v_hat = zeros(Complex{Tv}, N1)
       end
       u_hat[1:L[1]] .=      u
       v_hat[1:K[1]] .= conj(v)

       # Fourier transform
       fft!(u_hat)
       fft!(v_hat)
       u_hat .= u_hat .* v_hat
       ifft!(u_hat)


    elseif order == 2

       # check the size of input
       prod(L[1:2]) == length(u) || error(DimensionMismatch("check the length of u"))
       prod(K[1:2]) == length(v) || error(DimensionMismatch("check the length of v"))

       # the size of array
       N1 = L[1]+K[1]-1
       N2 = L[2]+K[2]-1

       # padding zeros
       if Tv <: Complex
          u_hat = zeros(Tv, N1, N2)
          v_hat = zeros(Tv, N1, N2)
       else
          u_hat = zeros(Complex{Tv}, N1, N2)
          v_hat = zeros(Complex{Tv}, N1, N2)
       end
       u_hat[1:L[1],1:L[2]] .=      reshape(u, L[1], L[2])
       v_hat[1:K[1],1:K[2]] .= conj(reshape(v, K[1], K[2]))

       # Fourier transform
       fft!(u_hat)
       fft!(v_hat)
       u_hat .= u_hat .* v_hat
       ifft!(u_hat)


    elseif order == 3

       # check the size of input
       prod(L[1:3]) == length(u) || error(DimensionMismatch("check the length of u"))
       prod(K[1:3]) == length(v) || error(DimensionMismatch("check the length of v"))

       # the size of array
       N1 = L[1]+K[1]-1
       N2 = L[2]+K[2]-1
       N3 = L[3]+K[3]-1

       # padding zeros
       if Tv <: Complex
          u_hat = zeros(Tv, N1, N2, N3)
          v_hat = zeros(Tv, N1, N2, N3)
       else
          u_hat = zeros(Complex{Tv}, N1, N2, N3)
          v_hat = zeros(Complex{Tv}, N1, N2, N3)
       end
       u_hat[1:L[1],1:L[2],1:L[3]] .=      reshape(u, L[1], L[2], L[3])
       v_hat[1:K[1],1:K[2],1:K[3]] .= conj(reshape(v, K[1], K[2], K[3]))

       # Fourier transform
       fft!(u_hat)
       fft!(v_hat)
       u_hat .= u_hat .* v_hat
       ifft!(u_hat)


    elseif order == 4

       # check the size of input
       prod(L[1:4]) == length(u) || error(DimensionMismatch("check the length of u"))
       prod(K[1:4]) == length(v) || error(DimensionMismatch("check the length of v"))

       # the size of array
       N1 = L[1]+K[1]-1
       N2 = L[2]+K[2]-1
       N3 = L[3]+K[3]-1
       N4 = L[4]+K[4]-1

       # padding zeros
       if Tv <: Complex
          u_hat = zeros(Tv, N1, N2, N3, N4)
          v_hat = zeros(Tv, N1, N2, N3, N4)
       else
          u_hat = zeros(Complex{Tv}, N1, N2, N3, N4)
          v_hat = zeros(Complex{Tv}, N1, N2, N3, N4)
       end
       u_hat[1:L[1],1:L[2],1:L[3],1:L[4]] .=      reshape(u, L[1], L[2], L[3], L[4])
       v_hat[1:K[1],1:K[2],1:K[3],1:K[4]] .= conj(reshape(v, K[1], K[2], K[3], K[4]))

       # Fourier transform
       fft!(u_hat)
       fft!(v_hat)
       u_hat .= u_hat .* v_hat
       ifft!(u_hat)

    elseif order == 5

       # check the size of input
       prod(L[1:5]) == length(u) || error(DimensionMismatch("check the length of u"))
       prod(K[1:5]) == length(v) || error(DimensionMismatch("check the length of v"))

       # the size of array
       N1 = L[1]+K[1]-1
       N2 = L[2]+K[2]-1
       N3 = L[3]+K[3]-1
       N4 = L[4]+K[4]-1
       N5 = L[5]+K[5]-1

       # padding zeros
       if Tv <: Complex
          u_hat = zeros(Tv, N1, N2, N3, N4, N5)
          v_hat = zeros(Tv, N1, N2, N3, N4, N5)
       else
          u_hat = zeros(Complex{Tv}, N1, N2, N3, N4, N5)
          v_hat = zeros(Complex{Tv}, N1, N2, N3, N4, N5)
       end
       u_hat[1:L[1],1:L[2],1:L[3],1:L[4],1:L[5]] .=      reshape(u, L[1], L[2], L[3], L[4], L[5])
       v_hat[1:K[1],1:K[2],1:K[3],1:K[4],1:K[5]] .= conj(reshape(v, K[1], K[2], K[3], K[4], K[5]))

       # Fourier transform
       fft!(u_hat)
       fft!(v_hat)
       u_hat .= u_hat .* v_hat
       ifft!(u_hat)

    else
       error("only support up-to fiveth order")
    end

    # convert the output to the same type as input
    if Tv <: Complex
       return u_hat
    else
       return real(u_hat)
    end

end


"""
   anti-diagonal average via building multi-level hankel matrix first.
"""
function anti_diagonal_summation_slow(u::Vector{Tv}, v::Vector{Tv},
                                      L::Union{Ti,Vector{Ti}},
                                      K::Union{Ti,Vector{Ti}}) where {Tv <: Number, Ti<:Int64}

    # order of hankel matrix
    order =  length(L)
    order == length(K) || error(DimensionMismatch("length of L mismatch length of K"))

    # outer product u * v'
    H = u * v'

    if order == 1

       # check the size of input
       L[1] == length(u) || error(DimensionMismatch("check the length of u"))
       K[1] == length(v) || error(DimensionMismatch("check the length of v"))

       # the size of array
       N1 = L[1]+K[1]-1

       # allocate memory for the result
       d  = zeros(Tv,N1)

       # first layer
       for j1 = 1 : K[1]
           for i1 = 1 : L[1]
               n1 = i1 + j1 - 1

               d[n1] += H[i1,j1]
           end
       end #first


    elseif order == 2

       # check the size of input
       prod(L[1:2]) == length(u) || error(DimensionMismatch("check the length of u"))
       prod(K[1:2]) == length(v) || error(DimensionMismatch("check the length of v"))

       # the size of array
       N1 = L[1]+K[1]-1
       N2 = L[2]+K[2]-1

       # allocate memory for the result
       d  = zeros(Tv, N1, N2)

       # second layer
       for j2 = 1 : K[2]
           c2 = (j2-1)*K[1]
           for i2 = 1 : L[2]
               r2 = (i2-1)*L[1]
               n2 = i2 + j2 - 1

               # first layer
               for j1 = 1 : K[1]
                   c1 = c2 + j1
                   for i1 = 1 : L[1]
                       r1 = r2 + i1
                       n1 = i1 + j1 - 1

                       d[n1,n2] += H[r1,c1]
                   end
               end #first
           end
       end #second


    elseif order == 3

       # check the size of input
       prod(L[1:3]) == length(u) || error(DimensionMismatch("check the length of u"))
       prod(K[1:3]) == length(v) || error(DimensionMismatch("check the length of v"))

       # the size of array
       N1 = L[1]+K[1]-1
       N2 = L[2]+K[2]-1
       N3 = L[3]+K[3]-1

       # allocate memory for the result
       d  = zeros(Tv, N1, N2, N3)

       # anti-diagonal summation
       for j3 = 1 : K[3]
           c3 = (j3-1)*K[2]*K[1]
           for i3 = 1 : L[3]
               r3 = (i3-1)*L[2]*L[1]
               n3 = i3 + j3 - 1

               # second layer
               for j2 = 1 : K[2]
                   c2 = (j2-1)*K[1]
                   for i2 = 1 : L[2]
                       r2 = (i2-1)*L[1]
                       n2 = i2 + j2 - 1

                       # first layer
                       for j1 = 1 : K[1]
                           c1 = c3 + c2 + j1
                           for i1 = 1 : L[1]
                               r1 = r3 + r2 + i1
                               n1 = i1 + j1 - 1

                               d[n1,n2,n3] += H[r1,c1]
                           end
                       end #first
                   end
               end #second
           end
       end #thirds


    elseif order == 4

       # check the size of input
       prod(L[1:4]) == length(u) || error(DimensionMismatch("check the length of u"))
       prod(K[1:4]) == length(v) || error(DimensionMismatch("check the length of v"))

       # the size of array
       N1 = L[1]+K[1]-1
       N2 = L[2]+K[2]-1
       N3 = L[3]+K[3]-1
       N4 = L[4]+K[4]-1

       # allocate memory for the result
       d  = zeros(Tv, N1, N2, N3, N4)

       # anti-diagonal summation
       for j4 = 1 : K[4]
           c4 = (j4-1)*K[3]*K[2]*K[1]
           for i4 = 1 : L[4]
               r4 = (i4-1)*L[3]*L[2]*L[1]
               n4 = i4 + j4 - 1

               # third layer
               for j3 = 1 : K[3]
                   c3 = (j3-1)*K[2]*K[1]
                   for i3 = 1 : L[3]
                       r3 = (i3-1)*L[2]*L[1]
                       n3 = i3 + j3 - 1

                       # second layer
                       for j2 = 1 : K[2]
                           c2 = (j2-1)*K[1]
                           for i2 = 1 : L[2]
                               r2 = (i2-1)*L[1]
                               n2 = i2 + j2 - 1

                               # first layer
                               for j1 = 1 : K[1]
                                   c1 = c4 + c3 + c2 + j1
                                   for i1 = 1 : L[1]
                                       r1 = r4 + r3 + r2 + i1
                                       n1 = i1 + j1 - 1

                                       d[n1,n2,n3,n4] += H[r1,c1]
                                   end
                               end #first
                           end
                       end #second
                   end
               end #third
           end
       end #fourth


    elseif order == 5

       # check the size of input
       prod(L[1:5]) == length(u) || error(DimensionMismatch("check the length of u"))
       prod(K[1:5]) == length(v) || error(DimensionMismatch("check the length of v"))

       # the size of array
       N1 = L[1]+K[1]-1
       N2 = L[2]+K[2]-1
       N3 = L[3]+K[3]-1
       N4 = L[4]+K[4]-1
       N5 = L[5]+K[5]-1

       # allocate memory for the result
       d  = zeros(Tv, N1, N2, N3, N4, N5)

       # anti-diagonal summation
       for j5 = 1 : K[5]
           c5 = (j5-1)*K[4]*K[3]*K[2]*K[1]
           for i5 = 1 : L[5]
               r5 = (i5-1)*L[4]*L[3]*L[2]*L[1]
               n5 = i5 + j5 - 1

               for j4 = 1 : K[4]
                   c4 = (j4-1)*K[3]*K[2]*K[1]
                   for i4 = 1 : L[4]
                       r4 = (i4-1)*L[3]*L[2]*L[1]
                       n4 = i4 + j4 - 1

                       # third layer
                       for j3 = 1 : K[3]
                           c3 = (j3-1)*K[2]*K[1]
                           for i3 = 1 : L[3]
                               r3 = (i3-1)*L[2]*L[1]
                               n3 = i3 + j3 - 1

                               # second layer
                               for j2 = 1 : K[2]
                                   c2 = (j2-1)*K[1]
                                   for i2 = 1 : L[2]
                                       r2 = (i2-1)*L[1]
                                       n2 = i2 + j2 - 1

                                       # first layer
                                       for j1 = 1 : K[1]
                                           c1 = c5 + c4 + c3 + c2 + j1
                                           for i1 = 1 : L[1]
                                               r1 = r5 + r4 + r3 + r2 + i1
                                               n1 = i1 + j1 - 1

                                               d[n1,n2,n3,n4,n5] += H[r1,c1]
                                           end
                                       end #first
                                   end
                               end #second
                           end
                       end #third
                   end
               end #fourth
           end
       end # fiveth


    else
       error("only support up-to fiveth order")
    end

    return d
end

Tv= Complex{Float32}

# 1D
L = [12];
K = [11];
u = randn(Tv, prod(L[1:1]));
v = randn(Tv, prod(K[1:1]));
@time d1 = anti_diagonal_summation1(u, v, L, K);
@time d2 = anti_diagonal_summation2(u, v, L, K);
norm(d1-d2) / norm(d1)

# 2D
L = [8,9];
K = [7,8];
u = randn(Tv, prod(L[1:2]));
v = randn(Tv, prod(K[1:2]));
@time d1 = anti_diagonal_summation1(u, v, L, K);
@time d2 = anti_diagonal_summation2(u, v, L, K);
norm(d1-d2) / norm(d1)

# 3D
L = [7,8,9];
K = [6,7,8];
u = randn(Tv, prod(L[1:3]));
v = randn(Tv, prod(K[1:3]));
@time d1 = anti_diagonal_summation1(u, v, L, K);
@time d2 = anti_diagonal_summation2(u, v, L, K);
norm(d1-d2) / norm(d1)

# 4D
L = [7,8,9,10];
K = [6,7,8,9 ];
u = randn(Tv, prod(L[1:4]));
v = randn(Tv, prod(K[1:4]));
@time d1 = anti_diagonal_summation1(u, v, L, K);
@time d2 = anti_diagonal_summation2(u, v, L, K);
norm(d1-d2) / norm(d1)

# 5D
L = [6,7,8,9,10];
K = [5,6,7,8,9 ];
u = randn(Tv, prod(L[1:5]));
v = randn(Tv, prod(K[1:5]));
@time d1 = anti_diagonal_summation1(u, v, L, K);
@time d2 = anti_diagonal_summation2(u, v, L, K);
norm(d1-d2) / norm(d1)



#
# # forward
# r1 = H * vec(v);
# r2 = hankel_multiplication(d, v; flag="forward");
# r2 = vec(r2);
# norm(r1-r2) / norm(r1)
#
# # adjoint
# v = randn(Tv, L1, L2, L3, L4, L5);
# r3 = H' * vec(v);
# r4 = hankel_multiplication(d, v; flag="adjoint");
# r4 = vec(r4);
# norm(r3-r4) / norm(r3)

# function count_hankel1(dims::Union{Ti,Vector{Ti}}) where {Ti<:Int64}
#
#     # convert to a vector when the input is a scalar
#     if typeof(dims) == Int64
#        dims = [dims]
#     end
#
#     # determine the dimensions of an array
#     N = length(dims)
#
#     # level 1
#     if N == 1
#
#        # allocate memory to record the number
#        count_num = zeros(Int64, dims[1])
#
#        # size of hankel matrix
#        L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
#
#
#        # build hankel
#        for j1 = 1 : K1
#            for i1 = 1 : L1
#                n1 = i1+j1-1
#
#                count_num[n1] += 1
#            end
#        end
#
#
#     # level 2
#     elseif N == 2
#
#        # allocate memory to record the number
#        count_num = zeros(Int64, dims[1], dims[2])
#
#        # size of hankel matrix
#        L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
#        L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
#
#        # build Hankel
#        for j2 = 1 : K2
#            c2 = (j2-1)*K1
#            for i2 = 1 : L2
#                r2 = (i2-1)*L1
#                n2 = i2 + j2 - 1
#
#                # first layer
#                for j1 = 1 : K1
#                    c1 = c2 + j1
#                    for i1 = 1 : L1
#                        r1 = r2 + i1
#                        n1 = i1 + j1 - 1
#
#                        count_num[n1,n2] += 1
#                    end
#                end #first
#            end
#        end #second
#
#
#     # level 3
#     elseif N == 3
#
#        # allocate memory to record the number
#        count_num = zeros(Int64, dims[1], dims[2], dims[3])
#
#        # size of hankel matrix
#        L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
#        L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
#        L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1
#
#        # build Hankel
#        for j3 = 1 : K3
#            c3 = (j3-1)*K2*K1
#            for i3 = 1 : L3
#                r3 = (i3-1)*L2*L1
#                n3 = i3 + j3 - 1
#
#                # second layer
#                for j2 = 1 : K2
#                    c2 = (j2-1)*K1
#                    for i2 = 1 : L2
#                        r2 = (i2-1)*L1
#                        n2 = i2 + j2 - 1
#
#                        # first layer
#                        for j1 = 1 : K1
#                            c1 = c3 + c2 + j1
#                            for i1 = 1 : L1
#                                r1 = r3 + r2 + i1
#                                n1 = i1 + j1 - 1
#
#                                count_num[n1,n2,n3] += 1
#                            end
#                        end #first
#                    end
#                end #second
#            end
#        end #third
#
#
#     # level 4
#     elseif N == 4
#
#        # allocate memory to record the number
#        count_num = zeros(Int64, dims[1], dims[2], dims[3], dims[4])
#
#        # size of hankel matrix
#        L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
#        L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
#        L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1
#        L4 = floor(Int64,dims[4]/2)+1; K4 = dims[4]-L4+1
#
#        # build Hankel
#        for j4 = 1 : K4
#            c4 = (j4-1)*K3*K2*K1
#            for i4 = 1 : L4
#                r4 = (i4-1)*L3*L2*L1
#                n4 = i4 + j4 - 1
#
#                # third layer
#                for j3 = 1 : K3
#                    c3 = (j3-1)*K2*K1
#                    for i3 = 1 : L3
#                        r3 = (i3-1)*L2*L1
#                        n3 = i3 + j3 - 1
#
#                        # second layer
#                        for j2 = 1 : K2
#                            c2 = (j2-1)*K1
#                            for i2 = 1 : L2
#                                r2 = (i2-1)*L1
#                                n2 = i2 + j2 - 1
#
#                                # first layer
#                                for j1 = 1 : K1
#                                    c1 = c4 + c3 + c2 + j1
#                                    for i1 = 1 : L1
#                                        r1 = r4 + r3 + r2 + i1
#                                        n1 = i1 + j1 - 1
#
#                                        count_num[n1,n2,n3,n4] += 1
#                                    end
#                                end #first
#                            end
#                        end #second
#                    end
#                end #third
#            end
#        end #fourth
#
#
#     # level 5
#     elseif N == 5
#
#        # allocate memory to record the number
#        count_num = zeros(Int64, dims[1], dims[2], dims[3], dims[4], dims[5])
#
#        # size of hankel matrix
#        L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
#        L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
#        L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1
#        L4 = floor(Int64,dims[4]/2)+1; K4 = dims[4]-L4+1
#        L5 = floor(Int64,dims[5]/2)+1; K5 = dims[5]-L5+1
#
#        # build Hankel
#        for j5 = 1 : K5
#            c5 = (j5-1)*K4*K3*K2*K1
#            for i5 = 1 : L5
#                r5 = (i5-1)*L4*L3*L2*L1
#                n5 = i5 + j5 - 1
#
#                for j4 = 1 : K4
#                    c4 = (j4-1)*K3*K2*K1
#                    for i4 = 1 : L4
#                        r4 = (i4-1)*L3*L2*L1
#                        n4 = i4 + j4 - 1
#
#                        # third layer
#                        for j3 = 1 : K3
#                            c3 = (j3-1)*K2*K1
#                            for i3 = 1 : L3
#                                r3 = (i3-1)*L2*L1
#                                n3 = i3 + j3 - 1
#
#                                # second layer
#                                for j2 = 1 : K2
#                                    c2 = (j2-1)*K1
#                                    for i2 = 1 : L2
#                                        r2 = (i2-1)*L1
#                                        n2 = i2 + j2 - 1
#
#                                        # first layer
#                                        for j1 = 1 : K1
#                                            c1 = c5 + c4 + c3 + c2 + j1
#                                            for i1 = 1 : L1
#                                                r1 = r5 + r4 + r3 + r2 + i1
#                                                n1 = i1 + j1 - 1
#
#                                                count_num[n1,n2,n3,n4,n5] += 1
#                                            end
#                                        end #first
#                                    end
#                                end #second
#                            end
#                        end #third
#                    end
#                end #fourth
#            end
#        end # fiveth
#
#
#     else
#        error("only support up to 5 dimensions")
#     end
#
#     return count_num
# end
