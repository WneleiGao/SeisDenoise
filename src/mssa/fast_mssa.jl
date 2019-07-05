"""
   build a circulant matrix from a vector
"""
function build_circulant_matrix(v::Vector{Tv}) where {Tv<:Number}

    # length of the vector
    N = length(v)

    # allocate space for circulant matrix N × N
    C = Matrix{Tv}(undef, N, N)

    # loop over column
    for i2 = 1 : N

        # loop over row, first part
        for i1 = 1 : i2-1
            idx = N - i2 + i1 + 1
            C[i1,i2] = v[idx]
        end

        # loop over row, second part
        for i1 = i2 : N
            idx = i1 - i2 + 1
            C[i1,i2] = v[idx]
        end
    end

    return C
end

"""
   build a Toeplitze matrix from a vector, this matrix is a square (length of v is odd)
or one more row than columns (length of v is even)
"""
function build_toeplitz_matrix(v::Vector{Tv}) where {Tv <: Number}

    # length of the vector
    N = length(v)

    # compute L, K
    L = floor(Int64, N/2) + 1
    K = N + 1 - L

    # allocate space for Toeplitze matrix
    T = Matrix{Tv}(undef, L, K)

    # loop over column
    for i2 = 1 : K

        # loop over row
        for i1 = 1 : L
            idx = K + i1 - i2
            T[i1,i2] = v[idx]
        end
    end

    return T
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
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[1]-L2+1

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
          size(v, 1) == K1 ||  DimensionMismatch()

          # reverse v and padding zeros
          v_hat = reverse_order(v; n1=dims[1])
          v_hat = fft(v_hat)

          # element-wise multiplication
          r = ifft(d_hat .* v_hat)

          # return the last L1 element
          return r[K1:dims[1]]

       elseif flag == "adjoint"

          # check the size of v
          size(v, 1) == L1 ||  DimensionMismatch()

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
          size(v, 1) == K1 || DimensionMismatch()
          size(v, 2) == K2 || DimensionMismatch()

          # reverse v and padding zeros
          v_hat = reverse_order(v; n1=dims[1], n2=dims[2])
          v_hat = fft(v_hat)

          # element-wise multiplication
          r = ifft(d_hat .* v_hat)

          # return the last L1 element
          return r[K1:dims[1],K2:dims[2]]

       elseif flag == "adjoint"

          # check the size of v
          size(v, 1) == L1 || DimensionMismatch()
          size(v, 2) == L2 || DimensionMismatch()

          # reverse v and padding zeros
          v_hat = reverse_order(v; n1=dims[1], n2=dims[2])
          v_hat = fft(v_hat)

          # conjugate property
          d_tilde = copy(d_hat)
          for i2 = 1 : dims[2]
              j2 = i2 == 1 ? 1 : dims[2]-i2+2

              for i1 = 1 : dims[1]
                  j1 = i1 == 1 ? 1 : dims[1]-i1+2

                  d_tilde[i1,i2,i3,i4,i5] = conj(d_hat[j1,j2,j3,j4,j5])
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
          size(v, 1) == K1 || DimensionMismatch()
          size(v, 2) == K2 || DimensionMismatch()
          size(v, 3) == K3 || DimensionMismatch()

          # reverse v and padding zeros
          v_hat = reverse_order(v; n1=dims[1], n2=dims[2], n3=dims[3])
          v_hat = fft(v_hat)

          # element-wise multiplication
          r = ifft(d_hat .* v_hat)

          # return the last L1 element
          return r[K1:dims[1],K2:dims[2],K3:dims[3]]

       elseif flag == "adjoint"

          # check the size of v
          size(v, 1) == L1 || DimensionMismatch()
          size(v, 2) == L2 || DimensionMismatch()
          size(v, 3) == L3 || DimensionMismatch()

          # reverse v and padding zeros
          v_hat = reverse_order(v; n1=dims[1], n2=dims[2], n3=dims[3])
          v_hat = fft(v_hat)

          # conjugate property
          d_tilde = copy(d_hat)
          for i3 = 1 : dims[1]
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
          size(v, 1) == K1 || DimensionMismatch()
          size(v, 2) == K2 || DimensionMismatch()
          size(v, 3) == K3 || DimensionMismatch()
          size(v, 4) == K4 || DimensionMismatch()

          # reverse v and padding zeros
          v_hat = reverse_order(v; n1=dims[1], n2=dims[2], n3=dims[3], n4=dims[4])
          v_hat = fft(v_hat)

          # element-wise multiplication
          r = ifft(d_hat .* v_hat)

          # return the last L1 element
          return r[K1:dims[1],K2:dims[2],K3:dims[3],K4:dims[4]]

       elseif flag == "adjoint"

          # check the size of v
          size(v, 1) == L1 || DimensionMismatch()
          size(v, 2) == L2 || DimensionMismatch()
          size(v, 3) == L3 || DimensionMismatch()
          size(v, 4) == L4 || DimensionMismatch()

          # reverse v and padding zeros
          v_hat = reverse_order(v; n1=dims[1], n2=dims[2], n3=dims[3], n4=dims[4])
          v_hat = fft(v_hat)

          # conjugate property
          d_tilde = copy(d_hat)
          for i4 = 1 : dims[4]
              j4 = i4 == 1 ? 1 : dims[4]-i4+2

              for i3 = 1 : dims[1]
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
          size(v, 1) == K1 || DimensionMismatch()
          size(v, 2) == K2 || DimensionMismatch()
          size(v, 3) == K3 || DimensionMismatch()
          size(v, 4) == K4 || DimensionMismatch()
          size(v, 5) == K5 || DimensionMismatch()

          # reverse v and padding zeros
          v_hat = reverse_order(v; n1=dims[1], n2=dims[2], n3=dims[3], n4=dims[4], n5=dims[5])
          v_hat = fft(v_hat)

          # element-wise multiplication
          r = ifft(d_hat .* v_hat)

          # return the last L1 element
          return r[K1:dims[1],K2:dims[2],K3:dims[3],K4:dims[4],K5:dims[5]]

       elseif flag == "adjoint"

          # check the size of v
          size(v, 1) == L1 || DimensionMismatch()
          size(v, 2) == L2 || DimensionMismatch()
          size(v, 3) == L3 || DimensionMismatch()
          size(v, 4) == L4 || DimensionMismatch()
          size(v, 5) == L5 || DimensionMismatch()

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
    end
end

# Hankel matrix or its Hermitian transpose times a vector
function hankel_times_vector(c::Vector{Tv}, v::Vector{Tv}; flag="forward") where {Tv<:Number}

    # length of vector
    N = length(c)

    # compute L, K
    L = floor(Int64, N/2) + 1
    K = N + 1 - L

    # reverse the order of v
    v_hat = reverse(v)

    # Fourier transform of c
    c_hat = fft(c)

    if flag == "forward"
       # check the length of v
       length(v) == K ||  DimensionMismatch()

       # padding zeros
       v_hat = vcat(v_hat, zeros(Tv, L-1))
       v_hat = fft(v_hat)

       # element-wise multiplication
       r = ifft(c_hat .* v_hat)

       # return the last L element
       return r[K:N]

    elseif flag == "adjoint"
       # check the length of v
       length(v) == L ||  DimensionMismatch()

       # padding zeros
       v_hat = vcat(v_hat, zeros(Tv, K-1))
       v_hat = fft(v_hat)

       # conjugate property
       c_tilde = copy(c_hat)
       c_tilde[1] = conj(c_hat[1])
       for i = 2 : N
           c_tilde[i] = conj(c_hat[N-i+2])
       end

       # element-wise multiplication
       r = ifft(c_tilde .* v_hat)

       # return the last L element
       return r[L:N]

    else
       error("non-surpported operation")
    end
end


# # test fast multiplication
# N = 112;
# s = randn(N);
# v = randn(N);
# C = build_circulant_matrix(s);
#
# r1 = C * v
# r2 = real(ifft( fft(s) .* fft(v) ))
# norm(r1-r2) / norm(r1)
#
# # test teoplitz matrix times a vector
# L = 51; K = 50;
# N = L+K-1;
# c = randn(N); v = randn(K);
# r = toeplitz_times_vector(c, v; embed_flag=2);
#
# T  = build_toeplitz_matrix(c);
# r1 = T * v;
# norm(r-r1) / norm(r)


# test Hankel matrix times a vector
# Tv = Complex{Float64}
# L = 51; K = 51; N = L+K-1;
#
# c = randn(Tv, N); v = randn(Tv, K);
# r1 = hankel_times_vector(c, v; flag="forward");
#
# H  = build_hankel_matrix(c);
# r2 = H * v;
# norm(r1-r2) / norm(r1)
#
# v = randn(Tv, L);
# r3 = hankel_times_vector(c, v; flag="adjoint");
#
# r4 = H' * v;
# norm(r3-r4) / norm(r3)


# # test the property of Fourier transform of complex conjugate of a vector
using FFTW, Algebra
N1 = 100
N2 = 122
a = rand(Complex{Float64}, N1, N2)
A = fft(a)

b = conj(a)
B = fft(b)

C = zeros(Complex{Float64}, N1, N2)
for i2 = 1 : N2
    j2 = i2 == 1 ? 1 : N2-i2+2

    for i1 = 1 : N1
        j1 = i1 == 1 ? 1 : N1-i1+2
        C[i1,i2] = conj(A[j1,j2])
    end
end

norm(C-B) / norm(C)
