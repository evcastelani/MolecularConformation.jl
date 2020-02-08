"""
  function align(x,y) -> xnew
  aligns two structures [sets of points in 3D space]. Solves
  the "Procrustes" problem. Structures are expected to be of the same size, and the 
  correspondence is assumed from the vector indices. 
 
  Returns x aligned, by performing the rigid body transformation [rotation
  and translation that minimizes the RMSD between x and y].
 
  x, y, and xnew (return) are matrices of dimensions (n,3) 
  (n is the number of points, 3 is the dimension of the space).
 
  L. Martinez, Institute of Chemistry - University of Campinas
  Jan 04, 2019
"""

using LinearAlgebra

function align( x :: Matrix{Float64}, y :: Matrix{Float64} )

  n = size(x,1)

  # Computing centroid

  cmx = zeros(3)
  cmy = zeros(3)
  for i in 1:n
    for j in 1:3
      cmx[j] = cmx[j] + x[i,j]
      cmy[j] = cmy[j] + y[i,j]
    end
  end
  cmx = cmx / n
  cmy = cmy / n

  # Translating both sets to the origin

  for i in 1:n
    for j in 1:3
      x[i,j] = x[i,j] - cmx[j]
      y[i,j] = y[i,j] - cmy[j]
    end
  end

  # Computing the quaternion matrix

  xm = Vector{Float64}(undef,n)
  ym = Vector{Float64}(undef,n)
  zm = Vector{Float64}(undef,n)
  xp = Vector{Float64}(undef,n)
  yp = Vector{Float64}(undef,n)
  zp = Vector{Float64}(undef,n)
  for i in 1:n
    xm[i] = y[i,1] - x[i,1]
    ym[i] = y[i,2] - x[i,2]
    zm[i] = y[i,3] - x[i,3]
    xp[i] = y[i,1] + x[i,1]
    yp[i] = y[i,2] + x[i,2]
    zp[i] = y[i,3] + x[i,3]
  end

  q = zeros(4,4)
  for i in 1:n
    q[1,1] = q[1,1] + xm[i]^2 + ym[i]^2 + zm[i]^2
    q[1,2] = q[1,2] + yp[i]*zm[i] - ym[i]*zp[i]
    q[1,3] = q[1,3] + xm[i]*zp[i] - xp[i]*zm[i]
    q[1,4] = q[1,4] + xp[i]*ym[i] - xm[i]*yp[i]
    q[2,2] = q[2,2] + yp[i]^2 + zp[i]^2 + xm[i]^2
    q[2,3] = q[2,3] + xm[i]*ym[i] - xp[i]*yp[i]
    q[2,4] = q[2,4] + xm[i]*zm[i] - xp[i]*zp[i]
    q[3,3] = q[3,3] + xp[i]^2 + zp[i]^2 + ym[i]^2
    q[3,4] = q[3,4] + ym[i]*zm[i] - yp[i]*zp[i]
    q[4,4] = q[4,4] + xp[i]^2 + yp[i]^2 + zm[i]^2
  end
  q[2,1] = q[1,2]
  q[3,1] = q[1,3]
  q[3,2] = q[2,3]
  q[4,1] = q[1,4]
  q[4,2] = q[2,4]
  q[4,3] = q[3,4]          

  # Computing the eigenvectors 'v' of the q matrix

  v = LinearAlgebra.eigvecs(q)

  # Compute rotation matrix
  
  u = Matrix{Float64}(undef,3,3)
  u[1,1] = v[1,1]^2 + v[2,1]^2 - v[3,1]^2 - v[4,1]^2
  u[1,2] = 2. * ( v[2,1]*v[3,1] + v[1,1]*v[4,1] )
  u[1,3] = 2. * ( v[2,1]*v[4,1] - v[1,1]*v[3,1] )
  u[2,1] = 2. * ( v[2,1]*v[3,1] - v[1,1]*v[4,1] )
  u[2,2] = v[1,1]^2 + v[3,1]^2 - v[2,1]^2 - v[4,1]^2
  u[2,3] = 2. * ( v[3,1]*v[4,1] + v[1,1]*v[2,1] )
  u[3,1] = 2. * ( v[2,1]*v[4,1] + v[1,1]*v[3,1] )
  u[3,2] = 2. * ( v[3,1]*v[4,1] - v[1,1]*v[2,1] )
  u[3,3] = v[1,1]^2 + v[4,1]^2 - v[2,1]^2 - v[3,1]^2      

  # Rotate vector x [will be stored in xnew], and restore y

  xnew = zeros(n,3)
  for i in 1:n
    for j in 1:3
      for k in 1:3
        xnew[i,j] = xnew[i,j] + u[j,k] * x[i,k]
      end 
    end
  end

  # Translate vector to the centroid of y [and restore x and y]

  for i in 1:n
    for j in 1:3
      xnew[i,j] = xnew[i,j] + cmy[j]
      y[i,j] = y[i,j] + cmy[j]
      x[i,j] = x[i,j] + cmx[j]
    end
  end

  return xnew

end

