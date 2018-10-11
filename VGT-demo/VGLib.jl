module VGLib
using PyCall,LsqFit
@pyimport numpy as np
@pyimport scipy.interpolate as si
sir=si.RectBivariateSpline;
typealias Cube{T} Array{T,3}
typealias Mat{T} Array{T,2}
typealias Vec{T} Array{T,1}
typealias Vecx{T} LinSpace{T}

export sban2d,getIC_PPV,sobel_conv_2d


function sban2d(Ax::Mat,Ay::Mat,dn)
	nx,ny=size(Ax)
	Ana=zeros(div(nx,dn),div(ny,dn));
	Ans=zeros(div(nx,dn),div(ny,dn));

	for  j in 1:div(ny,dn),i in 1:div(nx,dn) #try conv!
	   is=(i-1)*dn+1;
	   ie=i*dn;
	   js=(j-1)*dn+1;
	   je=j*dn;
	   Axx=Ax[is:ie,js:je];
	   Ayy=Ay[is:ie,js:je];
	   binsize=dn;
	   Apeak,Adisp=fit_gaussian_2do(Axx,Ayy,binsize);
	   Ana[i,j]=Apeak;
	   Ans[i,j]=Adisp;
	end
	return Ana,Ans
end

function fit_gaussian_2do(Ax::Mat,Ay::Mat,binsize)
	phi=atan(Ay./Ax)
	Gaus(x,p)=p[1]*exp(-(x-p[2]).^2*p[3]);
	ax,ac=hist(phi[~isnan(phi)][:],linspace(-pi/2,pi/2,binsize+1))
	ax=.5*(ax[1:end-1]+ax[2:end]);
	if (abs(ax[maxid(ac)])[1]<pi/4)
	fit=curve_fit(Gaus,ax,ac/sum(ac),[maximum(ac/sum(ac)),0.0,1.0])
	else
	ax=ax-pi/2;
	ac=fftshift(ac);
	fit=curve_fit(Gaus,ax,ac/sum(ac),[maximum(ac/sum(ac)),-pi/2,1.0]);
	end
	sigma=estimate_errors(fit,0.95);
	return fit.param[2],sigma[2];
end

function maxid(ax::Vec)
 return find(ax.==maximum(ax));
end

function sobel_conv_2d(A::Mat)
 Kx,Ky=sobel_kernel_2d(A);
 Ax=convoluting_kernel(A,Kx);
 Ay=convoluting_kernel(A,Ky);
 return Ax,Ay
end


function sobel_kernel_2d(A::Mat)
 nx,ny=size(A);
 Ax=zeros(size(A));
 Ay=zeros(size(A));
 vp=sobel_parallel(3);
 vl=sobel_perpendicular(3);
 Axx=zeros(3,3);
 Ayy=zeros(3,3);
 for j in 1:3, i in 1:3
  Axx[i,j]=vp[i]*vl[j];
  Ayy[i,j]=vl[i]*vp[j];
 end
 Ax[1:3,1:3]=circshift(Axx,(1,1));
 Ay[1:3,1:3]=circshift(Ayy,(1,1));
 Ax=circshift(Ax,(-1,-1));
 Ay=circshift(Ay,(-1,-1));
 return Ax,Ay
end

function sobel_perpendicular(size::Int)
 if (size>=3)
  v=zeros(size);
  v[2]=1
  v[1]=2
  v[end]=1
  return v
 else
  println("Size<3 not supported")
  return 0;
 end
end

function sobel_parallel(size::Int)
 if (size>=3)
  v=zeros(size);
  v[2]=-1
  v[end]=1
  return v
 else
  println("Size<3 not supported")
  return 0;
 end
end

function convoluting_kernel(A::Mat,B::Mat)
 Af=fft(A);
 Bf=fft(B);
 Cf=Af.*Bf
 C=real(ifft(Cf));
 return C
end

function getIC_PPV(d::Cube,v::Vec)
	nx,ny,nv=size(d);
	I,C=zeros((nx,ny)),zeros((nx,ny));
	for k in 1:nv, j in 1:ny, i in 1:nx
		I[i,j]+=d[i,j,k];
		C[i,j]+=d[i,j,k].*v[k];
	end
	C=C./I;
    return I,C
end


end