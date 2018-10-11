module CL03decompo
using StatsBase,HDF5, LsqFit, PyCall, FITSIO, Images,PyPlot
@pyimport numpy as np
typealias Cube{T} Array{T,3}
typealias Mat{T} Array{T,2}
typealias Vec{T} Array{T,1}
typealias Vecx{T} LinSpace{T}

#=
 KH Yuen @ Lazarian Technology
 Last Modified : Jan 06 2018
 Note:
 (1) Numerical setting:
	(a) my cube has sonic speed of 0.19195751
	(b) my cube's mean field is always z
	(c) B field is scaled up by sqrt(2pi), i.e. real B = b in numerical cube/2pi
 (2) CL02 has the correct cos^2 factor
 and both CL03 and CL02 are actually using "mean field" for decomposition
 (3) The method is not applicable to Ma>1
 (4) The method is very memory intensive
=#

export getmsma,pcurl,pdot,kvectors_complex,kvectors_original,get_vsvfva_CL03,get_vsvfva_merged

function pcurl(ax,ay,az,bx,by,bz)
  cx = ay.*bz-az.*by;
  cy = az.*bx-ax.*bz;
  cz = ax.*by-ay.*bx;
  return cx,cy,cz
end

function pdot(ax,ay,az,bx,by,bz)
	bb=sqrt(bx^2+by^2+bz^2);
	ab=ax*bx+ay*by+az*bz
	cx=ab*bx/bb;
	cy=ab*by/bb;
	cz=ab*bz/bb;
	return cx,cy,cz
end

function getmsma(d::Cube,iv::Cube,jv::Cube,kv::Cube,ib::Cube,jb::Cube,kb::Cube);
	# For my cube
	cs=0.19195751
	# Ms Ma are RMS value
	ms=(var(iv)+var(jv)+var(kv))/cs
	ivx=iv-mean(iv);
	jvx=jv-mean(jv);
	kvx=kv-mean(kv);
	v2=ivx.^2+jvx.^2+kvx.^2;
	b2=ib.^2+jb.^2+kb.^2
	ma=sqrt(mean(d.*v2./b2));
	return ms,ma
end

function kvectors_original(d::Cube,iv::Cube,jv::Cube,kv::Cube,ib::Cube,jb::Cube,kb::Cube)
	# Mean field assumption = z holds
	# The expression in CL03 is different from CL02
	# We stick with CL03 for this function EXCEPT for the cos^2 factor
	ms,ma=getmsma(d,iv,jv,kv,ib,jb,kb);
	beta=2*(ma/ms)^2;
	kfx=zeros(size(d));kfy=zeros(size(d));kfz=zeros(size(d));
	ksx=zeros(size(d));ksy=zeros(size(d));ksz=zeros(size(d));
	kax=zeros(size(d));kay=zeros(size(d));kaz=zeros(size(d));

	for i in 1:nx,j in 1:ny, k in 1:nz
		idx=i-div(nx,2);
		jdx=j-div(ny,2);
		kdx=k-div(nz,2);
		dd=sqrt(idx.^2+jdx.^2+kdx.^2);
		ddx=sqrt(idx.^2+jdx.^2);
		alpha=beta/2;
		ctheta=kdx./dd;
		if (dd.>0)
			# from CL02
			D=(1+alpha)^2-4*alpha*ctheta*ctheta;
			ksx[i,j,k]=(1+alpha-sqrt(D))*idx;
			ksy[i,j,k]=(1+alpha-sqrt(D))*jdx;
			ksz[i,j,k]=(-1+alpha-sqrt(D))*kdx;
			kfx[i,j,k]=(1+alpha+sqrt(D))*idx;
			kfy[i,j,k]=(1+alpha+sqrt(D))*jdx;
			kfz[i,j,k]=(-1+alpha+sqrt(D))*kdx;
		end
		if(ddx.>0)
			kax[i,j,k]=jdx./ddx;
			kay[i,j,k]=-idx./ddx;
		end
	end
	return ksx,ksy,ksz,kfx,kfy,kfz,kax,kay,kaz
end

function kvectors_complex(d::Cube,iv::Cube,jv::Cube,kv::Cube,ib::Cube,jb::Cube,kb::Cube)
	# Mean field assumption = z does not hold
	# Using the CL02 expression instead
	ms,ma=getmsma(d,iv,jv,kv,ib,jb,kb);
	beta=2*(ma/ms)^2;
	kfx=zeros(size(d));kfy=zeros(size(d));kfz=zeros(size(d));
	ksx=zeros(size(d));ksy=zeros(size(d));ksz=zeros(size(d));
	kax=zeros(size(d));kay=zeros(size(d));kaz=zeros(size(d));
	
	# The mean field vector
	ibm=mean(ib);
	jbm=mean(jb);
	kbm=mean(kb);
	bm=sqrt(ibm.^2+jbm.^2+kbm.^2);

	for i in 1:nx,j in 1:ny, k in 1:nz
		idx=i-div(nx,2);
		jdx=j-div(ny,2);
		kdx=k-div(nz,2);
		dd=sqrt(idx.^2+jdx.^2+kdx.^2);
		
		# k parallel: Dot product with mean field
		alpha=beta/2;
		ctheta=(idx*ibm+jdx*jbm+kdx*kbm)/(dd*bm);
		
		kllx=ibm/bb*ctheta*dd;
		klly=jbm/bb*ctheta*dd;
		klly=kbm/bb*ctheta*dd;
		kll=sqrt(kllx^2+klly^2+kllz^2);
		
		# k prep = k - k_ll
		kppx=idx-kllx;
		kppy=jdx-klly;
		kppz=kdx-kllz;
		kpp=sqrt(kppx^2+kppy^2+kppz^2);

		if (dd.>0)
			# from CL02
			D=(1+alpha)^2-4*alpha*ctheta*ctheta;
			#= CL02:
			z_s \propto k_ll + (1-sqrt(D)-alpha)/(1+sqrt(D)+alpha)(k_ll/K_pp)^2 k_pp
			=#
			X=(1-sqrt(D)-alpha)/(1+sqrt(D)+alpha);
			kll_dot_k=(kllx*idx+klly*jdx+kllz*kdx)/dd;
			kpp_dot_k=(kppx*idx+kppy*jdx+kppz*kdx)/dd;
			ksx[i,j,k]=(kll_dot_k+X*(kll/kpp)^2*kpp_dot_k)*idx/dd;
			ksy[i,j,k]=(kll_dot_k+X*(kll/kpp)^2*kpp_dot_k)*jdx/dd;
			ksz[i,j,k]=(kll_dot_k+X*(kll/kpp)^2*kpp_dot_k)*kdx/dd;
			#= CL02:
			z_f \propto k_pp + (1-sqrt(D)+alpha)/(1+sqrt(D)-alpha)(k_pp/k_ll)^2 k_ll
			=#
			X=(1-sqrt(D)+alpha)/(1+sqrt(D)-alpha);
			kfx[i,j,k]=(kll_dot_k+X*(kpp/kll)^2*kpp_dot_k)*idx/dd;
			kfy[i,j,k]=(kll_dot_k+X*(kpp/kll)^2*kpp_dot_k)*jdx/dd;
			kfz[i,j,k]=(kll_dot_k+X*(kpp/kll)^2*kpp_dot_k)*kdx/dd;
			#= CL02:
			z_a=k_ll x k_pp
			=#
			kax[i,j,k],kay[i,j,k],kaz[i,j,k]=pcurl(kllx,klly,kllz,kppx,kppy,kppz)
		end
	end
	return ksx,ksy,ksz,kfx,kfy,kfz,kax,kay,kaz
end

function get_vsvfva_CL03(d::Cube,iv::Cube,jv::Cube,kv::Cube,ib::Cube,jb::Cube,kb::Cube);
	nx,ny,nz=size(d);
	iv=iv-mean(iv);
	jv=jv-mean(jv);
	kv=kv-mean(kv);
	
	# The fftshift will force coordinate system to follow the kvectors system (cube-center)
	ikv=fftshift(fft(iv));
	jkv=fftshift(fft(jv));
	kkv=fftshift(fft(kv));
	iv=0;jv=0;kv=0;gc();
	# Depends on your taste
	# Use any kvectors function you want above
	ksx,ksy,ksz,kfx,kfy,kfz,kax,kay,kaz=kvectors_original(d,iv,jv,kv,ib,jb,kb)
	
	kf=sqrt(kfx.^2+kfy.^2+kfz.^2);
	ks=sqrt(ksx.^2+ksy.^2+ksz.^2);
	ka=sqrt(kax.^2+kay.^2+kaz.^2);
	vkf=(ikv.*kfx+jkv.*kfy+kkv.*kfz)./kf;
	vkfx=vkf.*kfx./kf;
	vkfy=vkf.*kfy./kf;
	vkfz=vkf.*kfz./kf;
	vkfx[isnan(vkfx)]=0;
	vkfx[isnan(vkfy)]=0;
	vkfx[isnan(vkfz)]=0;
	kfx=0;kfy=0;kfz=0;kf=0;gc()
	vks=(ikv.*ksx+jkv.*ksy+kkv.*ksz)./ks;
	vksx=vks.*ksx./ks;
	vksy=vks.*ksy./ks;
	vksz=vks.*ksz./ks;
	vksx[isnan(vksx)]=0;
	vksy[isnan(vksy)]=0;
	vksz[isnan(vksz)]=0;
	ksx=0;ksy=0;ksz=0;ks=0;gc()
	vka=(ikv.*kax+jkv.*kay+kkv.*kaz)./ka;
	vkax=vka.*kax./ka;
	vkay=vka.*kay./ka;
	vkaz=vka.*kaz./ka;
	vkax[isnan(vkax)]=0;
	vkay[isnan(vkay)]=0;
	vkaz[isnan(vkaz)]=0;
	kax=0;kay=0;kaz=0;ka=0;gc()

	vfx=real(ifft(ifftshift(vkfx)));
	vsx=real(ifft(ifftshift(vksx)));
	vax=real(ifft(ifftshift(vkax)));
	vfy=real(ifft(ifftshift(vkfy)));
	vsy=real(ifft(ifftshift(vksy)));
	vay=real(ifft(ifftshift(vkay)));
	vfz=real(ifft(ifftshift(vkfz)));
	vsz=real(ifft(ifftshift(vksz)));
	vaz=real(ifft(ifftshift(vkaz)));
	return vsx,vsy,vsz,vfx,vfy,vfz,vax,vay,vaz
end

function get_vsvfva_merged(d::Cube,iv::Cube,jv::Cube,kv::Cube,ib::Cube,jb::Cube,kb::Cube)
	# Mean field assumption = z does not hold
	# Using the CL02 expression instead
	ms,ma=getmsma(d,iv,jv,kv,ib,jb,kb);
	beta=2*(ma/ms)^2;
	nx,ny,nz=size(d);
	iv=iv-mean(iv);
	jv=jv-mean(jv);
	kv=kv-mean(kv);
	
	# The fftshift will force coordinate system to follow the kvectors system (cube-center)
	ikv=fftshift(fft(iv));
	jkv=fftshift(fft(jv));
	kkv=fftshift(fft(kv));
	iv=0;jv=0;kv=0;gc();
	vfx=zeros(size(d));vfy=zeros(size(d));vfz=zeros(size(d));
	vsx=zeros(size(d));vsy=zeros(size(d));vsz=zeros(size(d));
	vax=zeros(size(d));vay=zeros(size(d));vaz=zeros(size(d));
	
	# The mean field vector
	ibm=mean(ib);
	jbm=mean(jb);
	kbm=mean(kb);
	bm=sqrt(ibm.^2+jbm.^2+kbm.^2);

	for i in 1:nx,j in 1:ny, k in 1:nz
		idx=i-div(nx,2);
		jdx=j-div(ny,2);
		kdx=k-div(nz,2);
		dd=sqrt(idx.^2+jdx.^2+kdx.^2);
		
		# k parallel: Dot product with mean field
		alpha=beta/2;
		ctheta=(idx*ibm+jdx*jbm+kdx*kbm)/(dd*bm);
		
		kllx=ibm/bb*ctheta;
		klly=jbm/bb*ctheta;
		klly=kbm/bb*ctheta;
		kll=sqrt(kllx^2+klly^2+kllz^2);
		
		# k prep = k - k_ll
		kppx=idx-kllx;
		kppy=jdx-klly;
		kppz=kdx-kllz;
		kpp=sqrt(kppx^2+kppy^2+kppz^2);

		if (dd.>0)
			# from CL02
			D=(1+alpha)^2-4*alpha*ctheta*ctheta;
			#= CL02:
				z_s \propto k_ll + (1-sqrt(D)-alpha)/(1+sqrt(D)+alpha)(k_ll/K_pp)^2 k_pp
			=#
			X=(1-sqrt(D)-alpha)/(1+sqrt(D)+alpha);
			kll_dot_k=(kllx*idx+klly*jdx+kllz*kdx)/dd;
			kpp_dot_k=(kppx*idx+kppy*jdx+kppz*kdx)/dd;
			ksx=(kll_dot_k+X*(kll/kpp)^2*kpp_dot_k)*idx/dd;
			ksy=(kll_dot_k+X*(kll/kpp)^2*kpp_dot_k)*jdx/dd;
			ksz=(kll_dot_k+X*(kll/kpp)^2*kpp_dot_k)*kdx/dd;
			ks=sqrt(ksx^2+ksy^2+ksz^2);
			#= CL02:
				z_f \propto k_pp + (1-sqrt(D)+alpha)/(1+sqrt(D)-alpha)(k_pp/k_ll)^2 k_ll
			=#
			X=(1-sqrt(D)+alpha)/(1+sqrt(D)-alpha);
			kfx=(kll_dot_k+X*(kpp/kll)^2*kpp_dot_k)*idx/dd;
			kfy=(kll_dot_k+X*(kpp/kll)^2*kpp_dot_k)*idx/dd;
			kfz=(kll_dot_k+X*(kpp/kll)^2*kpp_dot_k)*idx/dd;
			kf=sqrt(kfx^2+kfy^2+kfz^2);
			#= CL02:
				z_a=k_ll x k_pp
			=#
			kax,kay,kaz=pcurl(kllx,klly,kllz,kppx,kppy,kppz);
			ka=sqrt(kax^2+kay^2+kaz^2);
			
			# Handle the stupid dot product by function to increase readability
			vfx[i,j,k],vfy[i,j,k],vfz[i,j,k]=pdot(ikv[i,j,k],jvk[i,j,k],kkv[i,j,k],kfx,kfy,kfz);
			vsx[i,j,k],vsy[i,j,k],vsz[i,j,k]=pdot(ikv[i,j,k],jvk[i,j,k],kkv[i,j,k],ksx,ksy,ksz);
			vax[i,j,k],vay[i,j,k],vaz[i,j,k]=pdot(ikv[i,j,k],jvk[i,j,k],kkv[i,j,k],kax,kay,kaz);
		end
	end
	vfx=real(ifft(ifftshift(vfx)));
	vsx=real(ifft(ifftshift(vsx)));
	vax=real(ifft(ifftshift(vax)));
	vfy=real(ifft(ifftshift(vfy)));
	vsy=real(ifft(ifftshift(vsy)));
	vay=real(ifft(ifftshift(vay)));
	vfz=real(ifft(ifftshift(vfz)));
	vsz=real(ifft(ifftshift(vsz)));
	vaz=real(ifft(ifftshift(vaz)));
	return vsx,vsy,vsz,vfx,vfy,vfz,vax,vay,vaz
end

end # module CL03decompo

