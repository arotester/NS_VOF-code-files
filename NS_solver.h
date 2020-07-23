BC1();
normals();
HF();
usum=0.0;
BC2();
//_________________________________________________________________STAGE I___________________________________________________________________________//
#pragma omp parallel for schedule(dynamic) private(i,j,k,r_ue,r_ve,r_we,r_uw,r_vw,r_ww,r_un,r_vn,r_wn,r_us,r_vs,r_ws,r_ut,r_vt,r_wt,r_ub,r_vb,r_wb,a_ue,a_ve,a_we,a_uw,a_vw,a_ww,a_un,a_vn,a_wn,a_us,a_vs,a_ws,a_ut,a_vt,a_wt,a_ub,a_vb,a_wb,limiter_ue,limiter_ve,limiter_we,limiter_uw,limiter_vw,limiter_ww,limiter_un,limiter_vn,limiter_wn,limiter_us,limiter_vs,limiter_ws,limiter_ut,limiter_vt,limiter_wt,limiter_ub,limiter_vb,limiter_wb,F_surfx,F_surfy,F_surfz,rhoe,rhow,rhon,rhos,rhot,rhob,Ke,Kn,Kt,Kw,Ks,Kb,viscw,visce,viscs,viscn,visct,viscb,ue,uw,vn,vs,wt,wb,aw_nav,ae_nav,as_nav,an_nav,ab_nav,at_nav,ap_nav,uwa,uea,una,usa,uta,uba,vwa,vea,vna,vsa,vta,vba,wwa,wea,wna,wsa,wta,wba,fxad,fyad,fzad,fxdf,fydf,fzdf)
	for(i=img;i<img+nx;i++)
		{
		for(j=img;j<img+ny;j++)
			{
			for(k=img;k<img+nz;k++)
				{
                    Ke=(Kurv[i+1][j][k]+Kurv[i][j][k])*0.5;
					Kn=(Kurv[i][j+1][k]+Kurv[i][j][k])*0.5;
					Kt=(Kurv[i][j][k+1]+Kurv[i][j][k])*0.5;
					Kw=(Kurv[i-1][j][k]+Kurv[i][j][k])*0.5;
					Ks=(Kurv[i][j-1][k]+Kurv[i][j][k])*0.5;
					Kb=(Kurv[i][j][k-1]+Kurv[i][j][k])*0.5;

                    //Face viscosity
					viscw=2.0*(visc[i][j][k]*visc[i-1][j][k])/(visc[i][j][k]+visc[i-1][j][k]);
					visce=2.0*(visc[i][j][k]*visc[i+1][j][k])/(visc[i][j][k]+visc[i+1][j][k]);
					viscn=2.0*(visc[i][j][k]*visc[i][j+1][k])/(visc[i][j][k]+visc[i][j+1][k]);
					viscs=2.0*(visc[i][j][k]*visc[i][j-1][k])/(visc[i][j][k]+visc[i][j-1][k]);
					visct=2.0*(visc[i][j][k]*visc[i][j][k+1])/(visc[i][j][k]+visc[i][j][k+1]);
					viscb=2.0*(visc[i][j][k]*visc[i][j][k-1])/(visc[i][j][k]+visc[i][j][k-1]);
                    //Face Density
					rhoe=(rho[i+1][j][k]+rho[i][j][k])/2.0;
					rhow=(rho[i-1][j][k]+rho[i][j][k])/2.0;
					rhon=(rho[i][j+1][k]+rho[i][j][k])/2.0;
					rhos=(rho[i][j-1][k]+rho[i][j][k])/2.0;
					rhot=(rho[i][j][k]+rho[i][j][k+1])/2.0;
					rhob=(rho[i][j][k]+rho[i][j][k-1])/2.0;
                    //Advecting face velocity
					ue=(u[i+1][j][k]+u[i][j][k])/2.0;
					uw=(u[i-1][j][k]+u[i][j][k])/2.0;
					vn=(v[i][j+1][k]+v[i][j][k])/2.0;
					vs=(v[i][j-1][k]+v[i][j][k])/2.0;
					wt=(w[i][j][k+1]+w[i][j][k])/2.0;
					wb=(w[i][j][k-1]+w[i][j][k])/2.0;
                    //Coefficient Matrix
					aw_nav=dy*dz/(dx);
					ae_nav=aw_nav;
					as_nav=dx*dz/(dy);
					an_nav=as_nav;
					ab_nav=dx*dy/(dz);
					at_nav=ab_nav;
					ap_nav=-(ae_nav*visce+aw_nav*viscw+an_nav*viscn+as_nav*viscs+at_nav*visct+ab_nav*viscb);

                    //Advection fluxes
                    Superbee();  //old

					fxad=((dt/dx)*(uwa*uw-uea*ue)+(dt/dy)*(usa*vs-una*vn)+(dt/dz)*(uba*wb-uta*wt));	// X-Advection
					fyad=((dt/dx)*(vwa*uw-vea*ue)+(dt/dy)*(vsa*vs-vna*vn)+(dt/dz)*(vba*wb-vta*wt));	// Y-Advection
					fzad=((dt/dx)*(wwa*uw-wea*ue)+(dt/dy)*(wsa*vs-wna*vn)+(dt/dz)*(wba*wb-wta*wt));	// Z-Advection

					// X-Diffusion
					fxdf=(dt/(dx*dy*dz*Re*rho[i][j][k]))*(aw_nav*viscw*u[i-1][j][k]+ae_nav*visce*u[i+1][j][k]+as_nav*viscs*u[i][j-1][k]+an_nav*viscn*u[i][j+1][k]+ab_nav*viscb*u[i][j][k-1]+at_nav*visct*u[i][j][k+1]);
					// Y-Diffusion
					fydf=(dt/(dx*dy*dz*Re*rho[i][j][k]))*(aw_nav*viscw*v[i-1][j][k]+ae_nav*visce*v[i+1][j][k]+as_nav*viscs*v[i][j-1][k]+an_nav*viscn*v[i][j+1][k]+ab_nav*viscb*v[i][j][k-1]+at_nav*visct*v[i][j][k+1]);
					// Z-Diffusion
					fzdf=(dt/(dx*dy*dz*Re*rho[i][j][k]))*(aw_nav*viscw*w[i-1][j][k]+ae_nav*visce*w[i+1][j][k]+as_nav*viscs*w[i][j-1][k]+an_nav*viscn*w[i][j+1][k]+ab_nav*viscb*w[i][j][k-1]+at_nav*visct*w[i][j][k+1]);

					// X-Surface Tension
					F_surfx=dt*(1.0/We)*0.5*(Ke*((phi[i+1][j][k]-phi[i][j][k])/dx)/rhoe+Kw*((phi[i][j][k]-phi[i-1][j][k])/dx)/rhow);
					// Y-Surface Tension
					F_surfy=dt*(1.0/We)*0.5*(Kn*((phi[i][j+1][k]-phi[i][j][k])/dy)/rhon+Ks*((phi[i][j][k]-phi[i][j-1][k])/dy)/rhos);
					// Z-Surface Tension
					F_surfz=dt*(1.0/We)*0.5*(Kt*((phi[i][j][k+1]-phi[i][j][k])/dz)/rhot+Kb*((phi[i][j][k]-phi[i][j][k-1])/dz)/rhob);


                    ut[i][j][k]=u[i][j][k]*(1.0+ap_nav*dt/(Re*rho[i][j][k]*dx*dy*dz))+fxad+fxdf+F_surfx;
			 		vt[i][j][k]=v[i][j][k]*(1.0+ap_nav*dt/(Re*rho[i][j][k]*dx*dy*dz))+fyad+fydf-(dt/Froude)+F_surfy;
					wwt[i][j][k]=w[i][j][k]*(1.0+ap_nav*dt/(Re*rho[i][j][k]*dx*dy*dz))+fzad+fzdf+F_surfz;
				}
			}
		}

BC2t();
#pragma omp parallel for schedule(dynamic) private(i,j,k,utw,ute,vtn,vts,wtt,wtb,rhoe,rhow,rhos,rhon,rhot,rhob,uwp,uep,vsp,vnp,wtp,wbp)
	for(i=img;i<img+nx;i++)
		{
		for(j=img;j<img+ny;j++)
			{
			for(k=img;k<img+nz;k++)
				{
					utw=(ut[i][j][k]+ut[i-1][j][k])/2.0;
					ute=(ut[i+1][j][k]+ut[i][j][k])/2.0;
					vtn=(vt[i][j][k]+vt[i][j+1][k])/2.0;
					vts=(vt[i][j][k]+vt[i][j-1][k])/2.0;
					wtt=(wwt[i][j][k]+wwt[i][j][k+1])/2.0;
					wtb=(wwt[i][j][k]+wwt[i][j][k-1])/2.0;

					//face density
					rhoe=(rho[i+1][j][k]+rho[i][j][k])/2.0;
					rhow=(rho[i-1][j][k]+rho[i][j][k])/2.0;
					rhon=(rho[i][j+1][k]+rho[i][j][k])/2.0;
					rhos=(rho[i][j-1][k]+rho[i][j][k])/2.0;
					rhot=(rho[i][j][k]+rho[i][j][k+1])/2.0;
					rhob=(rho[i][j][k]+rho[i][j][k-1])/2.0;

					uwp=utw-(dt/(dx*rhow))*(pr[i][j][k]-pr[i-1][j][k]);
					uep=ute-(dt/(dx*rhoe))*(pr[i+1][j][k]-pr[i][j][k]);
					vsp=vts-(dt/(dy*rhos))*(pr[i][j][k]-pr[i][j-1][k]);
					vnp=vtn-(dt/(dy*rhon))*(pr[i][j+1][k]-pr[i][j][k]);
					wtp=wtt-(dt/(dz*rhot))*(pr[i][j][k+1]-pr[i][j][k]);
					wbp=wtb-(dt/(dz*rhob))*(pr[i][j][k]-pr[i][j][k-1]);

					so[i][j][k]=uwp*dy*dz-uep*dy*dz-vnp*dx*dz+vsp*dx*dz-wtp*dx*dy+wbp*dx*dy;

				}
			}
		}

// 		-------------------------------------Pressure Correction calculations (Bi-CGSTAB)--------------------------------------
BC_pc();
int gcount=0;
alp_num=0.0;
#pragma omp parallel for schedule(dynamic) private(i,j,k,rhoe,rhow,rhos,rhon,rhot,rhob,apg,aeg,asg,awg,ang,atg,abg) reduction(+:alp_num)
for(i=1;i<=nx;i++)
		{
		for(j=1;j<=ny;j++)
			{
            for(k=1;k<=nz;k++)
                {
                    frho();
                    apg=(dt*dy*dz/(rhoe*dx))+(dt*dy*dz/(rhow*dx))+(dx*dz*dt/(rhon*dy))+(dt*dx*dz/(rhos*dy))+(dx*dy*dt/(rhot*dz))+(dt*dx*dy/(rhob*dz));
                    aeg=-(dt*dy*dz/(rhoe*dx));
                    awg=-(dt*dy*dz/(rhow*dx));
                    ang=-(dx*dz*dt/(rhon*dy));
                    asg=-(dt*dx*dz/(rhos*dy));
                    atg=-(dx*dy*dt/(rhot*dz));
                    abg=-(dt*dx*dy/(rhob*dz));
                    Res[i][j][k]=s[i][j][k]-(aeg*(pc[i+1][j][k])+awg*(pc[i-1][j][k])+ang*(pc[i][j+1][k])+asg*(pc[i][j-1][k])+atg*(pc[i][j][k+1])+abg*(pc[i][j][k-1])+apg*pc[i][j][k]);
                    direc[i][j][k]=Res[i][j][k];
                    Res_star[i][j][k]=Res[i][j][k];
                    alp_num=alp_num+Res[i][j][k]*Res_star[i][j][k];
                }
            }
		}

do
{
BC_pc();
// Matrix multiplication [A]*[Res_0*] for alp_denominator
alp_den=0.0;
#pragma omp parallel for schedule(dynamic) private(i,j,k,rhoe,rhow,rhos,rhon,rhot,rhob,apg,aeg,asg,awg,ang,atg,abg) reduction(+:alp_den)
	for(i=1;i<=nx;i++)
		{
		for(j=1;j<=ny;j++)
			{
            for(k=1;k<=nz;k++)
                {

                    frho();
                    apg=(dt*dy*dz/(rhoe*dx))+(dt*dy*dz/(rhow*dx))+(dx*dz*dt/(rhon*dy))+(dt*dx*dz/(rhos*dy))+(dx*dy*dt/(rhot*dz))+(dt*dx*dy/(rhob*dz));
                    aeg=-(dt*dy*dz/(rhoe*dx));
                    awg=-(dt*dy*dz/(rhow*dx));
                    ang=-(dx*dz*dt/(rhon*dy));
                    asg=-(dt*dx*dz/(rhos*dy));
                    atg=-(dx*dy*dt/(rhot*dz));
                    abg=-(dt*dx*dy/(rhob*dz));
                    //Vertices
                    if(i==1 && j==1 && k==1)
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+ang*direc[i][j+1][k]+atg*direc[i][j][k+1];
                    if(i==nx && j==1 && k==1)
                        Ad[i][j][k]=apg*direc[i][j][k]+awg*direc[i-1][j][k]+ang*direc[i][j+1][k]+atg*direc[i][j][k+1];
                    if(i==nx && j==ny && k==nz)
                        Ad[i][j][k]=apg*direc[i][j][k]+awg*direc[i-1][j][k]+asg*direc[i][j-1][k]+abg*direc[i][j][k-1];
                    if(i==1 && j==ny && k==1)
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+asg*direc[i][j-1][k]+atg*direc[i][j][k+1];
                    if(i==1 && j==1 && k==nz)
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+ang*direc[i][j+1][k]+abg*direc[i][j][k-1];
                    if(i==nx && j==ny && k==1)
                        Ad[i][j][k]=apg*direc[i][j][k]+awg*direc[i-1][j][k]+asg*direc[i][j-1][k]+atg*direc[i][j][k+1];
                    if(i==nx && j==1 && k==nz)
                        Ad[i][j][k]=apg*direc[i][j][k]+awg*direc[i-1][j][k]+ang*direc[i][j+1][k]+abg*direc[i][j][k-1];
                    if(i==1 && j==ny && k==nz)
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+asg*direc[i][j-1][k]+abg*direc[i][j][k-1];
                    //Edges
                    if(i==1 && (j>1 && j<ny) && k==1)
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+ang*direc[i][j+1][k]+asg*direc[i][j-1][k]+atg*direc[i][j][k+1];
                    if(i==nx && (j>1 && j<ny) && k==1)
                        Ad[i][j][k]=apg*direc[i][j][k]+awg*direc[i-1][j][k]+ang*direc[i][j+1][k]+asg*direc[i][j-1][k]+atg*direc[i][j][k+1];
                    if(i==nx && (j>1 && j<ny) && k==nz)
                        Ad[i][j][k]=apg*direc[i][j][k]+awg*direc[i-1][j][k]+ang*direc[i][j+1][k]+asg*direc[i][j-1][k]+abg*direc[i][j][k-1];
                    if(i==1 && (j>1 && j<ny) && k==nz)
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+ang*direc[i][j+1][k]+asg*direc[i][j-1][k]+abg*direc[i][j][k-1];
                    if((i>1 && i<nx) && j==ny && k==1)
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+awg*direc[i-1][j][k]+asg*direc[i][j-1][k]+atg*direc[i][j][k+1];
                    if((i>1 && i<nx) && j==ny && k==nz)
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+awg*direc[i-1][j][k]+asg*direc[i][j-1][k]+abg*direc[i][j][k-1];
                    if((i>1 && i<nx) && j==1 && k==1)
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+awg*direc[i-1][j][k]+ang*direc[i][j+1][k]+atg*direc[i][j][k+1];
                    if((i>1 && i<nx) && j==1 && k==nz)
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+awg*direc[i-1][j][k]+ang*direc[i][j+1][k]+abg*direc[i][j][k-1];
                    if(i==nx && j==ny && (k>1 && k<nz))
                        Ad[i][j][k]=apg*direc[i][j][k]+awg*direc[i-1][j][k]+asg*direc[i][j-1][k]+atg*direc[i][j][k+1]+abg*direc[i][j][k-1];
                    if(i==1 && j==ny && (k>1 && k<nz))
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+asg*direc[i][j-1][k]+atg*direc[i][j][k+1]+abg*direc[i][j][k-1];
                    if(i==1 && j==1 && (k>1 && k<nz))
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+ang*direc[i][j+1][k]+atg*direc[i][j][k+1]+abg*direc[i][j][k-1];
                    if(i==nx && j==1 && (k>1 && k<nz))
                        Ad[i][j][k]=apg*direc[i][j][k]+awg*direc[i-1][j][k]+ang*direc[i][j+1][k]+atg*direc[i][j][k+1]+abg*direc[i][j][k-1];
                    // Faces
                    if((i>1 && i<nx) && j==1 && (k>1 && k<nz))
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+awg*direc[i-1][j][k]+ang*direc[i][j+1][k]+atg*direc[i][j][k+1]+abg*direc[i][j][k-1];
                    if((i>1 && i<nx) && j==ny && (k>1 && k<nz))
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+awg*direc[i-1][j][k]+asg*direc[i][j-1][k]+atg*direc[i][j][k+1]+abg*direc[i][j][k-1];
                    if(i==nx && (j>1 && j<ny) && (k>1 && k<nz))
                        Ad[i][j][k]=apg*direc[i][j][k]+awg*direc[i-1][j][k]+asg*direc[i][j-1][k]+ang*direc[i][j+1][k]+atg*direc[i][j][k+1]+abg*direc[i][j][k-1];
                    if(i==1 && (j>1 && j<ny) && (k>1 && k<nz))
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+asg*direc[i][j-1][k]+ang*direc[i][j+1][k]+atg*direc[i][j][k+1]+abg*direc[i][j][k-1];
                    if((i>1 && i<nx) && (j>1 && j<ny) && k==nz)
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+awg*direc[i-1][j][k]+asg*direc[i][j-1][k]+ang*direc[i][j+1][k]+abg*direc[i][j][k-1];
                    if((i>1 && i<nx) && (j>1 && j<ny) && k==1)
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+awg*direc[i-1][j][k]+asg*direc[i][j-1][k]+ang*direc[i][j+1][k]+atg*direc[i][j][k+1];
                    //Interior
                    if((1<i<nx) && (1<j<ny) && (k>1 && k<nz))
                        Ad[i][j][k]=apg*direc[i][j][k]+awg*direc[i-1][j][k]+aeg*direc[i+1][j][k]+asg*direc[i][j-1][k]+ang*direc[i][j+1][k]+atg*direc[i][j][k+1]+abg*direc[i][j][k-1];
                    alp_den=alp_den+Res_star[i][j][k]*Ad[i][j][k];
                }
            }
		}


// Finding alpha (Relaxation factor)
Relaxation=alp_num/alp_den;
// Shadow element 2 calculation
#pragma omp parallel for schedule(dynamic) private(i,j,k)
for(i=1;i<=nx;i++)
		{
		for(j=1;j<=ny;j++)
			{
            for(k=1;k<=nz;k++)
                {
                    shad2[i][j][k]=Res[i][j][k]-Relaxation*Ad[i][j][k];	// Updated pressure correction value
                }
            }
		}

// Matrix multiplication [A]*[s] for alp_denominator
om_num=0.0;
om_den=0.0;
#pragma omp parallel for schedule(dynamic) private(i,j,k,rhoe,rhow,rhos,rhon,rhot,rhob,apg,aeg,asg,awg,ang,atg,abg) reduction(+:om_num,om_den)
	for(i=1;i<=nx;i++)
		{
		for(j=1;j<=ny;j++)
			{
            for(k=1;k<=nz;k++)
                {
                    frho();
                     apg=(dt*dy*dz/(rhoe*dx))+(dt*dy*dz/(rhow*dx))+(dx*dz*dt/(rhon*dy))+(dt*dx*dz/(rhos*dy))+(dx*dy*dt/(rhot*dz))+(dt*dx*dy/(rhob*dz));
                    aeg=-(dt*dy*dz/(rhoe*dx));
                    awg=-(dt*dy*dz/(rhow*dx));
                    ang=-(dx*dz*dt/(rhon*dy));
                    asg=-(dt*dx*dz/(rhos*dy));
                    atg=-(dx*dy*dt/(rhot*dz));
                    abg=-(dt*dx*dy/(rhob*dz));

                    //Vertices
                    if(i==1 && j==1 && k==1)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+ang*shad2[i][j+1][k]+atg*shad2[i][j][k+1];
                    if(i==nx && j==1 && k==1)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+awg*shad2[i-1][j][k]+ang*shad2[i][j+1][k]+atg*shad2[i][j][k+1];
                    if(i==nx && j==ny && k==nz)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+awg*shad2[i-1][j][k]+asg*shad2[i][j-1][k]+abg*shad2[i][j][k-1];
                    if(i==1 && j==ny && k==1)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+asg*shad2[i][j-1][k]+atg*shad2[i][j][k+1];
                    if(i==1 && j==1 && k==nz)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+ang*shad2[i][j+1][k]+abg*shad2[i][j][k-1];
                    if(i==nx && j==ny && k==1)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+awg*shad2[i-1][j][k]+asg*shad2[i][j-1][k]+atg*shad2[i][j][k+1];
                    if(i==nx && j==1 && k==nz)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+awg*shad2[i-1][j][k]+ang*shad2[i][j+1][k]+abg*shad2[i][j][k-1];
                    if(i==1 && j==ny && k==nz)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+asg*shad2[i][j-1][k]+abg*shad2[i][j][k-1];
                    //Edges
                    if(i==1 && (j>1 && j<ny) && k==1)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+ang*shad2[i][j+1][k]+asg*shad2[i][j-1][k]+atg*shad2[i][j][k+1];
                    if(i==nx && (j>1 && j<ny) && k==1)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+awg*shad2[i-1][j][k]+ang*shad2[i][j+1][k]+asg*shad2[i][j-1][k]+atg*shad2[i][j][k+1];
                    if(i==nx && (j>1 && j<ny) && k==nz)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+awg*shad2[i-1][j][k]+ang*shad2[i][j+1][k]+asg*shad2[i][j-1][k]+abg*shad2[i][j][k-1];
                    if(i==1 && (j>1 && j<ny) && k==nz)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+ang*shad2[i][j+1][k]+asg*shad2[i][j-1][k]+abg*shad2[i][j][k-1];
                    if((i>1 && i<nx) && j==ny && k==1)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+awg*shad2[i-1][j][k]+asg*shad2[i][j-1][k]+atg*shad2[i][j][k+1];
                    if((i>1 && i<nx) && j==ny && k==nz)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+awg*shad2[i-1][j][k]+asg*shad2[i][j-1][k]+abg*shad2[i][j][k-1];
                    if((i>1 && i<nx) && j==1 && k==1)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+awg*shad2[i-1][j][k]+ang*shad2[i][j+1][k]+atg*shad2[i][j][k+1];
                    if((i>1 && i<nx) && j==1 && k==nz)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+awg*shad2[i-1][j][k]+ang*shad2[i][j+1][k]+abg*shad2[i][j][k-1];
                    if(i==nx && j==ny && (k>1 && k<nz))
                        Ashad[i][j][k]=apg*shad2[i][j][k]+awg*shad2[i-1][j][k]+asg*shad2[i][j-1][k]+atg*shad2[i][j][k+1]+abg*shad2[i][j][k-1];
                    if(i==1 && j==ny && (k>1 && k<nz))
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+asg*shad2[i][j-1][k]+atg*shad2[i][j][k+1]+abg*shad2[i][j][k-1];
                    if(i==1 && j==1 && (k>1 && k<nz))
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+ang*shad2[i][j+1][k]+atg*shad2[i][j][k+1]+abg*shad2[i][j][k-1];
                    if(i==nx && j==1 && (k>1 && k<nz))
                        Ashad[i][j][k]=apg*shad2[i][j][k]+awg*shad2[i-1][j][k]+ang*shad2[i][j+1][k]+atg*shad2[i][j][k+1]+abg*shad2[i][j][k-1];
                    // Faces
                    if((i>1 && i<nx) && j==1 && (k>1 && k<nz))
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+awg*shad2[i-1][j][k]+ang*shad2[i][j+1][k]+atg*shad2[i][j][k+1]+abg*shad2[i][j][k-1];
                    if((i>1 && i<nx) && j==ny && (k>1 && k<nz))
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+awg*shad2[i-1][j][k]+asg*shad2[i][j-1][k]+atg*shad2[i][j][k+1]+abg*shad2[i][j][k-1];
                    if(i==nx && (j>1 && j<ny) && (k>1 && k<nz))
                        Ashad[i][j][k]=apg*shad2[i][j][k]+awg*shad2[i-1][j][k]+asg*shad2[i][j-1][k]+ang*shad2[i][j+1][k]+atg*shad2[i][j][k+1]+abg*shad2[i][j][k-1];
                    if(i==1 && (j>1 && j<ny) && (k>1 && k<nz))
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+asg*shad2[i][j-1][k]+ang*shad2[i][j+1][k]+atg*shad2[i][j][k+1]+abg*shad2[i][j][k-1];
                    if((i>1 && i<nx) && (j>1 && j<ny) && k==nz)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+awg*shad2[i-1][j][k]+asg*shad2[i][j-1][k]+ang*shad2[i][j+1][k]+abg*shad2[i][j][k-1];
                    if((i>1 && i<nx) && (j>1 && j<ny) && k==1)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+awg*shad2[i-1][j][k]+asg*shad2[i][j-1][k]+ang*shad2[i][j+1][k]+atg*shad2[i][j][k+1];
                    //Interior
                    if((1<i<nx) && (1<j<ny) && (k>1 && k<nz))
                        Ashad[i][j][k]=apg*shad2[i][j][k]+awg*shad2[i-1][j][k]+aeg*shad2[i+1][j][k]+asg*shad2[i][j-1][k]+ang*shad2[i][j+1][k]+atg*shad2[i][j][k+1]+abg*shad2[i][j][k-1];
                    om_num=om_num+Ashad[i][j][k]*shad2[i][j][k];
                    om_den=om_den+Ashad[i][j][k]*Ashad[i][j][k];
                }
            }
		}
//Calculating omega
omega=om_num/om_den;
#pragma omp parallel for schedule(dynamic) private(i,j,k)
for(i=1;i<=nx;i++)
		{
		for(j=1;j<=ny;j++)
			{
            for(k=1;k<=nz;k++)
                {
                    pc[i][j][k]=pc[i][j][k]+Relaxation*direc[i][j][k]+omega*(shad2[i][j][k]);	// Updated pressure correction value
                    Res[i][j][k]=shad2[i][j][k]-omega*Ashad[i][j][k];	// New Residue
                }
			}
		}

b_den=alp_num;
alp_num=0.0;
#pragma omp parallel for schedule(dynamic) private(i,j,k) reduction(+:alp_num)
for(i=1;i<=nx;i++)
		{
		for(j=1;j<=ny;j++)
			{
            for(k=1;k<=nz;k++)
                {
                    alp_num=alp_num+Res[i][j][k]*Res_star[i][j][k];		// updated value for next iteration
                }
			}
		}

direc_fac=(alp_num/b_den)*(Relaxation/omega);
#pragma omp parallel for schedule(dynamic) private(i,j,k)
for(i=1;i<=nx;i++)
		{
		for(j=1;j<=ny;j++)
			{
            for(k=1;k<=nz;k++)
                {
                    direc[i][j][k]=Res[i][j][k]+direc_fac*(direc[i][j][k]-omega*Ad[i][j][k]);	// updated value for the next iteration
                }
			}
		}
gRMS=sqrt((alp_num/(nx*ny*nz)));
gcount++;

}while(gRMS>=0.000001);

#pragma omp parallel for schedule(dynamic) private(i,j,k)
for(i=img;i<img+nx;i++)
		{
		for(j=img;j<img+ny;j++)
			{
			for(k=img;k<img+nz;k++)
				{
					pnew[i][j][k]=pr[i][j][k]+pc[i][j][k];
					pc[i][j][k]=0.0;
				}
			}
		}

BC_cor();
#pragma omp parallel for schedule(dynamic) private(i,j,k,diff,sqr,rhot,rhob,rhoe,rhow,rhos,rhon) reduction(+:usum)
	for(i=img;i<img+nx;i++)
		{
		for(j=img;j<img+ny;j++)
			{
			for(k=img;k<img+nz;k++)
				{

					//Face density
					rhoe=(rho[i+1][j][k]+rho[i][j][k])/2.0;
					rhow=(rho[i-1][j][k]+rho[i][j][k])/2.0;
					rhon=(rho[i][j+1][k]+rho[i][j][k])/2.0;
					rhos=(rho[i][j-1][k]+rho[i][j][k])/2.0;
					rhot=(rho[i][j][k]+rho[i][j][k+1])/2.0;
					rhob=(rho[i][j][k]+rho[i][j][k-1])/2.0;

					uc[i][j][k]=0.5*(dt/dx)*((pnew[i+1][j][k]/rhoe)-(pnew[i-1][j][k]/rhow)+pnew[i][j][k]*(1.0/rhow-1.0/rhoe));
					vc[i][j][k]=0.5*(dt/dy)*((pnew[i][j+1][k]/rhon)-(pnew[i][j-1][k]/rhos)+pnew[i][j][k]*(1.0/rhos-1.0/rhon));
					wc[i][j][k]=0.5*(dt/dz)*((pnew[i][j][k+1]/rhot)-(pnew[i][j][k-1]/rhob)+pnew[i][j][k]*(1.0/rhob-1.0/rhot));
					unew[i][j][k]=ut[i][j][k]-uc[i][j][k];
					vnew[i][j][k]=vt[i][j][k]-vc[i][j][k];
					wnew[i][j][k]=wwt[i][j][k]-wc[i][j][k];
					//RMS for steady state
					diff=unew[i][j][k]-u[i][j][k];
					sqr=diff*diff;
					usum=usum+sqr;

					//Updating the old values
					u1[i][j][k]=unew[i][j][k];
					v1[i][j][k]=vnew[i][j][k];
					w1[i][j][k]=wnew[i][j][k];
					pr1[i][j][k]=pnew[i][j][k];
				//	printf("i=%d  j=%d  uc=%f\tvc=%f\tp=%f\n",i,j,uc[i][j],vc[i][j],pnew[i][j]);
				}
			}
		}

//_________________________________________________________________STAGE II___________________________________________________________________________//
BC2();
BC_in();
#pragma omp parallel for schedule(dynamic) private(i,j,k,r_ue,r_ve,r_we,r_uw,r_vw,r_ww,r_un,r_vn,r_wn,r_us,r_vs,r_ws,r_ut,r_vt,r_wt,r_ub,r_vb,r_wb,a_ue,a_ve,a_we,a_uw,a_vw,a_ww,a_un,a_vn,a_wn,a_us,a_vs,a_ws,a_ut,a_vt,a_wt,a_ub,a_vb,a_wb,limiter_ue,limiter_ve,limiter_we,limiter_uw,limiter_vw,limiter_ww,limiter_un,limiter_vn,limiter_wn,limiter_us,limiter_vs,limiter_ws,limiter_ut,limiter_vt,limiter_wt,limiter_ub,limiter_vb,limiter_wb,F_surfx,F_surfy,F_surfz,rhoe,rhow,rhon,rhos,rhot,rhob,Ke,Kn,Kt,Kw,Ks,Kb,viscw,visce,viscs,viscn,visct,viscb,ue,uw,vn,vs,wt,wb,aw_nav,ae_nav,as_nav,an_nav,ab_nav,at_nav,ap_nav,uwa,uea,una,usa,uta,uba,vwa,vea,vna,vsa,vta,vba,wwa,wea,wna,wsa,wta,wba,fxad,fyad,fzad,fxdf,fydf,fzdf,uwa1,uea1,una1,usa1,uta1,uba1,vwa1,vea1,vna1,vsa1,vta1,vba1,wwa1,wea1,wna1,wsa1,wta1,wba1,fxad1,fyad1,fzad1,fxdf1,fydf1,fzdf1,ue1,uw1,vs1,,vn1,wb1,wt1)
	for(i=img;i<img+nx;i++)
		{
		for(j=img;j<img+ny;j++)
			{
			for(k=img;k<img+nz;k++)
				{
                    Ke=(Kurv[i+1][j][k]+Kurv[i][j][k])*0.5;
					Kn=(Kurv[i][j+1][k]+Kurv[i][j][k])*0.5;
					Kt=(Kurv[i][j][k+1]+Kurv[i][j][k])*0.5;
					Kw=(Kurv[i-1][j][k]+Kurv[i][j][k])*0.5;
					Ks=(Kurv[i][j-1][k]+Kurv[i][j][k])*0.5;
					Kb=(Kurv[i][j][k-1]+Kurv[i][j][k])*0.5;

                    //Face viscosity
					viscw=2.0*(visc[i][j][k]*visc[i-1][j][k])/(visc[i][j][k]+visc[i-1][j][k]);
					visce=2.0*(visc[i][j][k]*visc[i+1][j][k])/(visc[i][j][k]+visc[i+1][j][k]);
					viscn=2.0*(visc[i][j][k]*visc[i][j+1][k])/(visc[i][j][k]+visc[i][j+1][k]);
					viscs=2.0*(visc[i][j][k]*visc[i][j-1][k])/(visc[i][j][k]+visc[i][j-1][k]);
					visct=2.0*(visc[i][j][k]*visc[i][j][k+1])/(visc[i][j][k]+visc[i][j][k+1]);
					viscb=2.0*(visc[i][j][k]*visc[i][j][k-1])/(visc[i][j][k]+visc[i][j][k-1]);
                    //Face Density
					rhoe=(rho[i+1][j][k]+rho[i][j][k])/2.0;
					rhow=(rho[i-1][j][k]+rho[i][j][k])/2.0;
					rhon=(rho[i][j+1][k]+rho[i][j][k])/2.0;
					rhos=(rho[i][j-1][k]+rho[i][j][k])/2.0;
					rhot=(rho[i][j][k]+rho[i][j][k+1])/2.0;
					rhob=(rho[i][j][k]+rho[i][j][k-1])/2.0;
                    //Advecting face velocity
					ue=(u[i+1][j][k]+u[i][j][k])/2.0;
					uw=(u[i-1][j][k]+u[i][j][k])/2.0;
					vn=(v[i][j+1][k]+v[i][j][k])/2.0;
					vs=(v[i][j-1][k]+v[i][j][k])/2.0;
					wt=(w[i][j][k+1]+w[i][j][k])/2.0;
					wb=(w[i][j][k-1]+w[i][j][k])/2.0;
                    //Coefficient Matrix
					aw_nav=dy*dz/(dx);
					ae_nav=aw_nav;
					as_nav=dx*dz/(dy);
					an_nav=as_nav;
					ab_nav=dx*dy/(dz);
					at_nav=ab_nav;
					ap_nav=-(ae_nav*visce+aw_nav*viscw+an_nav*viscn+as_nav*viscs+at_nav*visct+ab_nav*viscb);

                    //Advection fluxes
                    Superbee();  //old
                    Superbee1(); //Inter
//_______________________________old__________________________________________________//
					fxad=((dt/dx)*(uwa*uw-uea*ue)+(dt/dy)*(usa*vs-una*vn)+(dt/dz)*(uba*wb-uta*wt));	// X-Advection
					fyad=((dt/dx)*(vwa*uw-vea*ue)+(dt/dy)*(vsa*vs-vna*vn)+(dt/dz)*(vba*wb-vta*wt));	// Y-Advection
					fzad=((dt/dx)*(wwa*uw-wea*ue)+(dt/dy)*(wsa*vs-wna*vn)+(dt/dz)*(wba*wb-wta*wt));	// Z-Advection

					// X-Diffusion
					fxdf=(dt/(dx*dy*dz*Re*rho[i][j][k]))*(aw_nav*viscw*u[i-1][j][k]+ae_nav*visce*u[i+1][j][k]+as_nav*viscs*u[i][j-1][k]+an_nav*viscn*u[i][j+1][k]+ab_nav*viscb*u[i][j][k-1]+at_nav*visct*u[i][j][k+1]);
					// Y-Diffusion
					fydf=(dt/(dx*dy*dz*Re*rho[i][j][k]))*(aw_nav*viscw*v[i-1][j][k]+ae_nav*visce*v[i+1][j][k]+as_nav*viscs*v[i][j-1][k]+an_nav*viscn*v[i][j+1][k]+ab_nav*viscb*v[i][j][k-1]+at_nav*visct*v[i][j][k+1]);
					// Z-Diffusion
					fzdf=(dt/(dx*dy*dz*Re*rho[i][j][k]))*(aw_nav*viscw*w[i-1][j][k]+ae_nav*visce*w[i+1][j][k]+as_nav*viscs*w[i][j-1][k]+an_nav*viscn*w[i][j+1][k]+ab_nav*viscb*w[i][j][k-1]+at_nav*visct*w[i][j][k+1]);

//_______________________________Inter__________________________________________________//
					fxad1=((dt/dx)*(uwa1*uw1-uea1*ue1)+(dt/dy)*(usa1*vs1-una1*vn1)+(dt/dz)*(uba1*wb1-uta1*wt1));	// X-Advection
					fyad1=((dt/dx)*(vwa1*uw1-vea1*ue1)+(dt/dy)*(vsa1*vs1-vna1*vn1)+(dt/dz)*(vba1*wb1-vta1*wt1));	// Y-Advection
					fzad1=((dt/dx)*(wwa1*uw1-wea1*ue1)+(dt/dy)*(wsa1*vs1-wna1*vn1)+(dt/dz)*(wba1*wb1-wta1*wt1));	// Z-Advection

					// X-Diffusion
					fxdf1=(dt/(dx*dy*dz*Re*rho[i][j][k]))*(aw_nav*viscw*u1[i-1][j][k]+ae_nav*visce*u1[i+1][j][k]+as_nav*viscs*u1[i][j-1][k]+an_nav*viscn*u1[i][j+1][k]+ab_nav*viscb*u1[i][j][k-1]+at_nav*visct*u1[i][j][k+1]);
					// Y-Diffusion
					fydf1=(dt/(dx*dy*dz*Re*rho[i][j][k]))*(aw_nav*viscw*v1[i-1][j][k]+ae_nav*visce*v1[i+1][j][k]+as_nav*viscs*v1[i][j-1][k]+an_nav*viscn*v1[i][j+1][k]+ab_nav*viscb*v1[i][j][k-1]+at_nav*visct*v1[i][j][k+1]);
					// Z-Diffusion
					fzdf1=(dt/(dx*dy*dz*Re*rho[i][j][k]))*(aw_nav*viscw*w1[i-1][j][k]+ae_nav*visce*w1[i+1][j][k]+as_nav*viscs*w1[i][j-1][k]+an_nav*viscn*w1[i][j+1][k]+ab_nav*viscb*w1[i][j][k-1]+at_nav*visct*w1[i][j][k+1]);

					// X-Surface Tension
					F_surfx=dt*(1.0/We)*0.5*(Ke*((phi[i+1][j][k]-phi[i][j][k])/dx)/rhoe+Kw*((phi[i][j][k]-phi[i-1][j][k])/dx)/rhow);
					// Y-Surface Tension
					F_surfy=dt*(1.0/We)*0.5*(Kn*((phi[i][j+1][k]-phi[i][j][k])/dy)/rhon+Ks*((phi[i][j][k]-phi[i][j-1][k])/dy)/rhos);
					// Z-Surface Tension
					F_surfz=dt*(1.0/We)*0.5*(Kt*((phi[i][j][k+1]-phi[i][j][k])/dz)/rhot+Kb*((phi[i][j][k]-phi[i][j][k-1])/dz)/rhob);

                    ut[i][j][k]=u[i][j][k]*(1.0+ap_nav*dt/(Re*rho[i][j][k]*dx*dy*dz))+0.5*(fxad+fxad1)+0.5*(fxdf+fxdf1)+F_surfx;
			 		vt[i][j][k]=v[i][j][k]*(1.0+ap_nav*dt/(Re*rho[i][j][k]*dx*dy*dz))+0.5*(fyad+fyad1)+0.5*(fydf+fydf1)-(dt/Froude)+F_surfy;
					wwt[i][j][k]=w[i][j][k]*(1.0+ap_nav*dt/(Re*rho[i][j][k]*dx*dy*dz))+0.5*(fzad+fzad1)+0.5*(fzdf+fzdf1)+F_surfz;
				}
			}
		}
BC_t_in();
#pragma omp parallel for schedule(dynamic) private(i,j,k,utw,ute,vtn,vts,wtt,wtb,rhoe,rhow,rhos,rhon,rhot,rhob,uwp,uep,vsp,vnp,wtp,wbp)
	for(i=img;i<img+nx;i++)
		{
		for(j=img;j<img+ny;j++)
			{
			for(k=img;k<img+nz;k++)
				{
					utw=(ut[i][j][k]+ut[i-1][j][k])/2.0;
					ute=(ut[i+1][j][k]+ut[i][j][k])/2.0;
					vtn=(vt[i][j][k]+vt[i][j+1][k])/2.0;
					vts=(vt[i][j][k]+vt[i][j-1][k])/2.0;
					wtt=(wwt[i][j][k]+wwt[i][j][k+1])/2.0;
					wtb=(wwt[i][j][k]+wwt[i][j][k-1])/2.0;

					//face density
					rhoe=(rho[i+1][j][k]+rho[i][j][k])/2.0;
					rhow=(rho[i-1][j][k]+rho[i][j][k])/2.0;
					rhon=(rho[i][j+1][k]+rho[i][j][k])/2.0;
					rhos=(rho[i][j-1][k]+rho[i][j][k])/2.0;
					rhot=(rho[i][j][k]+rho[i][j][k+1])/2.0;
					rhob=(rho[i][j][k]+rho[i][j][k-1])/2.0;

					uwp=utw-(dt/(dx*rhow))*(0.5*(pr[i][j][k]-pr[i-1][j][k])+0.5*(pr1[i][j][k]-pr1[i-1][j][k]));
					uep=ute-(dt/(dx*rhoe))*(0.5*(pr[i+1][j][k]-pr[i][j][k])+0.5*(pr1[i+1][j][k]-pr1[i][j][k]));
					vsp=vts-(dt/(dy*rhos))*(0.5*(pr[i][j][k]-pr[i][j-1][k])+0.5*(pr1[i][j][k]-pr1[i][j-1][k]));
					vnp=vtn-(dt/(dy*rhon))*(0.5*(pr[i][j+1][k]-pr[i][j][k])+0.5*(pr1[i][j+1][k]-pr1[i][j][k]));
					wtp=wtt-(dt/(dz*rhot))*(0.5*(pr[i][j][k+1]-pr[i][j][k])+0.5*(pr1[i][j][k+1]-pr1[i][j][k]));
					wbp=wtb-(dt/(dz*rhob))*(0.5*(pr[i][j][k]-pr[i][j][k-1])+0.5*(pr1[i][j][k]-pr1[i][j][k-1]));

					so[i][j][k]=uwp*dy*dz-uep*dy*dz-vnp*dx*dz+vsp*dx*dz-wtp*dx*dy+wbp*dx*dy;

			//		printf("source=%f %f %f %f\n",utw,ute,vtn,vts);
				}
			}
		}

// 		-------------------------------------Pressure Correction calculations (Bi-CGSTAB)--------------------------------------
gcount=0;
for(i=img;i<img+nx;i++)
		{
		for(j=img;j<img+ny;j++)
			{
			for(k=img;k<img+nz;k++)
				{
                    pc[i][j][k]=0.0;
				}
			}
		}

BC_pc()
alp_num=0.0;
#pragma omp parallel for schedule(dynamic) private(i,j,k,rhoe,rhow,rhos,rhon,rhot,rhob,apg,aeg,asg,awg,ang,atg,abg) reduction(+:alp_num)
for(i=1;i<=nx;i++)
		{
		for(j=1;j<=ny;j++)
			{
            for(k=1;k<=nz;k++)
                {
                    frho();
                    apg=(dt*dy*dz/(rhoe*dx))+(dt*dy*dz/(rhow*dx))+(dx*dz*dt/(rhon*dy))+(dt*dx*dz/(rhos*dy))+(dx*dy*dt/(rhot*dz))+(dt*dx*dy/(rhob*dz));
                    aeg=-(dt*dy*dz/(rhoe*dx));
                    awg=-(dt*dy*dz/(rhow*dx));
                    ang=-(dx*dz*dt/(rhon*dy));
                    asg=-(dt*dx*dz/(rhos*dy));
                    atg=-(dx*dy*dt/(rhot*dz));
                    abg=-(dt*dx*dy/(rhob*dz));
                    Res[i][j][k]=s[i][j][k]-(aeg*(pc[i+1][j][k])+awg*(pc[i-1][j][k])+ang*(pc[i][j+1][k])+asg*(pc[i][j-1][k])+atg*(pc[i][j][k+1])+abg*(pc[i][j][k-1])+apg*pc[i][j][k]);
                    direc[i][j][k]=Res[i][j][k];
                    Res_star[i][j][k]=Res[i][j][k];
                    alp_num=alp_num+Res[i][j][k]*Res_star[i][j][k];
                }
            }
		}

do
{
BC_pc();

// Matrix multiplication [A]*[Res_0*] for alp_denominator
alp_den=0.0;
#pragma omp parallel for schedule(dynamic) private(i,j,k,rhoe,rhow,rhos,rhon,rhot,rhob,apg,aeg,asg,awg,ang,atg,abg) reduction(+:alp_den)
	for(i=1;i<=nx;i++)
		{
		for(j=1;j<=ny;j++)
			{
            for(k=1;k<=nz;k++)
                {

                    frho();
                    apg=(dt*dy*dz/(rhoe*dx))+(dt*dy*dz/(rhow*dx))+(dx*dz*dt/(rhon*dy))+(dt*dx*dz/(rhos*dy))+(dx*dy*dt/(rhot*dz))+(dt*dx*dy/(rhob*dz));
                    aeg=-(dt*dy*dz/(rhoe*dx));
                    awg=-(dt*dy*dz/(rhow*dx));
                    ang=-(dx*dz*dt/(rhon*dy));
                    asg=-(dt*dx*dz/(rhos*dy));
                    atg=-(dx*dy*dt/(rhot*dz));
                    abg=-(dt*dx*dy/(rhob*dz));
                    //Vertices
                    if(i==1 && j==1 && k==1)
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+ang*direc[i][j+1][k]+atg*direc[i][j][k+1];
                    if(i==nx && j==1 && k==1)
                        Ad[i][j][k]=apg*direc[i][j][k]+awg*direc[i-1][j][k]+ang*direc[i][j+1][k]+atg*direc[i][j][k+1];
                    if(i==nx && j==ny && k==nz)
                        Ad[i][j][k]=apg*direc[i][j][k]+awg*direc[i-1][j][k]+asg*direc[i][j-1][k]+abg*direc[i][j][k-1];
                    if(i==1 && j==ny && k==1)
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+asg*direc[i][j-1][k]+atg*direc[i][j][k+1];
                    if(i==1 && j==1 && k==nz)
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+ang*direc[i][j+1][k]+abg*direc[i][j][k-1];
                    if(i==nx && j==ny && k==1)
                        Ad[i][j][k]=apg*direc[i][j][k]+awg*direc[i-1][j][k]+asg*direc[i][j-1][k]+atg*direc[i][j][k+1];
                    if(i==nx && j==1 && k==nz)
                        Ad[i][j][k]=apg*direc[i][j][k]+awg*direc[i-1][j][k]+ang*direc[i][j+1][k]+abg*direc[i][j][k-1];
                    if(i==1 && j==ny && k==nz)
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+asg*direc[i][j-1][k]+abg*direc[i][j][k-1];
                    //Edges
                    if(i==1 && (j>1 && j<ny) && k==1)
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+ang*direc[i][j+1][k]+asg*direc[i][j-1][k]+atg*direc[i][j][k+1];
                    if(i==nx && (j>1 && j<ny) && k==1)
                        Ad[i][j][k]=apg*direc[i][j][k]+awg*direc[i-1][j][k]+ang*direc[i][j+1][k]+asg*direc[i][j-1][k]+atg*direc[i][j][k+1];
                    if(i==nx && (j>1 && j<ny) && k==nz)
                        Ad[i][j][k]=apg*direc[i][j][k]+awg*direc[i-1][j][k]+ang*direc[i][j+1][k]+asg*direc[i][j-1][k]+abg*direc[i][j][k-1];
                    if(i==1 && (j>1 && j<ny) && k==nz)
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+ang*direc[i][j+1][k]+asg*direc[i][j-1][k]+abg*direc[i][j][k-1];
                    if((i>1 && i<nx) && j==ny && k==1)
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+awg*direc[i-1][j][k]+asg*direc[i][j-1][k]+atg*direc[i][j][k+1];
                    if((i>1 && i<nx) && j==ny && k==nz)
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+awg*direc[i-1][j][k]+asg*direc[i][j-1][k]+abg*direc[i][j][k-1];
                    if((i>1 && i<nx) && j==1 && k==1)
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+awg*direc[i-1][j][k]+ang*direc[i][j+1][k]+atg*direc[i][j][k+1];
                    if((i>1 && i<nx) && j==1 && k==nz)
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+awg*direc[i-1][j][k]+ang*direc[i][j+1][k]+abg*direc[i][j][k-1];
                    if(i==nx && j==ny && (k>1 && k<nz))
                        Ad[i][j][k]=apg*direc[i][j][k]+awg*direc[i-1][j][k]+asg*direc[i][j-1][k]+atg*direc[i][j][k+1]+abg*direc[i][j][k-1];
                    if(i==1 && j==ny && (k>1 && k<nz))
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+asg*direc[i][j-1][k]+atg*direc[i][j][k+1]+abg*direc[i][j][k-1];
                    if(i==1 && j==1 && (k>1 && k<nz))
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+ang*direc[i][j+1][k]+atg*direc[i][j][k+1]+abg*direc[i][j][k-1];
                    if(i==nx && j==1 && (k>1 && k<nz))
                        Ad[i][j][k]=apg*direc[i][j][k]+awg*direc[i-1][j][k]+ang*direc[i][j+1][k]+atg*direc[i][j][k+1]+abg*direc[i][j][k-1];
                    // Faces
                    if((i>1 && i<nx) && j==1 && (k>1 && k<nz))
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+awg*direc[i-1][j][k]+ang*direc[i][j+1][k]+atg*direc[i][j][k+1]+abg*direc[i][j][k-1];
                    if((i>1 && i<nx) && j==ny && (k>1 && k<nz))
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+awg*direc[i-1][j][k]+asg*direc[i][j-1][k]+atg*direc[i][j][k+1]+abg*direc[i][j][k-1];
                    if(i==nx && (j>1 && j<ny) && (k>1 && k<nz))
                        Ad[i][j][k]=apg*direc[i][j][k]+awg*direc[i-1][j][k]+asg*direc[i][j-1][k]+ang*direc[i][j+1][k]+atg*direc[i][j][k+1]+abg*direc[i][j][k-1];
                    if(i==1 && (j>1 && j<ny) && (k>1 && k<nz))
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+asg*direc[i][j-1][k]+ang*direc[i][j+1][k]+atg*direc[i][j][k+1]+abg*direc[i][j][k-1];
                    if((i>1 && i<nx) && (j>1 && j<ny) && k==nz)
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+awg*direc[i-1][j][k]+asg*direc[i][j-1][k]+ang*direc[i][j+1][k]+abg*direc[i][j][k-1];
                    if((i>1 && i<nx) && (j>1 && j<ny) && k==1)
                        Ad[i][j][k]=apg*direc[i][j][k]+aeg*direc[i+1][j][k]+awg*direc[i-1][j][k]+asg*direc[i][j-1][k]+ang*direc[i][j+1][k]+atg*direc[i][j][k+1];
                    //Interior
                    if((1<i<nx) && (1<j<ny) && (k>1 && k<nz))
                        Ad[i][j][k]=apg*direc[i][j][k]+awg*direc[i-1][j][k]+aeg*direc[i+1][j][k]+asg*direc[i][j-1][k]+ang*direc[i][j+1][k]+atg*direc[i][j][k+1]+abg*direc[i][j][k-1];
                    alp_den=alp_den+Res_star[i][j][k]*Ad[i][j][k];
                }
            }
		}


// Finding alpha (Relaxation factor)
Relaxation=alp_num/alp_den;
// Shadow element 2 calculation
#pragma omp parallel for schedule(dynamic) private(i,j,k)
for(i=1;i<=nx;i++)
		{
		for(j=1;j<=ny;j++)
			{
            for(k=1;k<=nz;k++)
                {
                    shad2[i][j][k]=Res[i][j][k]-Relaxation*Ad[i][j][k];	// Updated pressure correction value
                }
            }
		}

// Matrix multiplication [A]*[s] for alp_denominator
om_num=0.0;
om_den=0.0;
#pragma omp parallel for schedule(dynamic) private(i,j,k,rhoe,rhow,rhos,rhon,rhot,rhob,apg,aeg,asg,awg,ang,atg,abg) reduction(+:om_num,om_den)
	for(i=1;i<=nx;i++)
		{
		for(j=1;j<=ny;j++)
			{
            for(k=1;k<=nz;k++)
                {
                    frho();
                     apg=(dt*dy*dz/(rhoe*dx))+(dt*dy*dz/(rhow*dx))+(dx*dz*dt/(rhon*dy))+(dt*dx*dz/(rhos*dy))+(dx*dy*dt/(rhot*dz))+(dt*dx*dy/(rhob*dz));
                    aeg=-(dt*dy*dz/(rhoe*dx));
                    awg=-(dt*dy*dz/(rhow*dx));
                    ang=-(dx*dz*dt/(rhon*dy));
                    asg=-(dt*dx*dz/(rhos*dy));
                    atg=-(dx*dy*dt/(rhot*dz));
                    abg=-(dt*dx*dy/(rhob*dz));

                    //Vertices
                    if(i==1 && j==1 && k==1)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+ang*shad2[i][j+1][k]+atg*shad2[i][j][k+1];
                    if(i==nx && j==1 && k==1)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+awg*shad2[i-1][j][k]+ang*shad2[i][j+1][k]+atg*shad2[i][j][k+1];
                    if(i==nx && j==ny && k==nz)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+awg*shad2[i-1][j][k]+asg*shad2[i][j-1][k]+abg*shad2[i][j][k-1];
                    if(i==1 && j==ny && k==1)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+asg*shad2[i][j-1][k]+atg*shad2[i][j][k+1];
                    if(i==1 && j==1 && k==nz)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+ang*shad2[i][j+1][k]+abg*shad2[i][j][k-1];
                    if(i==nx && j==ny && k==1)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+awg*shad2[i-1][j][k]+asg*shad2[i][j-1][k]+atg*shad2[i][j][k+1];
                    if(i==nx && j==1 && k==nz)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+awg*shad2[i-1][j][k]+ang*shad2[i][j+1][k]+abg*shad2[i][j][k-1];
                    if(i==1 && j==ny && k==nz)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+asg*shad2[i][j-1][k]+abg*shad2[i][j][k-1];
                    //Edges
                    if(i==1 && (j>1 && j<ny) && k==1)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+ang*shad2[i][j+1][k]+asg*shad2[i][j-1][k]+atg*shad2[i][j][k+1];
                    if(i==nx && (j>1 && j<ny) && k==1)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+awg*shad2[i-1][j][k]+ang*shad2[i][j+1][k]+asg*shad2[i][j-1][k]+atg*shad2[i][j][k+1];
                    if(i==nx && (j>1 && j<ny) && k==nz)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+awg*shad2[i-1][j][k]+ang*shad2[i][j+1][k]+asg*shad2[i][j-1][k]+abg*shad2[i][j][k-1];
                    if(i==1 && (j>1 && j<ny) && k==nz)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+ang*shad2[i][j+1][k]+asg*shad2[i][j-1][k]+abg*shad2[i][j][k-1];
                    if((i>1 && i<nx) && j==ny && k==1)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+awg*shad2[i-1][j][k]+asg*shad2[i][j-1][k]+atg*shad2[i][j][k+1];
                    if((i>1 && i<nx) && j==ny && k==nz)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+awg*shad2[i-1][j][k]+asg*shad2[i][j-1][k]+abg*shad2[i][j][k-1];
                    if((i>1 && i<nx) && j==1 && k==1)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+awg*shad2[i-1][j][k]+ang*shad2[i][j+1][k]+atg*shad2[i][j][k+1];
                    if((i>1 && i<nx) && j==1 && k==nz)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+awg*shad2[i-1][j][k]+ang*shad2[i][j+1][k]+abg*shad2[i][j][k-1];
                    if(i==nx && j==ny && (k>1 && k<nz))
                        Ashad[i][j][k]=apg*shad2[i][j][k]+awg*shad2[i-1][j][k]+asg*shad2[i][j-1][k]+atg*shad2[i][j][k+1]+abg*shad2[i][j][k-1];
                    if(i==1 && j==ny && (k>1 && k<nz))
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+asg*shad2[i][j-1][k]+atg*shad2[i][j][k+1]+abg*shad2[i][j][k-1];
                    if(i==1 && j==1 && (k>1 && k<nz))
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+ang*shad2[i][j+1][k]+atg*shad2[i][j][k+1]+abg*shad2[i][j][k-1];
                    if(i==nx && j==1 && (k>1 && k<nz))
                        Ashad[i][j][k]=apg*shad2[i][j][k]+awg*shad2[i-1][j][k]+ang*shad2[i][j+1][k]+atg*shad2[i][j][k+1]+abg*shad2[i][j][k-1];
                    // Faces
                    if((i>1 && i<nx) && j==1 && (k>1 && k<nz))
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+awg*shad2[i-1][j][k]+ang*shad2[i][j+1][k]+atg*shad2[i][j][k+1]+abg*shad2[i][j][k-1];
                    if((i>1 && i<nx) && j==ny && (k>1 && k<nz))
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+awg*shad2[i-1][j][k]+asg*shad2[i][j-1][k]+atg*shad2[i][j][k+1]+abg*shad2[i][j][k-1];
                    if(i==nx && (j>1 && j<ny) && (k>1 && k<nz))
                        Ashad[i][j][k]=apg*shad2[i][j][k]+awg*shad2[i-1][j][k]+asg*shad2[i][j-1][k]+ang*shad2[i][j+1][k]+atg*shad2[i][j][k+1]+abg*shad2[i][j][k-1];
                    if(i==1 && (j>1 && j<ny) && (k>1 && k<nz))
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+asg*shad2[i][j-1][k]+ang*shad2[i][j+1][k]+atg*shad2[i][j][k+1]+abg*shad2[i][j][k-1];
                    if((i>1 && i<nx) && (j>1 && j<ny) && k==nz)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+awg*shad2[i-1][j][k]+asg*shad2[i][j-1][k]+ang*shad2[i][j+1][k]+abg*shad2[i][j][k-1];
                    if((i>1 && i<nx) && (j>1 && j<ny) && k==1)
                        Ashad[i][j][k]=apg*shad2[i][j][k]+aeg*shad2[i+1][j][k]+awg*shad2[i-1][j][k]+asg*shad2[i][j-1][k]+ang*shad2[i][j+1][k]+atg*shad2[i][j][k+1];
                    //Interior
                    if((1<i<nx) && (1<j<ny) && (k>1 && k<nz))
                        Ashad[i][j][k]=apg*shad2[i][j][k]+awg*shad2[i-1][j][k]+aeg*shad2[i+1][j][k]+asg*shad2[i][j-1][k]+ang*shad2[i][j+1][k]+atg*shad2[i][j][k+1]+abg*shad2[i][j][k-1];
                    om_num=om_num+Ashad[i][j][k]*shad2[i][j][k];
                    om_den=om_den+Ashad[i][j][k]*Ashad[i][j][k];
                }
            }
		}
//Calculating omega
omega=om_num/om_den;
#pragma omp parallel for schedule(dynamic) private(i,j,k)
for(i=1;i<=nx;i++)
		{
		for(j=1;j<=ny;j++)
			{
            for(k=1;k<=nz;k++)
                {
                    pc[i][j][k]=pc[i][j][k]+Relaxation*direc[i][j][k]+omega*(shad2[i][j][k]);	// Updated pressure correction value
                    Res[i][j][k]=shad2[i][j][k]-omega*Ashad[i][j][k];	// New Residue
                }
			}
		}

b_den=alp_num;
alp_num=0.0;
#pragma omp parallel for schedule(dynamic) private(i,j,k) reduction(+:alp_num)
for(i=1;i<=nx;i++)
		{
		for(j=1;j<=ny;j++)
			{
            for(k=1;k<=nz;k++)
                {
                    alp_num=alp_num+Res[i][j][k]*Res_star[i][j][k];		// updated value for next iteration
                }
			}
		}

direc_fac=(alp_num/b_den)*(Relaxation/omega);
#pragma omp parallel for schedule(dynamic) private(i,j,k)
for(i=1;i<=nx;i++)
		{
		for(j=1;j<=ny;j++)
			{
            for(k=1;k<=nz;k++)
                {
                    direc[i][j][k]=Res[i][j][k]+direc_fac*(direc[i][j][k]-omega*Ad[i][j][k]);	// updated value for the next iteration
                }
			}
		}
gRMS=sqrt((alp_num/(nx*ny*nz)));
gcount++;

}while(gRMS>=0.000001);

#pragma omp parallel for schedule(dynamic) private(i,j,k)
for(i=img;i<img+nx;i++)
		{
		for(j=img;j<img+ny;j++)
			{
			for(k=img;k<img+nz;k++)
				{
					pnew[i][j][k]=0.5*(pr[i][j][k]+pr1[i][j][k])+pc[i][j][k];
					pc[i][j][k]=0.0;
				}
			}
		}
BC_cor();
#pragma omp parallel for schedule(dynamic) private(i,j,k,diff,sqr,rhot,rhob,rhoe,rhow,rhos,rhon) reduction(+:usum)
	for(i=img;i<img+nx;i++)
		{
		for(j=img;j<img+ny;j++)
			{
			for(k=img;k<img+nz;k++)
				{

					//Face density
					rhoe=(rho[i+1][j][k]+rho[i][j][k])/2.0;
					rhow=(rho[i-1][j][k]+rho[i][j][k])/2.0;
					rhon=(rho[i][j+1][k]+rho[i][j][k])/2.0;
					rhos=(rho[i][j-1][k]+rho[i][j][k])/2.0;
					rhot=(rho[i][j][k]+rho[i][j][k+1])/2.0;
					rhob=(rho[i][j][k]+rho[i][j][k-1])/2.0;

					uc[i][j][k]=0.5*(dt/dx)*((pnew[i+1][j][k]/rhoe)-(pnew[i-1][j][k]/rhow)+pnew[i][j][k]*(1.0/rhow-1.0/rhoe));
					vc[i][j][k]=0.5*(dt/dy)*((pnew[i][j+1][k]/rhon)-(pnew[i][j-1][k]/rhos)+pnew[i][j][k]*(1.0/rhos-1.0/rhon));
					wc[i][j][k]=0.5*(dt/dz)*((pnew[i][j][k+1]/rhot)-(pnew[i][j][k-1]/rhob)+pnew[i][j][k]*(1.0/rhob-1.0/rhot));
					unew[i][j][k]=ut[i][j][k]-uc[i][j][k];
					vnew[i][j][k]=vt[i][j][k]-vc[i][j][k];
					wnew[i][j][k]=wwt[i][j][k]-wc[i][j][k];
					//RMS for steady state
					diff=unew[i][j][k]-u[i][j][k];
					sqr=diff*diff;
					usum=usum+sqr;

					//Updating the old values
					u[i][j][k]=unew[i][j][k];
					v[i][j][k]=vnew[i][j][k];
					w[i][j][k]=wnew[i][j][k];
					pr[i][j][k]=pnew[i][j][k];
				//	printf("i=%d  j=%d  uc=%f\tvc=%f\tp=%f\n",i,j,uc[i][j],vc[i][j],pnew[i][j]);
				}
			}
		}

		vmax=0.0;
		umax=0.0;
		wmax=0.0;
		pmax=0.0;
		for(i=img;i<nx+img;i++)
				{
				for(j=img;j<ny+img;j++)
					{
					for(k=img;k<nz+img;k++)
						{
							if(fabs(v[i][j][k])>vmax)
								{
									vmax=v[i][j][k];

								}
							if(fabs(u[i][j][k])>umax)
								{
									umax=u[i][j][k];
								}
							if(fabs(w[i][j][k])>wmax)
								{
									wmax=w[i][j][k];
								}
							if(fabs(pr[i][j][k])>pmax)
							{
								pmax=pr[i][j][k];
							}
						}
					}
				}
double co;
		//		------------------------------Courant number calculations-----------------------------------------
		if(vmax>umax && vmax>wmax)
		co=vmax*dt/dy;
		else if (umax>vmax && umax>wmax)
		co=umax*dt/dx;
		else
		co=wmax*dt/dz;

		if(co>0.010)
			dt=0.05*dt;
		else if(co<0.01)
			dt=1.02*dt;


    	uRMS=sqrt(usum/(nx*ny*nz)); // u-vel RMS
   	printf("\tCourant=%f\tdt (iter= %d)=%f\n",co,ct,dt);
	printf("\tu-max=%f\t\tv-max=%f\tw-max=%f\t p-max=%f\t uRMS=%f\n\n",umax,vmax,wmax,pmax,uRMS);
	printf("Time= %f \n",t);

