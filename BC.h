void BC1()
{
// East and West Face
for(j=img;j<ny+img;j++)
	{
		for(k=img;k<nz+img;k++)
		{

			Kurv[img-1][j][k]=Kurv[img][j][k];
			Kurv[img+nx][j][k]=Kurv[img+nx-1][j][k];
			phi[img-1][j][k]=phi[img][j][k];
			phi[img+nx][j][k]=phi[img+nx-1][j][k];
			phi_t[img-1][j][k]=phi_t[img][j][k];
			phi_t[img+nx][j][k]=phi_t[img+nx-1][j][k];
		}
	}
// Top and Bottom Face
for(i=img;i<nx+img;i++)
{
	for(k=img;k<nz+img;k++)
	{
		Kurv[i][img-1][k]=Kurv[i][img][k];
		Kurv[i][ny+img][k]=Kurv[i][ny+img-1][k];

	if((i>=img+n_L1+1 && i<=(img+n_L1+n_nozx)) &&(k>=img+n_L2+1 && k<=(img+n_L2+n_nozz))) //Inlet
		{
			phi[i][img+ny][k]=2.0-phi[i][img+ny-1][k];
			phi_t[i][img+ny][k]=2.0-phi_t[i][img+ny-1][k];
		}
		else
		{
			phi[i][img+ny][k]=phi[i][ny+img-1][k];
			phi_t[i][img+ny][k]=phi_t[i][ny+img-1][k];
		}
		phi[i][img-1][k]=phi[i][img][k];
		phi_t[i][img-1][k]=phi_t[i][img][k];

	}
}

// Front and Back Face
for(i=img;i<nx+img;i++)
{
	for(j=img;j<ny+img;j++)
	{
		Kurv[i][j][img-1]=Kurv[i][j][img];
		Kurv[i][j][nz+img]=Kurv[i][j][nz+img-1];

		phi[i][j][img-1]=phi[i][j][img];
		phi[i][j][nz+img]=phi[i][j][nz+img-1];
		phi_t[i][j][img-1]=phi_t[i][j][img];
		phi_t[i][j][nz+img]=phi_t[i][j][nz+img-1];
	}
}
}


void BC2()
{
	for(j=img;j<ny+img;j++)
			{
				for(k=img;k<nz+img;k++)
				{
					u[img-1][j][k]=-u[img][j][k];
					u[img-2][j][k]=2.0*u[img-1][j][k]-u[img][j][k];
					v[img-1][j][k]=v[img][j][k];
					v[img-2][j][k]=2.0*v[img-1][j][k]-v[img][j][k];
					w[img-1][j][k]=w[img][j][k];
					w[img-2][j][k]=2.0*w[img-1][j][k]-w[img][j][k];

					u[nx+img][j][k]=-u[nx+img-1][j][k];
					u[nx+img+1][j][k]=2.0*u[nx+img][j][k]-u[nx+img-1][j][k];
					v[nx+img][j][k]=v[nx+img-1][j][k];
					v[nx+img+1][j][k]=2.0*v[nx+img][j][k]-v[nx+img-1][j][k];
					w[nx+img][j][k]=w[nx+img-1][j][k];
					w[nx+img+1][j][k]=2.0*w[nx+img][j][k]-w[nx+img-1][j][k];
					rho[img-1][j][k]=rho[img][j][k];
					visc[img-1][j][k]=visc[img][j][k];

					rho[nx+img][j][k]=rho[nx+img-1][j][k];
					visc[nx+img][j][k]=visc[nx+img-1][j][k];
				}
			}

			for(i=img;i<nx+img;i++)
			{
				for(k=img;k<nz+img;k++)
				{

					if((i>=img+n_L1+1 && i<=(img+n_L1+n_nozx)) &&(k>=img+n_L2+1 && k<=(img+n_L2+n_nozz)))//Inlet
					{
						v0=-1.5*v_m*(1.0-pow((x[i]-L/2.0),2.0)/(noz*noz));
						v[i][img+ny][k]=2.0*(v0)-v[i][ny+img-1][k];
						u[i][img+ny][k]=2.0*Amp*cos(freq*t)-u[i][img+ny-1][k];
					}
					else
					{
						v[i][img+ny][k]=-v[i][ny+img-1][k];
						u[i][img+ny][k]=-u[i][ny+img-1][k]; //Top Wall
					}

					u[i][img-1][k]=u[i][img][k];
					u[i][img-2][k]=2.0*u[i][img-1][k]-u[i][img][k];
					v[i][img-1][k]=v[i][img][k];
					v[i][img-2][k]=2.0*v[i][img-1][k]-v[i][img][k];
					w[i][img-1][k]=w[i][img][k];
					w[i][img-2][k]=2.0*w[i][img-1][k]-w[i][img][k];

					u[i][img+ny+1][k]=2.0*u[i][img+ny][k]-u[i][img+ny-1][k];
					v[i][img+ny+1][k]=2.0*v[i][img+ny][k]-v[i][img+ny-1][k];
					w[i][img+ny][k]=-w[i][img+ny-1][k];
					w[i][img+ny+1][k]=2.0*w[i][img+ny][k]-w[i][img+ny-1][k];

					rho[i][img-1][k]=rho[i][img][k];
					visc[i][img-1][k]=visc[i][img][k];

					rho[i][img+ny][k]=rho[i][img+ny-1][k];
					visc[i][img+ny][k]=visc[i][img+ny-1][k];
				}
			}

			for(i=img;i<img+nx;i++)
			{
				for(j=img;j<img+ny;j++)
				{
					u[i][j][img-1]=u[i][j][img];
					u[i][j][img-2]=2.0*u[i][j][img-1]-u[i][j][img];
					v[i][j][img-1]=v[i][j][img];
					v[i][j][img-2]=2.0*v[i][j][img-1]-v[i][j][img];
					w[i][j][img-1]=-w[i][j][img];
					w[i][j][img-2]=2.0*w[i][j][img-1]-w[i][j][img];

					u[i][j][img+nz]=u[i][j][img+nz-1];
					u[i][j][img+nz+1]=2.0*u[i][j][nz+img]-u[i][j][nz+img-1];
					v[i][j][img+nz]=v[i][j][img+nz-1];
					v[i][j][img+nz+1]=2.0*v[i][j][nz+img]-v[i][j][nz+img-1];
					w[i][j][img+nz]=-w[i][j][img+nz-1];
					w[i][j][img+nz+1]=2.0*w[i][j][nz+img]-w[i][j][nz+img-1];

					rho[i][j][img-1]=rho[i][j][img];
					visc[i][j][img-1]=visc[i][j][img];

					rho[i][j][nz+img]=rho[i][j][img+nz-1];
					visc[i][j][nz+img]=visc[i][j][img+nz-1];
				}
			}


}
void BC_in()
{

for(j=img;j<ny+img;j++)
			{
				for(k=img;k<nz+img;k++)
				{
					u1[img-1][j][k]=-u1[img][j][k];
					u1[img-2][j][k]=2.0*u1[img-1][j][k]-u1[img][j][k];
					v1[img-1][j][k]=v1[img][j][k];
					v1[img-2][j][k]=2.0*v1[img-1][j][k]-v1[img][j][k];
					w1[img-1][j][k]=w1[img][j][k];
					w1[img-2][j][k]=2.0*w1[img-1][j][k]-w1[img][j][k];

					u1[nx+img][j][k]=-u1[nx+img-1][j][k];
					u1[nx+img+1][j][k]=2.0*u1[nx+img][j][k]-u1[nx+img-1][j][k];
					v1[nx+img][j][k]=v1[nx+img-1][j][k];
					v1[nx+img+1][j][k]=2.0*v1[nx+img][j][k]-v1[nx+img-1][j][k];
					w1[nx+img][j][k]=w1[nx+img-1][j][k];
					w1[nx+img+1][j][k]=2.0*w1[nx+img][j][k]-w1[nx+img-1][j][k];
					rho[img-1][j][k]=rho[img][j][k];
					visc[img-1][j][k]=visc[img][j][k];

					rho[nx+img][j][k]=rho[nx+img-1][j][k];
					visc[nx+img][j][k]=visc[nx+img-1][j][k];
				}
			}

			for(i=img;i<nx+img;i++)
			{
				for(k=img;k<nz+img;k++)
				{

					if((i>=img+n_L1+1 && i<=(img+n_L1+n_nozx)) &&(k>=img+n_L2+1 && k<=(img+n_L2+n_nozz)))
					{
						v0=-1.5*v_m*(1.0-pow((x[i]-L/2.0),2.0)/(noz*noz));
						v1[i][img+ny][k]=2.0*(v0)-v1[i][ny+img-1][k];
						u1[i][img+ny][k]=2.0*Amp*cos(freq*t)-u1[i][img+ny-1][k];
					}
					else
					{
						v1[i][img+ny][k]=-v1[i][ny+img-1][k];
						u1[i][img+ny][k]=-u1[i][ny+img-1][k];
					}

					u1[i][img-1][k]=u1[i][img][k];
					u1[i][img-2][k]=2.0*u1[i][img-1][k]-u1[i][img][k];
					v1[i][img-1][k]=v1[i][img][k];
					v1[1i][img-2][k]=2.0*v1[i][img-1][k]-v1[i][img][k];
					w1[i][img-1][k]=w1[i][img][k];
					w1[i][img-2][k]=2.0*w1[i][img-1][k]-w1[i][img][k];

					u1[i][img+ny+1][k]=2.0*u1[i][img+ny][k]-u1[i][img+ny-1][k];
					v1[i][img+ny+1][k]=2.0*v1[i][img+ny][k]-v1[i][img+ny-1][k];
					w1[i][img+ny][k]=-w1[i][img+ny-1][k];
					w1[i][img+ny+1][k]=2.0*w1[i][img+ny][k]-w1[i][img+ny-1][k];

					rho[i][img-1][k]=rho[i][img][k];
					visc[i][img-1][k]=visc[i][img][k];

					rho[i][img+ny][k]=rho[i][img+ny-1][k];
					visc[i][img+ny][k]=visc[i][img+ny-1][k];
				}
			}

			for(i=img;i<img+nx;i++)
			{
				for(j=img;j<img+ny;j++)
				{

					u1[i][j][img-1]=u1[i][j][img];
					u1[i][j][img-2]=2.0*u1[i][j][img-1]-u1[i][j][img];
					v1[i][j][img-1]=v1[i][j][img];
					v1[i][j][img-2]=2.0*v1[i][j][img-1]-v1[i][j][img];
					w1[i][j][img-1]=-w1[i][j][img];
					w1[i][j][img-2]=2.0*w1[i][j][img-1]-w1[i][j][img];

					u1[i][j][img+nz]=u1[i][j][img+nz-1];
					u1[i][j][img+nz+1]=2.0*u1[i][j][nz+img]-u1[i][j][nz+img-1];
					v1[i][j][img+nz]=v1[i][j][img+nz-1];
					v1[i][j][img+nz+1]=2.0*v1[i][j][nz+img]-v1[i][j][nz+img-1];
					w1[i][j][img+nz]=-w1[i][j][img+nz-1];
					w1[i][j][img+nz+1]=2.0*w1[i][j][nz+img]-w1[i][j][nz+img-1];

					rho[i][j][img-1]=rho[i][j][img];
					visc[i][j][img-1]=visc[i][j][img];

					rho[i][j][nz+img]=rho[i][j][img+nz-1];
					visc[i][j][nz+img]=visc[i][j][img+nz-1];
				}
			}

}

void BC2t()
{



for(j=img;j<ny+img;j++)
{
	for(k=img;k<nz+img;k++)
	{
		ut[img-1][j][k]=-ut[img][j][k];
		vt[img-1][j][k]=vt[img][j][k];
		wwt[img-1][j][k]=wwt[img][j][k];

		ut[nx+img][j][k]=-ut[nx+img-1][j][k];
		vt[nx+img][j][k]=vt[nx+img-1][j][k];
		wwt[nx+img][j][k]=wwt[nx+img-1][j][k];

		rho[img-1][j][k]=rho[img][j][k];
		rho[nx+img][j][k]=rho[nx+img-1][j][k];
		pr[img-1][j][k]=pr[img][j][k];
		pr[nx+img][j][k]=pr[nx+img-1][j][k];
	}
}

for(i=img;i<nx+img;i++)
{
	for(k=img;k<nz+img;k++)
	{
		if((i>=img+n_L1+1 && i<=(img+n_L1+n_nozx)) &&(k>=img+n_L2+1 && k<=(img+n_L2+n_nozz)))
		{
			v0=-1.5*v_m*(1.0-pow((x[i]-L/2.0),2.0)/(noz*noz));
			vt[i][img+ny][k]=2.0*(v0)-vt[i][ny+img-1][k];
			ut[i][img+ny][k]=2.0*Amp*cos(freq*t)-ut[i][img+ny-1][k];

		}
		else
		{
			vt[i][img+ny][k]=-vt[i][ny+img-1][k];
			ut[i][img+ny][k]=-ut[i][ny+img-1][k]; //Top Wall
		}

		ut[i][img-1][k]=ut[i][img][k];
		vt[i][img-1][k]=vt[i][img][k];
		wwt[i][img-1][k]=wwt[i][img][k];

		wwt[i][img+ny][k]=-wwt[i][img+ny-1][k];
		pr[i][img+ny][k]=pr[i][img+ny-1][k];
		rho[i][img-1][k]=rho[i][img][k];
		rho[i][img+ny][k]=rho[i][img+ny-1][k];
		pr[i][img-1][k]=-pr[i][img][k];
	}
}

for(i=img;i<img+nx;i++)
{
	for(j=img;j<img+ny;j++)
	{
		ut[i][j][img-1]=ut[i][j][img];
		vt[i][j][img-1]=vt[i][j][img];
		wwt[i][j][img-1]=-wwt[i][j][img];

		ut[i][j][img+nz]=ut[i][j][img+nz-1];
		vt[i][j][img+nz]=vt[i][j][img+nz-1];
		wwt[i][j][img+nz]=-wwt[i][j][img+nz-1];

		rho[i][j][img-1]=rho[i][j][img];
		rho[i][j][nz+img]=rho[i][j][img+nz-1];
		pr[i][j][img-1]=pr[i][j][img];
		pr[i][j][nz+img]=pr[i][j][nz+img-1];
	}
}

}
void BC_t_in()
{


for(j=img;j<ny+img;j++)
{
	for(k=img;k<nz+img;k++)
	{
		ut[img-1][j][k]=-ut[img][j][k];
		vt[img-1][j][k]=vt[img][j][k];
		wwt[img-1][j][k]=wwt[img][j][k];

		ut[nx+img][j][k]=-ut[nx+img-1][j][k];
		vt[nx+img][j][k]=vt[nx+img-1][j][k];
		wwt[nx+img][j][k]=wwt[nx+img-1][j][k];

		rho[img-1][j][k]=rho[img][j][k];
		rho[nx+img][j][k]=rho[nx+img-1][j][k];
		pr[img-1][j][k]=pr[img][j][k];
		pr[nx+img][j][k]=pr[nx+img-1][j][k];
        pr1[img-1][j][k]=pr1[img][j][k];
		pr1[nx+img][j][k]=pr1[nx+img-1][j][k];
	}
}

for(i=img;i<nx+img;i++)
{
	for(k=img;k<nz+img;k++)
	{
		if((i>=img+n_L1+1 && i<=(img+n_L1+n_nozx)) &&(k>=img+n_L2+1 && k<=(img+n_L2+n_nozz)))
		{
			v0=-1.5*v_m*(1.0-pow((x[i]-L/2.0),2.0)/(noz*noz));
			vt[i][img+ny][k]=2.0*(v0)-vt[i][ny+img-1][k];
			ut[i][img+ny][k]=2.0*Amp*cos(freq*t)-ut[i][img+ny-1][k];

		}
		else
		{
			vt[i][img+ny][k]=-vt[i][ny+img-1][k];
			ut[i][img+ny][k]=-ut[i][ny+img-1][k];
		}

		ut[i][img-1][k]=ut[i][img][k];
		vt[i][img-1][k]=vt[i][img][k];
		wwt[i][img-1][k]=wwt[i][img][k];

		wwt[i][img+ny][k]=-wwt[i][img+ny-1][k];
		pr[i][img+ny][k]=pr[i][img+ny-1][k];
        pr1[i][img+ny][k]=pr1[i][img+ny-1][k];
		rho[i][img-1][k]=rho[i][img][k];
		rho[i][img+ny][k]=rho[i][img+ny-1][k];
		pr[i][img-1][k]=-pr[i][img][k];
        pr1[i][img-1][k]=-pr1[i][img][k];
	}
}

for(i=img;i<img+nx;i++)
{
	for(j=img;j<img+ny;j++)
	{

		ut[i][j][img-1]=ut[i][j][img];
		vt[i][j][img-1]=vt[i][j][img];
		wwt[i][j][img-1]=-wwt[i][j][img];

		ut[i][j][img+nz]=ut[i][j][img+nz-1];
		vt[i][j][img+nz]=vt[i][j][img+nz-1];
		wwt[i][j][img+nz]=-wwt[i][j][img+nz-1];

		rho[i][j][img-1]=rho[i][j][img];
		rho[i][j][nz+img]=rho[i][j][img+nz-1];
		pr[i][j][img-1]=pr[i][j][img];
		pr[i][j][nz+img]=pr[i][j][nz+img-1];
        pr1[i][j][img-1]=pr1[i][j][img];
		pr1[i][j][nz+img]=pr1[i][j][nz+img-1];
	}
}

}
void BC_pc()
{

for(j=img;j<ny+img;j++)
{
	for(k=img;k<nz+img;k++)
	{
		rho[img-1][j][k]=rho[img][j][k];
		rho[nx+img][j][k]=rho[nx+img-1][j][k];
		pc[img-1][j][k]=pc[img][j][k];
		pc[nx+img][j][k]=pc[nx+img-1][j][k];
	}
}
for(i=img;i<nx+img;i++)
{
	for(k=img;k<nz+img;k++)
	{
		rho[i][img-1][k]=rho[i][img][k];
		rho[i][img+ny][k]=rho[i][img+ny-1][k];
		pc[i][img-1][k]=-pc[i][img][k];
		pc[i][img+ny][k]=pc[i][img+ny-1][k];
	}
}
for(i=img;i<img+nx;i++)
{
	for(j=img;j<img+ny;j++)
	{
		rho[i][j][img-1]=rho[i][j][img];
		rho[i][j][nz+img]=rho[i][j][img+nz-1];
		pc[i][j][img-1]=pc[i][j][img];
		pc[i][j][nz+img]=pc[i][j][nz+img-1];
	}
}
}

void BC_cor()
{

for(j=img;j<ny+img;j++)
{
	for(k=img;k<nz+img;k++)
	{
		uc[img-1][j][k]=-uc[img][j][k];
		vc[img-1][j][k]=vc[img][j][k];
		wc[img-1][j][k]=wc[img][j][k];
		uc[nx+img][j][k]=-uc[nx+img-1][j][k];
		vc[nx+img][j][k]=vc[nx+img-1][j][k];
		wc[nx+img][j][k]=wc[nx+img-1][j][k];

		rho[img-1][j][k]=rho[img][j][k];
		rho[nx+img][j][k]=rho[nx+img-1][j][k];
		pnew[img-1][j][k]=pnew[img][j][k];
		pnew[nx+img][j][k]=pnew[nx+img-1][j][k];
	}
}

for(i=img;i<nx+img;i++)
{
	for(k=img;k<nz+img;k++)
	{
		uc[i][img-1][k]=uc[i][img][k];
		vc[i][img-1][k]=vc[i][img][k];
		wc[i][img-1][k]=wc[i][img][k];

		uc[i][img+ny][k]=-uc[i][img+ny-1][k];
		vc[i][img+ny][k]=-vc[i][img+ny-1][k];
		wc[i][img+ny][k]=-wc[i][img+ny-1][k];

		rho[i][img-1][k]=rho[i][img][k];
		rho[i][img+ny][k]=rho[i][img+ny-1][k];
		pnew[i][img-1][k]=-pnew[i][img][k];
		pnew[i][img+ny][k]=pnew[i][img+ny-1][k];
	}
}

for(i=img;i<img+nx;i++)
{
	for(j=img;j<img+ny;j++)
	{
		uc[i][j][img-1]=uc[i][j][img];
		vc[i][j][img-1]=vc[i][j][img];
		wc[i][j][img-1]=-wc[i][j][img];

		uc[i][j][img+nz]=uc[i][j][img+nz-1];
		vc[i][j][img+nz]=vc[i][j][img+nz-1];
		wc[i][j][img+nz]=-wc[i][j][img+nz-1];

		rho[i][j][img-1]=rho[i][j][img];
		rho[i][j][nz+img]=rho[i][j][img+nz-1];
		pnew[i][j][img-1]=pnew[i][j][img];
		pnew[i][j][nz+img]=pnew[i][j][nz+img-1];
	}
}
}
