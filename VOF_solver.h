void VOF()
{
//_______________________________VOF_______________________________________________//
		do
		{
			sum=0.0;
			//.........................boundary condition for velocity and vof (old and new).................................//
			//--------------------------- in x direction --------------------------------------------//
			for(j=img;j<img+ny;j++)
			{
				for(k=img;k<img+nz;k++)
				{
					//----------------------- west face ----------------------------------//
					u[img-1][j][k]=-u[img][j][k];
					v[img-1][j][k]=v[img][j][k];
					w[img-1][j][k]=w[img][j][k];
					phi[img-1][j][k]=phi[img][j][k];
					phi_t[img-1][j][k]=phi_t[img][j][k];
					//----------------------- east face -----------------------------------//
					u[nx+img][j][k]=-u[nx+img-1][j][k];
					v[nx+img][j][k]=v[nx+img-1][j][k];
					w[nx+img][j][k]=w[nx+img-1][j][k];
					phi[nx+img][j][k]=phi[nx+img-1][j][k];
					phi_t[nx+img][j][k]=phi_t[nx+img-1][j][k];
				}
			}
			//-------------------------------- in y direction --------------------------------//
			for(i=img;i<img+nx;i++)
				{
				for(k=img;k<img+nz;k++)
				{
					//---------------------- north face -----------------------------//
					if((i>=img+n_L1+1 && i<=(img+n_L1+n_nozx)) &&(k>=img+n_L2+1 && k<=(img+n_L2+n_nozz)))
					{
						v0=-2.0*v_m*(1.0-pow((x[i]-L/2.0),2.0)/(noz*noz));
						v[i][img+ny][k]=2.0*(v0)-v[i][ny+img-1][k];   // With Sinusoidal perturbations
						phi[i][img+ny][k]=2.0-phi[i][img+ny-1][k];
						phi_t[i][img+ny][k]=2.0-phi_t[i][img+ny-1][k];
						u[i][img+ny][k]=2.0*Amp*cos(freq*t)-u[i][img+ny-1][k]; //Top Wall
					}
					else
					{
						phi[i][img+ny][k]=phi[i][ny+img-1][k];
						phi_t[i][img+ny][k]=phi_t[i][ny+img-1][k];
						v[i][img+ny][k]=-v[i][ny+img-1][k];
						u[i][img+ny][k]=-u[i][ny+img-1][k]; //Top Wall
					}
					w[i][img+ny][k]=-w[i][ny+img-1][k];

					//--------------------- south face ------------------------------//
					w[i][img-1][k]=w[i][img][k];
					v[i][img-1][k]=v[i][img][k];
					u[i][img-1][k]=u[i][img][k];
					phi[i][img-1][k]=phi[i][img][k];
					phi_t[i][img-1][k]=phi_t[i][img][k];
				}
			}

			//----------------------------- in z direction ----------------------------------------//
			for(i=img;i<img+nx;i++)
				{
				for(j=img;j<img+ny;j++)
				{
					//---------------------------- bottom face -----------------------------//
					u[i][j][img-1]=u[i][j][img];
					v[i][j][img-1]=v[i][j][img];
					w[i][j][img-1]=-w[i][j][img];
					phi[i][j][img-1]=phi[i][j][img];
					phi_t[i][j][img-1]=phi_t[i][j][img];
					//----------------------------- top face ------------------------------//
					u[i][j][img+nz]=u[i][j][img+nz-1];
					v[i][j][img+nz]=v[i][j][img+nz-1];
					w[i][j][img+nz]=-w[i][j][img+nz-1];
					phi[i][j][img+nz]=phi[i][j][img+nz-1];
					phi_t[i][j][img+nz]=phi_t[i][j][img+nz-1];
				}
			}
#pragma omp parallel for schedule(dynamic) private(i,j,k,atd,atfc,atfu,theta,yy,b,atil,ue,uw,vn,vs,wt,wb,ap,ad,aa,au,ae,aw,an,as,at,ab,diva,mdiva,be,bw,bs,bn,bt,bb,ase,asw,asn,ass,ast,asb,res) reduction(+:sum)
			for(i=img;i<nx+img;i++)
			{
				for(j=img;j<ny+img;j++)
				{
					for(k=img;k<nz+img;k++)
					{
                        // Advecting velocity at faces
						ue=0.5*(u[i][j][k]+u[i+1][j][k]);
						uw=0.5*(u[i][j][k]+u[i-1][j][k]);
						vn=0.5*(v[i][j][k]+v[i][j+1][k]);
						vs=0.5*(v[i][j][k]+v[i][j-1][k]);
						wt=0.5*(w[i][j][k]+w[i][j][k+1]);
						wb=0.5*(w[i][j][k]+w[i][j][k-1]);
						ap=0.0;

						if(ue>0.0)
						{
                            //Coefficients Matrix
							ad=phi[i][j][k];
							aa=phi[i+1][j][k];
							au=phi[i-1][j][k];
							ae=0.5*(phi[i][j][k]+phi[i+1][j][k]);
							aw=0.5*(phi[i][j][k]+phi[i-1][j][k]);
							an=0.5*(phi[i][j][k]+phi[i][j+1][k]);
							as=0.5*(phi[i][j][k]+phi[i][j-1][k]);
							at=0.5*(phi[i][j][k]+phi[i][j][k+1]);
							ab=0.5*(phi[i][j][k]+phi[i][j][k-1]);
							diva=ae-aw;
							mdiva=sqrt(((ae-aw)/dx)*((ae-aw)/dx)+((an-as)/dy)*((an-as)/dy)+((at-ab)/dz)*((at-ab)/dz))*sqrt(dx*dx);

						//	courant(i,j,k);
                           if(ad!=1.0 || fabs(ad)>eps)
                            {
                                //****************tilda donor calculation***********************//
                                if(fabs(aa-au)<eps)
                                {
                                    atd=ad;
                                }
                                else
                                {
                                    atd=(ad-au)/(aa-au);
                                }
                                //****************SAISH******************************//
                                //****************convective boundedness scheme*****************//
                                if(atd>=0.0 && atd<=1.0/4.0)
                                    atfc=4.0*atd;
                                else if(atd>1.0/4.0  &&  atd<=1.0)
                                    atfc=1.0;
                                else
                                    atfc=atd;

                                //****************Hybrid HR******************************//
                                if(atd>=0.0 && atd<=1.0/6.0)
                                    atfu=11.0/6.0*atd;
                                else if(atd>1.0/6.0 && atd<=1.0/2.0)
                                     atfu=atd*(2.0-atd);			//CLAM
                                else if(atd>1.0/2.0 && atd<=3.0/4.0)
                                    atfu=atd+0.25;					// Fromm's
                                else if(atd>3.0/4.0 && atd<=1.0)
                                    atfu=1.0;
                                else
                                    atfu=atd;

                                theta=(fabs(diva/mdiva));
                                yy=fmin(pow(theta,2.0),1.0);
                                //****************tilda face calculation************************//
                                atil=yy*atfc+(1-yy)*atfu;
                                //******************************weighting factor calculation****************//
                                if(fabs(1.0-atd)<eps)
                                {
                                    b=0.5;
                                }
                                else
                                {
                                    b=(atil-atd)/(1.0-atd);
                                }
                            }
                            else
                            {
                                b=0.0;
                            }
							be=b;
							ase=((1.0-be)*(phi[i][j][k]+phi_t[i][j][k])*0.5)+(0.5*be*(phi[i+1][j][k]+phi_t[i+1][j][k]));
							ap=ap-(0.5*(1.0-be)*ue*(dt/dx));
						}
						else
						{
							ad=phi[i+1][j][k];
							aa=phi[i][j][k];
							au=phi[i+2][j][k];
							ae=0.5*(phi[i+1][j][k]+phi[i+2][j][k]);
							aw=0.5*(phi[i+1][j][k]+phi[i][j][k]);
							an=0.5*(phi[i+1][j][k]+phi[i+1][j+1][k]);
							as=0.5*(phi[i+1][j][k]+phi[i+1][j-1][k]);
							at=0.5*(phi[i+1][j][k]+phi[i+1][j][k+1]);
							ab=0.5*(phi[i+1][j][k]+phi[i+1][j][k-1]);
							diva=ae-aw;
							mdiva=sqrt(((ae-aw)/dx)*((ae-aw)/dx)+((an-as)/dy)*((an-as)/dy)+((at-ab)/dz)*((at-ab)/dz))*sqrt(dx*dx);

						//	courant(i+1,j,k);
						  if(ad!=1.0 || fabs(ad)>eps)
                            {
                                //****************tilda donor calculation***********************//
                                if(fabs(aa-au)<eps)
                                {
                                    atd=ad;
                                }
                                else
                                {
                                    atd=(ad-au)/(aa-au);
                                }
                                //****************SAISH******************************//
                                //****************convective boundedness scheme*****************//
                                if(atd>=0.0 && atd<=1.0/4.0)
                                    atfc=4.0*atd;
                                else if(atd>1.0/4.0  &&  atd<=1.0)
                                    atfc=1.0;
                                else
                                    atfc=atd;

                                //****************Hybrid HR******************************//
                                if(atd>=0.0 && atd<=1.0/6.0)
                                    atfu=11.0/6.0*atd;
                                else if(atd>1.0/6.0 && atd<=1.0/2.0)
                                     atfu=atd*(2.0-atd);			//CLAM
                                else if(atd>1.0/2.0 && atd<=3.0/4.0)
                                    atfu=atd+0.25;					// Fromm's
                                else if(atd>3.0/4.0 && atd<=1.0)
                                    atfu=1.0;
                                else
                                    atfu=atd;

                                theta=(fabs(diva/mdiva));
                                yy=fmin(pow(theta,2.0),1.0);
                                //****************tilda face calculation************************//
                                atil=yy*atfc+(1-yy)*atfu;
                                //******************************weighting factor calculation****************//
                                if(fabs(1.0-atd)<eps)
                                {
                                    b=0.5;
                                }
                                else
                                {
                                    b=(atil-atd)/(1.0-atd);
                                }
                            }
                            else
                            {
                                b=0.0;
                            }
							be=b;
							ase=((1.0-be)*0.5*(phi[i+1][j][k]+phi_t[i+1][j][k]))+(0.5*be*(phi[i][j][k]+phi_t[i][j][k]));
							ap=ap-(0.5*be*ue*(dt/dx));
						}
						if(uw>0.0)
						{
							ad=phi[i-1][j][k];
							aa=phi[i][j][k];
							au=phi[i-2][j][k];
							ae=0.5*(phi[i-1][j][k]+phi[i][j][k]);
							aw=0.5*(phi[i-1][j][k]+phi[i-2][j][k]);
							an=0.5*(phi[i-1][j][k]+phi[i-1][j+1][k]);
							as=0.5*(phi[i-1][j][k]+phi[i-1][j-1][k]);
							at=0.5*(phi[i-1][j][k]+phi[i-1][j][k+1]);
							ab=0.5*(phi[i-1][j][k]+phi[i-1][j][k-1]);
							diva=ae-aw;
							mdiva=sqrt(((ae-aw)/dx)*((ae-aw)/dx)+((an-as)/dy)*((an-as)/dy)+((at-ab)/dz)*((at-ab)/dz))*sqrt(dx*dx);

						//	courant(i-1,j,k);
						  if(ad!=1.0 || fabs(ad)>eps)
                            {
                                //****************tilda donor calculation***********************//
                                if(fabs(aa-au)<eps)
                                {
                                    atd=ad;
                                }
                                else
                                {
                                    atd=(ad-au)/(aa-au);
                                }
                                //****************SAISH******************************//
                                //****************convective boundedness scheme*****************//
                                if(atd>=0.0 && atd<=1.0/4.0)
                                    atfc=4.0*atd;
                                else if(atd>1.0/4.0  &&  atd<=1.0)
                                    atfc=1.0;
                                else
                                    atfc=atd;

                                //****************Hybrid HR******************************//
                                if(atd>=0.0 && atd<=1.0/6.0)
                                    atfu=11.0/6.0*atd;
                                else if(atd>1.0/6.0 && atd<=1.0/2.0)
                                     atfu=atd*(2.0-atd);			//CLAM
                                else if(atd>1.0/2.0 && atd<=3.0/4.0)
                                    atfu=atd+0.25;					// Fromm's
                                else if(atd>3.0/4.0 && atd<=1.0)
                                    atfu=1.0;
                                else
                                    atfu=atd;

                                theta=(fabs(diva/mdiva));
                                yy=fmin(pow(theta,2.0),1.0);
                                //****************tilda face calculation************************//
                                atil=yy*atfc+(1-yy)*atfu;
                                //******************************weighting factor calculation****************//
                                if(fabs(1.0-atd)<eps)
                                {
                                    b=0.5;
                                }
                                else
                                {
                                    b=(atil-atd)/(1.0-atd);
                                }
                            }
                            else
                            {
                                b=0.0;
                            }
							bw=b;
							asw=((1.0-bw)*(phi[i-1][j][k]+phi_t[i-1][j][k])*0.5)+(0.5*bw*(phi[i][j][k]+phi_t[i][j][k]));
							ap=ap-(0.5*bw*uw*(-dt/dx));
						}
						else
						{
							ad=phi[i][j][k];
							aa=phi[i-1][j][k];
							au=phi[i+1][j][k];
							ae=0.5*(phi[i][j][k]+phi[i+1][j][k]);
							aw=0.5*(phi[i][j][k]+phi[i-1][j][k]);
							an=0.5*(phi[i][j][k]+phi[i][j+1][k]);
							as=0.5*(phi[i][j][k]+phi[i][j-1][k]);
							at=0.5*(phi[i][j][k]+phi[i][j][k+1]);
							ab=0.5*(phi[i][j][k]+phi[i][j][k-1]);
							diva=ae-aw;
							mdiva=sqrt(((ae-aw)/dx)*((ae-aw)/dx)+((an-as)/dy)*((an-as)/dy)+((at-ab)/dz)*((at-ab)/dz))*sqrt(dx*dx);

						//	courant(i,j,k);
						  if(ad!=1.0 || fabs(ad)>eps)
                            {
                                //****************tilda donor calculation***********************//
                                if(fabs(aa-au)<eps)
                                {
                                    atd=ad;
                                }
                                else
                                {
                                    atd=(ad-au)/(aa-au);
                                }
                                //****************SAISH******************************//
                                //****************convective boundedness scheme*****************//
                                if(atd>=0.0 && atd<=1.0/4.0)
                                    atfc=4.0*atd;
                                else if(atd>1.0/4.0  &&  atd<=1.0)
                                    atfc=1.0;
                                else
                                    atfc=atd;

                                //****************Hybrid HR******************************//
                                if(atd>=0.0 && atd<=1.0/6.0)
                                    atfu=11.0/6.0*atd;
                                else if(atd>1.0/6.0 && atd<=1.0/2.0)
                                     atfu=atd*(2.0-atd);			//CLAM
                                else if(atd>1.0/2.0 && atd<=3.0/4.0)
                                    atfu=atd+0.25;					// Fromm's
                                else if(atd>3.0/4.0 && atd<=1.0)
                                    atfu=1.0;
                                else
                                    atfu=atd;

                                theta=(fabs(diva/mdiva));
                                yy=fmin(pow(theta,2.0),1.0);
                                //****************tilda face calculation************************//
                                atil=yy*atfc+(1-yy)*atfu;
                                //******************************weighting factor calculation****************//
                                if(fabs(1.0-atd)<eps)
                                {
                                    b=0.5;
                                }
                                else
                                {
                                    b=(atil-atd)/(1.0-atd);
                                }
                            }
                            else
                            {
                                b=0.0;
                            }
							bw=b;
							asw=((1.0-bw)*(phi[i][j][k]+phi_t[i][j][k])*0.5)+(0.5*bw*(phi[i-1][j][k]+phi_t[i-1][j][k]));
							ap=ap-(0.5*(1.0-bw)*uw*(-dt/dx));
						}
						if(vn>0.0)
						{
							ad=phi[i][j][k];
							aa=phi[i][j+1][k];
							au=phi[i][j-1][k];
							ae=0.5*(phi[i][j][k]+phi[i+1][j][k]);
							aw=0.5*(phi[i][j][k]+phi[i-1][j][k]);
							an=0.5*(phi[i][j][k]+phi[i][j+1][k]);
							as=0.5*(phi[i][j][k]+phi[i][j-1][k]);
							at=0.5*(phi[i][j][k]+phi[i][j][k+1]);
							ab=0.5*(phi[i][j][k]+phi[i][j][k-1]);
							diva=an-as;
							mdiva=sqrt(((ae-aw)/dx)*((ae-aw)/dx)+((an-as)/dy)*((an-as)/dy)+((at-ab)/dz)*((at-ab)/dz))*sqrt(dy*dy);

						//	courant(i,j,k);
						  if(ad!=1.0 || fabs(ad)>eps)
                            {
                                //****************tilda donor calculation***********************//
                                if(fabs(aa-au)<eps)
                                {
                                    atd=ad;
                                }
                                else
                                {
                                    atd=(ad-au)/(aa-au);
                                }
                                //****************SAISH******************************//
                                //****************convective boundedness scheme*****************//
                                if(atd>=0.0 && atd<=1.0/4.0)
                                    atfc=4.0*atd;
                                else if(atd>1.0/4.0  &&  atd<=1.0)
                                    atfc=1.0;
                                else
                                    atfc=atd;

                                //****************Hybrid HR******************************//
                                if(atd>=0.0 && atd<=1.0/6.0)
                                    atfu=11.0/6.0*atd;
                                else if(atd>1.0/6.0 && atd<=1.0/2.0)
                                     atfu=atd*(2.0-atd);			//CLAM
                                else if(atd>1.0/2.0 && atd<=3.0/4.0)
                                    atfu=atd+0.25;					// Fromm's
                                else if(atd>3.0/4.0 && atd<=1.0)
                                    atfu=1.0;
                                else
                                    atfu=atd;

                                theta=(fabs(diva/mdiva));
                                yy=fmin(pow(theta,2.0),1.0);
                                //****************tilda face calculation************************//
                                atil=yy*atfc+(1-yy)*atfu;
                                //******************************weighting factor calculation****************//
                                if(fabs(1.0-atd)<eps)
                                {
                                    b=0.5;
                                }
                                else
                                {
                                    b=(atil-atd)/(1.0-atd);
                                }
                            }
                            else
                            {
                                b=0.0;
                            }
							bn=b;
							asn=((1.0-bn)*(phi[i][j][k]+phi_t[i][j][k])*0.5)+(0.5*bn*(phi[i][j+1][k]+phi_t[i][j+1][k]));
							ap=ap-(0.5*(1.0-bn)*vn*(dt/dy));
						}
						else
						{
							ad=phi[i][j+1][k];
							aa=phi[i][j][k];
							au=phi[i][j+2][k];
							ae=0.5*(phi[i][j+1][k]+phi[i+1][j+1][k]);
							aw=0.5*(phi[i][j+1][k]+phi[i-1][j+1][k]);
							an=0.5*(phi[i][j+1][k]+phi[i][j+2][k]);
							as=0.5*(phi[i][j+1][k]+phi[i][j][k]);
							at=0.5*(phi[i][j+1][k]+phi[i][j+1][k+1]);
							ab=0.5*(phi[i][j+1][k]+phi[i][j+1][k-1]);
							diva=an-as;
							mdiva=sqrt(((ae-aw)/dx)*((ae-aw)/dx)+((an-as)/dy)*((an-as)/dy)+((at-ab)/dz)*((at-ab)/dz))*sqrt(dy*dy);

						//	courant(i,j+1,k);
						  if(ad!=1.0 || fabs(ad)>eps)
                            {
                                //****************tilda donor calculation***********************//
                                if(fabs(aa-au)<eps)
                                {
                                    atd=ad;
                                }
                                else
                                {
                                    atd=(ad-au)/(aa-au);
                                }
                                //****************SAISH******************************//
                                //****************convective boundedness scheme*****************//
                                if(atd>=0.0 && atd<=1.0/4.0)
                                    atfc=4.0*atd;
                                else if(atd>1.0/4.0  &&  atd<=1.0)
                                    atfc=1.0;
                                else
                                    atfc=atd;

                                //****************Hybrid HR******************************//
                                if(atd>=0.0 && atd<=1.0/6.0)
                                    atfu=11.0/6.0*atd;
                                else if(atd>1.0/6.0 && atd<=1.0/2.0)
                                     atfu=atd*(2.0-atd);			//CLAM
                                else if(atd>1.0/2.0 && atd<=3.0/4.0)
                                    atfu=atd+0.25;					// Fromm's
                                else if(atd>3.0/4.0 && atd<=1.0)
                                    atfu=1.0;
                                else
                                    atfu=atd;

                                theta=(fabs(diva/mdiva));
                                yy=fmin(pow(theta,2.0),1.0);
                                //****************tilda face calculation************************//
                                atil=yy*atfc+(1-yy)*atfu;
                                //******************************weighting factor calculation****************//
                                if(fabs(1.0-atd)<eps)
                                {
                                    b=0.5;
                                }
                                else
                                {
                                    b=(atil-atd)/(1.0-atd);
                                }
                            }
                            else
                            {
                                b=0.0;
                            }
							bn=b;
							asn=((1.0-bn)*(phi[i][j+1][k]+phi_t[i][j+1][k])*0.5)+(0.5*bn*(phi[i][j][k]+phi_t[i][j][k]));
							ap=ap-(0.5*bn*vn*(dt/dy));
						}
						if(vs>0.0)
						{
							ad=phi[i][j-1][k];
							aa=phi[i][j][k];
							au=phi[i][j-2][k];
							ae=0.5*(phi[i][j-1][k]+phi[i+1][j-1][k]);
							aw=0.5*(phi[i][j-1][k]+phi[i-1][j-1][k]);
							an=0.5*(phi[i][j-1][k]+phi[i][j][k]);
							as=0.5*(phi[i][j-1][k]+phi[i][j-2][k]);
							at=0.5*(phi[i][j-1][k]+phi[i][j-1][k+1]);
							ab=0.5*(phi[i][j-1][k]+phi[i][j-1][k-1]);
							diva=an-as;
							mdiva=sqrt(((ae-aw)/dx)*((ae-aw)/dx)+((an-as)/dy)*((an-as)/dy)+((at-ab)/dz)*((at-ab)/dz))*sqrt(dy*dy);


						//	courant(i,j-1,k);
						  if(ad!=1.0 || fabs(ad)>eps)
                            {
                                //****************tilda donor calculation***********************//
                                if(fabs(aa-au)<eps)
                                {
                                    atd=ad;
                                }
                                else
                                {
                                    atd=(ad-au)/(aa-au);
                                }
                                //****************SAISH******************************//
                                //****************convective boundedness scheme*****************//
                                if(atd>=0.0 && atd<=1.0/4.0)
                                    atfc=4.0*atd;
                                else if(atd>1.0/4.0  &&  atd<=1.0)
                                    atfc=1.0;
                                else
                                    atfc=atd;

                                //****************Hybrid HR******************************//
                                if(atd>=0.0 && atd<=1.0/6.0)
                                    atfu=11.0/6.0*atd;
                                else if(atd>1.0/6.0 && atd<=1.0/2.0)
                                     atfu=atd*(2.0-atd);			//CLAM
                                else if(atd>1.0/2.0 && atd<=3.0/4.0)
                                    atfu=atd+0.25;					// Fromm's
                                else if(atd>3.0/4.0 && atd<=1.0)
                                    atfu=1.0;
                                else
                                    atfu=atd;

                                theta=(fabs(diva/mdiva));
                                yy=fmin(pow(theta,2.0),1.0);
                                //****************tilda face calculation************************//
                                atil=yy*atfc+(1-yy)*atfu;
                                //******************************weighting factor calculation****************//
                                if(fabs(1.0-atd)<eps)
                                {
                                    b=0.5;
                                }
                                else
                                {
                                    b=(atil-atd)/(1.0-atd);
                                }
                            }
                            else
                            {
                                b=0.0;
                            }
							bs=b;
							ass=((1.0-bs)*(phi[i][j-1][k]+phi_t[i][j-1][k])*0.5)+(0.5*bs*(phi[i][j][k]+phi_t[i][j][k]));
							ap=ap-(0.5*bs*vs*(-dt/dy));
						}
						else
						{
							ad=phi[i][j][k];
							aa=phi[i][j-1][k];
							au=phi[i][j+1][k];
							ae=0.5*(phi[i][j][k]+phi[i+1][j][k]);
							aw=0.5*(phi[i][j][k]+phi[i-1][j][k]);
							an=0.5*(phi[i][j][k]+phi[i][j+1][k]);
							as=0.5*(phi[i][j][k]+phi[i][j-1][k]);
							at=0.5*(phi[i][j][k]+phi[i][j][k+1]);
							ab=0.5*(phi[i][j][k]+phi[i][j][k-1]);
							diva=an-as;
							mdiva=sqrt(((ae-aw)/dx)*((ae-aw)/dx)+((an-as)/dy)*((an-as)/dy)+((at-ab)/dz)*((at-ab)/dz))*sqrt(dy*dy);

						//	courant(i,j,k);
						  if(ad!=1.0 || fabs(ad)>eps)
                            {
                                //****************tilda donor calculation***********************//
                                if(fabs(aa-au)<eps)
                                {
                                    atd=ad;
                                }
                                else
                                {
                                    atd=(ad-au)/(aa-au);
                                }
                                //****************SAISH******************************//
                                //****************convective boundedness scheme*****************//
                                if(atd>=0.0 && atd<=1.0/4.0)
                                    atfc=4.0*atd;
                                else if(atd>1.0/4.0  &&  atd<=1.0)
                                    atfc=1.0;
                                else
                                    atfc=atd;

                                //****************Hybrid HR******************************//
                                if(atd>=0.0 && atd<=1.0/6.0)
                                    atfu=11.0/6.0*atd;
                                else if(atd>1.0/6.0 && atd<=1.0/2.0)
                                     atfu=atd*(2.0-atd);			//CLAM
                                else if(atd>1.0/2.0 && atd<=3.0/4.0)
                                    atfu=atd+0.25;					// Fromm's
                                else if(atd>3.0/4.0 && atd<=1.0)
                                    atfu=1.0;
                                else
                                    atfu=atd;

                                theta=(fabs(diva/mdiva));
                                yy=fmin(pow(theta,2.0),1.0);
                                //****************tilda face calculation************************//
                                atil=yy*atfc+(1-yy)*atfu;
                                //******************************weighting factor calculation****************//
                                if(fabs(1.0-atd)<eps)
                                {
                                    b=0.5;
                                }
                                else
                                {
                                    b=(atil-atd)/(1.0-atd);
                                }
                            }
                            else
                            {
                                b=0.0;
                            }
							bs=b;
							ass=((1.0-bs)*(phi[i][j][k]+phi_t[i][j][k])*0.5)+(0.5*bs*(phi[i][j-1][k]+phi_t[i][j-1][k]));
							ap=ap-(0.5*(1.0-bs)*vs*(-dt/dy));
						}
						if(wt>0.0)
						{
							ad=phi[i][j][k];
							aa=phi[i][j][k+1];
							au=phi[i][j][k-1];
							ae=0.5*(phi[i][j][k]+phi[i+1][j][k]);
							aw=0.5*(phi[i][j][k]+phi[i-1][j][k]);
							an=0.5*(phi[i][j][k]+phi[i][j+1][k]);
							as=0.5*(phi[i][j][k]+phi[i][j-1][k]);
							at=0.5*(phi[i][j][k]+phi[i][j][k+1]);
							ab=0.5*(phi[i][j][k]+phi[i][j][k-1]);
							diva=at-ab;
							mdiva=sqrt(((ae-aw)/dx)*((ae-aw)/dx)+((an-as)/dy)*((an-as)/dy)+((at-ab)/dz)*((at-ab)/dz))*sqrt(dz*dz);


						//	courant(i,j,k);
						  if(ad!=1.0 || fabs(ad)>eps)
                            {
                                //****************tilda donor calculation***********************//
                                if(fabs(aa-au)<eps)
                                {
                                    atd=ad;
                                }
                                else
                                {
                                    atd=(ad-au)/(aa-au);
                                }
                                //****************SAISH******************************//
                                //****************convective boundedness scheme*****************//
                                if(atd>=0.0 && atd<=1.0/4.0)
                                    atfc=4.0*atd;
                                else if(atd>1.0/4.0  &&  atd<=1.0)
                                    atfc=1.0;
                                else
                                    atfc=atd;

                                //****************Hybrid HR******************************//
                                if(atd>=0.0 && atd<=1.0/6.0)
                                    atfu=11.0/6.0*atd;
                                else if(atd>1.0/6.0 && atd<=1.0/2.0)
                                     atfu=atd*(2.0-atd);			//CLAM
                                else if(atd>1.0/2.0 && atd<=3.0/4.0)
                                    atfu=atd+0.25;					// Fromm's
                                else if(atd>3.0/4.0 && atd<=1.0)
                                    atfu=1.0;
                                else
                                    atfu=atd;

                                theta=(fabs(diva/mdiva));
                                yy=fmin(pow(theta,2.0),1.0);
                                //****************tilda face calculation************************//
                                atil=yy*atfc+(1-yy)*atfu;
                                //******************************weighting factor calculation****************//
                                if(fabs(1.0-atd)<eps)
                                {
                                    b=0.5;
                                }
                                else
                                {
                                    b=(atil-atd)/(1.0-atd);
                                }
                            }
                            else
                            {
                                b=0.0;
                            }
							bt=b;
							ast=((1.0-bt)*(phi[i][j][k]+phi_t[i][j][k])*0.5)+(0.5*bt*(phi[i][j][k+1]+phi_t[i][j][k+1]));
							ap=ap-(0.5*(1.0-bt)*wt*(dt/dz));
						}
						else
						{
							ad=phi[i][j][k+1];
							aa=phi[i][j][k];
							au=phi[i][j][k+2];
							ae=0.5*(phi[i][j][k+1]+phi[i+1][j][k+1]);
							aw=0.5*(phi[i][j][k+1]+phi[i-1][j][k+1]);
							an=0.5*(phi[i][j][k+1]+phi[i][j+1][k+1]);
							as=0.5*(phi[i][j][k+1]+phi[i][j-1][k+1]);
							at=0.5*(phi[i][j][k+1]+phi[i][j][k+2]);
							ab=0.5*(phi[i][j][k+1]+phi[i][j][k]);
							diva=at-ab;
							mdiva=sqrt(((ae-aw)/dx)*((ae-aw)/dx)+((an-as)/dy)*((an-as)/dy)+((at-ab)/dz)*((at-ab)/dz))*sqrt(dz*dz);

						//	courant(i,j,k+1);
						  if(ad!=1.0 || fabs(ad)>eps)
                            {
                                //****************tilda donor calculation***********************//
                                if(fabs(aa-au)<eps)
                                {
                                    atd=ad;
                                }
                                else
                                {
                                    atd=(ad-au)/(aa-au);
                                }
                                //****************SAISH******************************//
                                //****************convective boundedness scheme*****************//
                                if(atd>=0.0 && atd<=1.0/4.0)
                                    atfc=4.0*atd;
                                else if(atd>1.0/4.0  &&  atd<=1.0)
                                    atfc=1.0;
                                else
                                    atfc=atd;

                                //****************Hybrid HR******************************//
                                if(atd>=0.0 && atd<=1.0/6.0)
                                    atfu=11.0/6.0*atd;
                                else if(atd>1.0/6.0 && atd<=1.0/2.0)
                                     atfu=atd*(2.0-atd);			//CLAM
                                else if(atd>1.0/2.0 && atd<=3.0/4.0)
                                    atfu=atd+0.25;					// Fromm's
                                else if(atd>3.0/4.0 && atd<=1.0)
                                    atfu=1.0;
                                else
                                    atfu=atd;

                                theta=(fabs(diva/mdiva));
                                yy=fmin(pow(theta,2.0),1.0);
                                //****************tilda face calculation************************//
                                atil=yy*atfc+(1-yy)*atfu;
                                //******************************weighting factor calculation****************//
                                if(fabs(1.0-atd)<eps)
                                {
                                    b=0.5;
                                }
                                else
                                {
                                    b=(atil-atd)/(1.0-atd);
                                }
                            }
                            else
                            {
                                b=0.0;
                            }
							bt=b;
							ast=((1.0-bt)*(phi[i][j][k+1]+phi_t[i][j][k+1])*0.5)+(0.5*bt*(phi[i][j][k]+phi_t[i][j][k]));
							ap=ap-(0.5*bt*wt*(dt/dz));
						}
						if(wb>0.0)
						{
							ad=phi[i][j][k-1];
							aa=phi[i][j][k];
							au=phi[i][j][k-2];
							ae=0.5*(phi[i][j][k-1]+phi[i+1][j][k-1]);
							aw=0.5*(phi[i][j][k-1]+phi[i-1][j][k-1]);
							an=0.5*(phi[i][j][k-1]+phi[i][j+1][k-1]);
							as=0.5*(phi[i][j][k-1]+phi[i][j-1][k-1]);
							at=0.5*(phi[i][j][k-1]+phi[i][j][k]);
							ab=0.5*(phi[i][j][k-1]+phi[i][j][k-2]);
							diva=at-ab;
							mdiva=sqrt(((ae-aw)/dx)*((ae-aw)/dx)+((an-as)/dy)*((an-as)/dy)+((at-ab)/dz)*((at-ab)/dz))*sqrt(dz*dz);


						//	courant(i,j,k-1);
						  if(ad!=1.0 || fabs(ad)>eps)
                            {
                                //****************tilda donor calculation***********************//
                                if(fabs(aa-au)<eps)
                                {
                                    atd=ad;
                                }
                                else
                                {
                                    atd=(ad-au)/(aa-au);
                                }
                                //****************SAISH******************************//
                                //****************convective boundedness scheme*****************//
                                if(atd>=0.0 && atd<=1.0/4.0)
                                    atfc=4.0*atd;
                                else if(atd>1.0/4.0  &&  atd<=1.0)
                                    atfc=1.0;
                                else
                                    atfc=atd;

                                //****************Hybrid HR******************************//
                                if(atd>=0.0 && atd<=1.0/6.0)
                                    atfu=11.0/6.0*atd;
                                else if(atd>1.0/6.0 && atd<=1.0/2.0)
                                     atfu=atd*(2.0-atd);			//CLAM
                                else if(atd>1.0/2.0 && atd<=3.0/4.0)
                                    atfu=atd+0.25;					// Fromm's
                                else if(atd>3.0/4.0 && atd<=1.0)
                                    atfu=1.0;
                                else
                                    atfu=atd;

                                theta=(fabs(diva/mdiva));
                                yy=fmin(pow(theta,2.0),1.0);
                                //****************tilda face calculation************************//
                                atil=yy*atfc+(1-yy)*atfu;
                                //******************************weighting factor calculation****************//
                                if(fabs(1.0-atd)<eps)
                                {
                                    b=0.5;
                                }
                                else
                                {
                                    b=(atil-atd)/(1.0-atd);
                                }
                            }
                            else
                            {
                                b=0.0;
                            }
							bb=b;
							asb=((1.0-bb)*(phi[i][j][k-1]+phi_t[i][j][k-1])*0.5)+(0.5*bb*(phi[i][j][k]+phi_t[i][j][k]));
							ap=ap-(0.5*bb*wb*(-dt/dz));
						}
						else
						{
							ad=phi[i][j][k];
							aa=phi[i][j][k-1];
							au=phi[i][j][k+1];
							ae=0.5*(phi[i][j][k]+phi[i+1][j][k]);
							aw=0.5*(phi[i][j][k]+phi[i-1][j][k]);
							an=0.5*(phi[i][j][k]+phi[i][j+1][k]);
							as=0.5*(phi[i][j][k]+phi[i][j-1][k]);
							at=0.5*(phi[i][j][k]+phi[i][j][k+1]);
							ab=0.5*(phi[i][j][k]+phi[i][j][k-1]);
							diva=at-ab;
							mdiva=sqrt(((ae-aw)/dx)*((ae-aw)/dx)+((an-as)/dy)*((an-as)/dy)+((at-ab)/dz)*((at-ab)/dz))*sqrt(dz*dz);

						//	courant(i,j,k);
						  if(ad!=1.0 || fabs(ad)>eps)
                            {
                                //****************tilda donor calculation***********************//
                                if(fabs(aa-au)<eps)
                                {
                                    atd=ad;
                                }
                                else
                                {
                                    atd=(ad-au)/(aa-au);
                                }
                                //****************SAISH******************************//
                                //****************convective boundedness scheme*****************//
                                if(atd>=0.0 && atd<=1.0/4.0)
                                    atfc=4.0*atd;
                                else if(atd>1.0/4.0  &&  atd<=1.0)
                                    atfc=1.0;
                                else
                                    atfc=atd;

                                //****************Hybrid HR******************************//
                                if(atd>=0.0 && atd<=1.0/6.0)
                                    atfu=11.0/6.0*atd;
                                else if(atd>1.0/6.0 && atd<=1.0/2.0)
                                     atfu=atd*(2.0-atd);			//CLAM
                                else if(atd>1.0/2.0 && atd<=3.0/4.0)
                                    atfu=atd+0.25;					// Fromm's
                                else if(atd>3.0/4.0 && atd<=1.0)
                                    atfu=1.0;
                                else
                                    atfu=atd;

                                theta=(fabs(diva/mdiva));
                                yy=fmin(pow(theta,2.0),1.0);
                                //****************tilda face calculation************************//
                                atil=yy*atfc+(1-yy)*atfu;
                                //******************************weighting factor calculation****************//
                                if(fabs(1.0-atd)<eps)
                                {
                                    b=0.5;
                                }
                                else
                                {
                                    b=(atil-atd)/(1.0-atd);
                                }
                            }
                            else
                            {
                                b=0.0;
                            }
							bb=b;
							asb=((1.0-bb)*(phi[i][j][k]+phi_t[i][j][k])*0.5)+(0.5*bb*(phi[i][j][k-1]+phi_t[i][j][k-1]));
							ap=ap-(0.5*(1.0-bb)*wb*(-dt/dz));
						}

						ap=-(ap-1.0);
						res=-((ase*ue*dy*dz)-(asw*uw*dy*dz)+(asn*vn*dx*dz)-(ass*vs*dx*dz)+(ast*wt*dx*dy)-(asb*wb*dx*dy))*(dt/(dx*dy*dz))+phi[i][j][k]-phi_t[i][j][k];
						phi_t[i][j][k]=phi_t[i][j][k]+(res*alpha/ap);
						sum=sum+(res*res);
					}
				}
			}
			rms=sqrt(sum/(nx*ny*nz));
			l=l+1;
		}while(rms>=1e-15);

        redistribution();

//**************************new value assigned to old value**************************/
#pragma omp parallel for schedule(dynamic) private(i,j,k)
		for(i=img;i<nx+img;i++)
	    	{
			for(j=img;j<ny+img;j++)
				{
			    		for(k=img;k<nx+img;k++)
					{
						phi[i][j][k]=phi_t[i][j][k];
					}
				}
	    	}
	//****************************VOF ends here**************************************//
}

//Redistribution of volume fraction values
void redistribution()
{
	do
	{
		forlorn=0;
//#pragma omp parallel for collapse(3) schedule(static) private(I,min,max,al,bl,cl)
	for(i=img;i<nx+img;i++)
    	{
        for(j=img;j<ny+img;j++)
        	{
		        for(k=img;k<nz+img;k++)
				{
						if(phi_t[i][j][k]<-eps)//undershoot
                        {
                                I=1;
								do
								{
									min=10.0;
									for(s=i-I;s<=I+i;s++)
									{
										for(m=j-I;m<=j+I;m++)
										{
											for(n=k-I;n<=k+I;n++)
											{
												if((s==i&&m==j&&n==k)||(s<1)||(m<1)||(n<1)||(s>nx)||(m>ny)||(n>nz))
													continue;
												else
													{
														if(phi_t[s][m][n]<=min&&phi_t[s][m][n]!=0.0)
																{
																	min=phi_t[s][m][n];
																	al=s;
																	bl=m;
																	cl=n;
																}
													}
											}
										}
									}
												if(phi_t[al][bl][cl]==0.0)
													{
														I=I+1;
													}
												phi_temp=phi_t[al][bl][cl]+phi_t[i][j][k];
												phi_t[al][bl][cl]=fmax(phi_temp,0.0);
												phi_t[i][j][k]=fmin(phi_temp,0.0);
							}while(phi_t[i][j][k]<-eps);
						}else if (phi_t[i][j][k]>1.0+eps)//Overshoot
							{
								I=1;
										do
										{
											max=-10.0;
											for(s=i-I;s<=I+i;s++)
											{
												for(m=j-I;m<=j+I;m++)
												{
													for(n=k-I;n<=k+I;n++)
													{
														if((s==i&&m==j&&n==k)||(s<1)||(m<1)||(n<1)||(s>nx)||(m>ny)||(n>nz))
														continue;
														else
															{
																if(phi_t[s][m][n]>=max&&phi_t[s][m][n]!=1.0)
																	{
																		max=phi_t[s][m][n];
																		al=s;
																		bl=m;
																		cl=n;
																	}
															}
													}
												}
											}
											if(phi_t[al][bl][cl]==1.0)
											{
												I=I+1;
											}
											phi_temp=phi_t[al][bl][cl]+phi_t[i][j][k]-1.0;
											phi_t[al][bl][cl]=fmin(phi_temp,1.0);
											phi_t[i][j][k]=fmax(phi_temp,1.0);


										}while(phi_t[i][j][k]>1.0+eps);
							}else
                            continue;
					}
				}
			}

#pragma omp parallel for schedule(dynamic) private(i,j,k)
for(i=img;i<nx+img;i++)
    	{
        for(j=img;j<ny+img;j++)
        	{
		        for(k=img;k<nz+img;k++)
				{
							if((phi_t[i][j][k]>1.0+eps)||(phi_t[i][j][k]<-eps))
							{
							forlorn=2;
							break;
							}else
								continue;

					}
				}
			}
	}while(forlorn>0);
}


void normals()
{
#pragma omp parallel for schedule(dynamic) private(i,j,k,nx_temp,ny_temp,nz_temp)
for(i=img;i<nx+img;i++)
	{
	for(j=img;j<ny+img;j++)
		{
		for(k=img;k<nz+img;k++)
			{

				Nx[i][j][k]=(phi[i+1][j+1][k+1]+2.0*phi[i+1][j+1][k]+phi[i+1][j+1][k-1]+2.0*phi[i+1][j][k+1]+4.0*phi[i+1][j][k]+2.0*phi[i+1][j][k-1]+phi[i+1][j-1][k+1]+2.0*phi[i+1][j-1][k]+phi[i+1][j-1][k-1]-phi[i-1][j+1][k+1]-2.0*phi[i-1][j+1][k]-phi[i-1][j+1][k-1]-2.0*phi[i-1][j][k+1]-4.0*phi[i-1][j][k]-2.0*phi[i-1][j][k-1]-phi[i-1][j-1][k+1]-2.0*phi[i-1][j-1][k]-phi[i-1][j-1][k-1])/(32.0*dx);
				Ny[i][j][k]=(phi[i+1][j+1][k+1]+2.0*phi[i+1][j+1][k]+phi[i+1][j+1][k-1]+2.0*phi[i][j+1][k+1]+4.0*phi[i][j+1][k]+2.0*phi[i][j+1][k-1]+phi[i-1][j+1][k+1]+2.0*phi[i-1][j+1][k]+phi[i-1][j+1][k-1]-phi[i+1][j-1][k+1]-2.0*phi[i+1][j-1][k]-phi[i+1][j-1][k-1]-2.0*phi[i][j-1][k+1]-4.0*phi[i][j-1][k]-2.0*phi[i][j-1][k-1]-phi[i-1][j-1][k+1]-2.0*phi[i-1][j-1][k]-phi[i-1][j-1][k-1])/(32.0*dy);
				Nz[i][j][k]=(phi[i+1][j+1][k+1]+2.0*phi[i+1][j][k+1]+phi[i+1][j-1][k+1]+2.0*phi[i][j+1][k+1]+4.0*phi[i][j][k+1]+2.0*phi[i][j-1][k+1]+phi[i-1][j+1][k+1]+2.0*phi[i-1][j][k+1]+phi[i-1][j-1][k+1]-phi[i+1][j+1][k-1]-2.0*phi[i+1][j][k-1]-phi[i+1][j-1][k-1]-2.0*phi[i][j+1][k-1]-4.0*phi[i][j][k-1]-2.0*phi[i][j-1][k-1]-phi[i-1][j+1][k-1]-2.0*phi[i-1][j][k-1]-phi[i-1][j-1][k-1])/(32.0*dz);

				fabs(Nx[i][j][k])<eps?Nx[i][j][k]=0.0:0;
				fabs(Ny[i][j][k])<eps?Ny[i][j][k]=0.0:0;
				fabs(Nz[i][j][k])<eps?Nz[i][j][k]=0.0:0;

				nx_temp=Nx[i][j][k];
				ny_temp=Ny[i][j][k];
				nz_temp=Nz[i][j][k];

				if(sqrt(nx_temp*nx_temp+ny_temp*ny_temp+nz_temp*nz_temp)!=0.0)
				{
					Nx[i][j][k]=Nx[i][j][k]/sqrt(nx_temp*nx_temp+ny_temp*ny_temp+nz_temp*nz_temp);
					Ny[i][j][k]=Ny[i][j][k]/sqrt(nx_temp*nx_temp+ny_temp*ny_temp+nz_temp*nz_temp);
					Nz[i][j][k]=Nz[i][j][k]/sqrt(nx_temp*nx_temp+ny_temp*ny_temp+nz_temp*nz_temp);
				}

			}
		}
	}
}
void HF()
{
tup=3;
tlow=-3;
#pragma omp parallel for schedule(dynamic) private(i,j,k,h1,h2,h3,h4,h5,h6,h7,h8,h9,ha,hb,haa,hbb,hab)
for(i=img;i<nx+img;i++)
	{
	for(j=img;j<ny+img;j++)
		{
		for(k=img;k<nz+img;k++)
			{
                rt=-1.0;
                h1=0.0; 	h2=0.0;		h3=0.0; 	h4=0.0;		h5=0.0; 	h6=0.0;		h7=0.0; 	h8=0.0; 	h9=0.0;
				if((fabs(Nx[i][j][k])>=fabs(Ny[i][j][k])) && (fabs(Nx[i][j][k])>=fabs(Nz[i][j][k]))) // Nx is maximum
					{

							for(to=i+tlow;to<=i+tup;to++)
								{
										h1+=phi[to][j+1][k+1]*dx;
										h2+=phi[to][j+1][k]*dx;
										h3+=phi[to][j+1][k-1]*dx;
										h4+=phi[to][j][k+1]*dx;
										h5+=phi[to][j][k]*dx;
										h6+=phi[to][j][k-1]*dx;
										h7+=phi[to][j-1][k+1]*dx;
										h8+=phi[to][j-1][k]*dx;
										h9+=phi[to][j-1][k-1]*dx;
								}
								haa=(h4-2.0*h5+h6)/(dz*dz);
								hbb=(h2-2.0*h5+h8)/(dy*dy);
								ha=((h4-h6))/(2.0*dz);
								hb=((h2-h8))/(2.0*dy);
								hab=(h1-h3-h7+h9)/(4.0*dy*dz);


					h1=0.0;  h2=0.0;	h3=0.0; 	h4=0.0;		h5=0.0;	 	h6=0.0;		h7=0.0;	 	h8=0.0;	 	h9=0.0;
                    Kurv[i][j][k]=rt*(haa+hbb+haa*hb*hb+hbb*ha*ha-2.0*hab*ha*hb)/pow((1.0+ha*ha+hb*hb),1.5);
					}else if((fabs(Ny[i][j][k])>=fabs(Nx[i][j][k])) && (fabs(Ny[i][j][k])>=fabs(Nz[i][j][k]))) // Ny is maximum
						{

								for(to=j+tlow;to<=j+tup;to++)
								{
									h1+=phi[i+1][to][k+1]*dy;
									h2+=phi[i+1][to][k]*dy;
									h3+=phi[i+1][to][k-1]*dy;
									h4+=phi[i][to][k+1]*dy;
									h5+=phi[i][to][k]*dy;
									h6+=phi[i][to][k-1]*dy;
									h7+=phi[i-1][to][k+1]*dy;
									h8+=phi[i-1][to][k]*dy;
									h9+=phi[i-1][to][k-1]*dy;
								}


								haa=(h4-2.0*h5+h6)/(dz*dz);
								hbb=(h2-2.0*h5+h8)/(dx*dx);
								ha=((h4-h6))/(2.0*dz);
								hb=((h2-h8))/(2.0*dx);
								hab=(h1-h3-h7+h9)/(4.0*dx*dz);
							h1=0.0; 	h2=0.0;		h3=0.0; 	h4=0.0;		h5=0.0; 	h6=0.0;		h7=0.0; 	h8=0.0; 	h9=0.0;
							Kurv[i][j][k]=rt*(haa+hbb+haa*hb*hb+hbb*ha*ha-2.0*hab*ha*hb)/pow((1.0+ha*ha+hb*hb),1.5);
						}else if((fabs(Nz[i][j][k])>=fabs(Ny[i][j][k])) && (fabs(Nz[i][j][k])>=fabs(Nx[i][j][k])))	// Nz is maximum
							{

									for(to=k+tlow;to<=k+tup;to++)
									{
										h1+=phi[i+1][j+1][to]*dz;
										h2+=phi[i+1][j][to]*dz;
										h3+=phi[i+1][j-1][to]*dz;
										h4+=phi[i][j+1][to]*dz;
										h5+=phi[i][j][to]*dz;
										h6+=phi[i][j-1][to]*dz;
										h7+=phi[i-1][j+1][to]*dz;
										h8+=phi[i-1][j][to]*dz;
										h9+=phi[i-1][j-1][to]*dz;
									}

									haa=(h4-2.0*h5+h6)/(dy*dy);
									hbb=(h2-2.0*h5+h8)/(dx*dx);
									ha=((h4-h6))/(2.0*dy);
									hb=((h2-h8))/(2.0*dx);
									hab=(h1-h3-h7+h9)/(4.0*dy*dx);
									h1=0.0; 	h2=0.0;		h3=0.0; 	h4=0.0;		h5=0.0; 	h6=0.0;		h7=0.0; 	h8=0.0; 	h9=0.0;
									Kurv[i][j][k]=rt*(haa+hbb+haa*hb*hb+hbb*ha*ha-2.0*hab*ha*hb)/pow((1.0+ha*ha+hb*hb),1.5);
							}

			}
		}
	}
}

void cvr()
{
		//				Cell Vertices
				// 		4(top left)		3(top right)
 				//		1(bottom left)	 	2(bottom right)
	xv[1]=xp[i];		yv[1]=yp[j];	zv[1]=zp[k];
	xv[2]=xp[i+1];		yv[2]=yp[j];	zv[2]=zp[k];
	xv[3]=xp[i+1];		yv[3]=yp[j+1];	zv[3]=zp[k];
	xv[4]=xp[i];		yv[4]=yp[j+1];	zv[4]=zp[k];

	xv[5]=xp[i];		yv[5]=yp[j];	zv[5]=zp[k+1];
	xv[6]=xp[i+1];		yv[6]=yp[j];	zv[6]=zp[k+1];
	xv[7]=xp[i+1];		yv[7]=yp[j+1];	zv[7]=zp[k+1];
	xv[8]=xp[i];		yv[8]=yp[j+1];	zv[8]=zp[k+1];
}

