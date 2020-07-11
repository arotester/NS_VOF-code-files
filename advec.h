void Superbee()
{

       //TVD face velocities
		if(ue>0.0)
			{
                r_ue=(u[i][j][k]-u[i-1][j][k])/(u[i+1][j][k]-u[i][j][k]);
                r_ve=(v[i][j][k]-v[i-1][j][k])/(v[i+1][j][k]-v[i][j][k]);
                r_we=(w[i][j][k]-w[i-1][j][k])/(w[i+1][j][k]-w[i][j][k]);

                a_ue=fmin(1.0,2.0*r_ue);
                b_ue=fmin(2.0,r_ue);
                c_ue=fmax(a_ue,b_ue);
                limiter_ue=fmax(0.0,c_ue);

                a_ve=fmin(1.0,2.0*r_ve);
                b_ve=fmin(2.0,r_ve);
                c_ve=fmax(a_ve,b_ve);
                limiter_ve=fmax(0.0,c_ve);

                a_we=fmin(1.0,2.0*r_we);
                b_we=fmin(2.0,r_we);
                c_we=fmax(a_we,b_we);
                limiter_we=fmax(0.0,c_we);

				uea=u[i][j][k]+0.5*limiter_ue*(u[i+1][j][k]-u[i][j][k]);
				vea=v[i][j][k]+0.5*limiter_ve*(v[i+1][j][k]-v[i][j][k]);
				wea=w[i][j][k]+0.5*limiter_we*(w[i+1][j][k]-w[i][j][k]);
			}else
				{
                    r_ue=(u[i+1][j][k]-u[i][j][k])/(u[i][j][k]-u[i-1][j][k]);
                    r_ve=(v[i+1][j][k]-v[i][j][k])/(v[i][j][k]-v[i-1][j][k]);
                    r_we=(w[i+1][j][k]-w[i][j][k])/(w[i][j][k]-w[i-1][j][k]);

                    a_ue=fmin(1.0,2.0*r_ue);
                    b_ue=fmin(2.0,r_ue);
                    c_ue=fmax(a_ue,b_ue);
                    limiter_ue=fmax(0.0,c_ue);

                    a_ve=fmin(1.0,2.0*r_ve);
                    b_ve=fmin(2.0,r_ve);
                    c_ve=fmax(a_ve,b_ve);
                    limiter_ve=fmax(0.0,c_ve);

                    a_we=fmin(1.0,2.0*r_we);
                    b_we=fmin(2.0,r_we);
                    c_we=fmax(a_we,b_we);
                    limiter_we=fmax(0.0,c_we);

					uea=u[i+1][j][k]+0.5*limiter_ue*(u[i][j][k]-u[i+1][j][k]);
					vea=v[i+1][j][k]+0.5*limiter_ve*(v[i][j][k]-v[i+1][j][k]);
					wea=w[i+1][j][k]+0.5*limiter_we*(w[i][j][k]-w[i+1][j][k]);
				}
		if(uw>0.0)
			{

                r_uw=(u[i-1][j][k]-u[i-2][j][k])/(u[i][j][k]-u[i-1][j][k]);
                r_vw=(v[i-1][j][k]-v[i-2][j][k])/(v[i][j][k]-v[i-1][j][k]);
                r_ww=(w[i-1][j][k]-w[i-2][j][k])/(w[i][j][k]-w[i-1][j][k]);

                a_uw=fmin(1.0,2.0*r_uw);
                b_uw=fmin(2.0,r_uw);
                c_uw=fmax(a_uw,b_uw);
                limiter_uw=fmax(0.0,c_uw);

                a_vw=fmin(1.0,2.0*r_vw);
                b_vw=fmin(2.0,r_vw);
                c_vw=fmax(a_vw,b_vw);
                limiter_vw=fmax(0.0,c_vw);

                a_ww=fmin(1.0,2.0*r_ww);
                b_ww=fmin(2.0,r_ww);
                c_ww=fmax(a_ww,b_ww);
                limiter_ww=fmax(0.0,c_ww);

				uwa=u[i-1][j][k]+0.5*limiter_uw*(u[i][j][k]-u[i-1][j][k]);
				vwa=v[i-1][j][k]+0.5*limiter_vw*(v[i][j][k]-v[i-1][j][k]);
				wwa=w[i-1][j][k]+0.5*limiter_ww*(w[i][j][k]-w[i-1][j][k]);

			}else
				{
                    r_uw=(u[i][j][k]-u[i-1][j][k])/(u[i-1][j][k]-u[i-2][j][k]);
                    r_vw=(v[i][j][k]-v[i-1][j][k])/(v[i-1][j][k]-v[i-2][j][k]);
                    r_ww=(w[i][j][k]-w[i-1][j][k])/(w[i-1][j][k]-w[i-2][j][k]);

                    a_uw=fmin(1.0,2.0*r_uw);
                    b_uw=fmin(2.0,r_uw);
                    c_uw=fmax(a_uw,b_uw);
                    limiter_uw=fmax(0.0,c_uw);

                    a_vw=fmin(1.0,2.0*r_vw);
                    b_vw=fmin(2.0,r_vw);
                    c_vw=fmax(a_vw,b_vw);
                    limiter_vw=fmax(0.0,c_vw);

                    a_ww=fmin(1.0,2.0*r_ww);
                    b_ww=fmin(2.0,r_ww);
                    c_ww=fmax(a_ww,b_ww);
                    limiter_ww=fmax(0.0,c_ww);

					uwa=u[i][j][k]+0.5*limiter_uw*(u[i-1][j][k]-u[i][j][k]);
					vwa=v[i][j][k]+0.5*limiter_vw*(v[i-1][j][k]-v[i][j][k]);
                    wwa=w[i][j][k]+0.5*limiter_ww*(w[i-1][j][k]-w[i][j][k]);
				}
		if(vn>0.0)
			{
                r_un=(u[i][j][k]-u[i][j-1][k])/(u[i][j+1][k]-u[i][j][k]);
                r_vn=(v[i][j][k]-v[i][j-1][k])/(v[i][j+1][k]-v[i][j][k]);
                r_wn=(w[i][j][k]-w[i][j-1][k])/(w[i][j+1][k]-w[i][j][k]);

                a_un=fmin(1.0,2.0*r_un);
                b_un=fmin(2.0,r_un);
                c_un=fmax(a_un,b_un);
                limiter_un=fmax(0.0,c_un);

                a_vn=fmin(1.0,2.0*r_vn);
                b_vn=fmin(2.0,r_vn);
                c_vn=fmax(a_vn,b_vn);
                limiter_vn=fmax(0.0,c_vn);

                a_wn=fmin(1.0,2.0*r_wn);
                b_wn=fmin(2.0,r_wn);
                c_wn=fmax(a_wn,b_wn);
                limiter_wn=fmax(0.0,c_wn);


				una=u[i][j][k]+0.5*limiter_un*(u[i][j+1][k]-u[i][j][k]);
				vna=v[i][j][k]+0.5*limiter_vn*(v[i][j+1][k]-v[i][j][k]);
                wna=w[i][j][k]+0.5*limiter_wn*(w[i][j+1][k]-w[i][j][k]);
			}else
				{
                    r_un=(u[i][j+1][k]-u[i][j][k])/(u[i][j][k]-u[i][j-1][k]);
                    r_vn=(v[i][j+1][k]-v[i][j][k])/(v[i][j][k]-v[i][j-1][k]);
                    r_wn=(w[i][j+1][k]-w[i][j][k])/(w[i][j][k]-w[i][j-1][k]);

                    a_un=fmin(1.0,2.0*r_un);
                    b_un=fmin(2.0,r_un);
                    c_un=fmax(a_un,b_un);
                    limiter_un=fmax(0.0,c_un);

                    a_vn=fmin(1.0,2.0*r_vn);
                    b_vn=fmin(2.0,r_vn);
                    c_vn=fmax(a_vn,b_vn);
                    limiter_vn=fmax(0.0,c_vn);

                    a_wn=fmin(1.0,2.0*r_wn);
                    b_wn=fmin(2.0,r_wn);
                    c_wn=fmax(a_wn,b_wn);
                    limiter_wn=fmax(0.0,c_wn);

					una=u[i][j+1][k]+0.5*limiter_un*(u[i][j][k]-u[i][j+1][k]);
					vna=v[i][j+1][k]+0.5*limiter_vn*(v[i][j][k]-v[i][j+1][k]);
                    wna=w[i][j+1][k]+0.5*limiter_wn*(w[i][j][k]-w[i][j+1][k]);
				}
		if(vs>0.0)
			{

                r_us=(u[i][j-1][k]-u[i][j-2][k])/(u[i][j][k]-u[i][j-1][k]);
                r_vs=(v[i][j-1][k]-v[i][j-2][k])/(v[i][j][k]-v[i][j-1][k]);
                r_ws=(w[i][j-1][k]-w[i][j-2][k])/(w[i][j][k]-w[i][j-1][k]);

                a_us=fmin(1.0,2.0*r_us);
                b_us=fmin(2.0,r_us);
                c_us=fmax(a_un,b_us);
                limiter_us=fmax(0.0,c_us);

                a_vs=fmin(1.0,2.0*r_vs);
                b_vs=fmin(2.0,r_vs);
                c_vs=fmax(a_vs,b_vs);
                limiter_vs=fmax(0.0,c_vs);

                a_ws=fmin(1.0,2.0*r_ws);
                b_ws=fmin(2.0,r_ws);
                c_ws=fmax(a_ws,b_ws);
                limiter_ws=fmax(0.0,c_ws);

				usa=u[i][j-1][k]+0.5*limiter_us*(u[i][j][k]-u[i][j-1][k]);
				vsa=v[i][j-1][k]+0.5*limiter_vs*(v[i][j][k]-v[i][j-1][k]);
				wsa=w[i][j-1][k]+0.5*limiter_ws*(w[i][j][k]-w[i][j-1][k]);
			}else
				{
                    r_us=(u[i][j][k]-u[i][j-1][k])/(u[i][j-1][k]-u[i][j-2][k]);
                    r_vs=(v[i][j][k]-v[i][j-1][k])/(v[i][j-1][k]-v[i][j-2][k]);
                    r_ws=(w[i][j][k]-w[i][j-1][k])/(w[i][j-1][k]-w[i][j-2][k]);

                    a_us=fmin(1.0,2.0*r_us);
                    b_us=fmin(2.0,r_us);
                    c_us=fmax(a_un,b_us);
                    limiter_us=fmax(0.0,c_us);

                    a_vs=fmin(1.0,2.0*r_vs);
                    b_vs=fmin(2.0,r_vs);
                    c_vs=fmax(a_vs,b_vs);
                    limiter_vs=fmax(0.0,c_vs);

                    a_ws=fmin(1.0,2.0*r_ws);
                    b_ws=fmin(2.0,r_ws);
                    c_ws=fmax(a_ws,b_ws);
                    limiter_ws=fmax(0.0,c_ws);

					usa=u[i][j][k]+0.5*limiter_us*(u[i][j-1][k]-u[i][j][k]);
					vsa=v[i][j][k]+0.5*limiter_vs*(v[i][j-1][k]-v[i][j][k]);
                    wsa=w[i][j][k]+0.5*limiter_ws*(w[i][j-1][k]-w[i][j][k]);
				}
		if(wt>0.0)
			{
                r_ut=(u[i][j][k]-u[i][j][k-1])/(u[i][j][k+1]-u[i][j][k]);
                r_vt=(v[i][j][k]-v[i][j][k-1])/(v[i][j][k+1]-v[i][j][k]);
                r_wt=(w[i][j][k]-w[i][j][k-1])/(w[i][j][k+1]-w[i][j][k]);

                a_ut=fmin(1.0,2.0*r_ut);
                b_ut=fmin(2.0,r_ut);
                c_ut=fmax(a_ut,b_ut);
                limiter_ut=fmax(0.0,c_ut);

                a_vt=fmin(1.0,2.0*r_vt);
                b_vt=fmin(2.0,r_vt);
                c_vt=fmax(a_vt,b_vt);
                limiter_vt=fmax(0.0,c_vt);

                a_wt=fmin(1.0,2.0*r_wt);
                b_wt=fmin(2.0,r_wt);
                c_wt=fmax(a_wt,b_wt);
                limiter_wt=fmax(0.0,c_wt);


				uta=u[i][j][k]+0.5*limiter_ut*(u[i][j][k+1]-u[i][j][k]);
				vta=v[i][j][k]+0.5*limiter_vt*(v[i][j][k+1]-v[i][j][k]);
                wta=w[i][j][k]+0.5*limiter_wt*(w[i][j][k+1]-w[i][j][k]);

			}else
				{

                    r_ut=(u[i][j][k+1]-u[i][j][k])/(u[i][j][k]-u[i][j][k-1]);
                    r_vt=(v[i][j][k+1]-v[i][j][k])/(v[i][j][k]-v[i][j][k-1]);
                    r_wt=(w[i][j][k+1]-w[i][j][k])/(w[i][j][k]-w[i][j][k-1]);


                    a_ut=fmin(1.0,2.0*r_ut);
                    b_ut=fmin(2.0,r_ut);
                    c_ut=fmax(a_ut,b_ut);
                    limiter_ut=fmax(0.0,c_ut);

                    a_vt=fmin(1.0,2.0*r_vt);
                    b_vt=fmin(2.0,r_vt);
                    c_vt=fmax(a_vt,b_vt);
                    limiter_vt=fmax(0.0,c_vt);

                    a_wt=fmin(1.0,2.0*r_wt);
                    b_wt=fmin(2.0,r_wt);
                    c_wt=fmax(a_wt,b_wt);
                    limiter_wt=fmax(0.0,c_wt);

					uta=u[i][j][k+1]+0.5*limiter_ut*(u[i][j][k]-u[i][j][k+1]);
					vta=v[i][j][k+1]+0.5*limiter_vt*(v[i][j][k]-v[i][j][k+1]);
                    wta=w[i][j][k+1]+0.5*limiter_wt*(w[i][j][k]-w[i][j][k+1]);
				}
		if(wb>0.0)
			{

                r_ub=(u[i][j][k-1]-u[i][j][k-2])/(u[i][j][k]-u[i][j][k-1]);
                r_vb=(v[i][j][k-1]-v[i][j][k-2])/(v[i][j][k]-v[i][j][k-1]);
                r_wb=(w[i][j][k-1]-w[i][j][k-2])/(w[i][j][k]-w[i][j][k-1]);


                a_ub=fmin(1.0,2.0*r_ub);
                b_ub=fmin(2.0,r_ub);
                c_ub=fmax(a_ub,b_ub);
                limiter_ub=fmax(0.0,c_ub);

                a_vb=fmin(1.0,2.0*r_vb);
                b_vb=fmin(2.0,r_vb);
                c_vb=fmax(a_vb,b_vb);
                limiter_vb=fmax(0.0,c_vb);

                a_wb=fmin(1.0,2.0*r_wb);
                b_wb=fmin(2.0,r_wb);
                c_wb=fmax(a_wb,b_wb);
                limiter_wb=fmax(0.0,c_wb);

				uba=u[i][j][k-1]+0.5*limiter_ub*(u[i][j][k]-u[i][j][k-1]);
				vba=v[i][j][k-1]+0.5*limiter_vb*(v[i][j][k]-v[i][j][k-1]);
				wba=w[i][j][k-1]+0.5*limiter_wb*(w[i][j][k]-w[i][j][k-1]);
			}else
				{
                    r_ub=(u[i][j][k]-u[i][j][k-1])/(u[i][j][k-1]-u[i][j][k-2]);
                    r_vb=(v[i][j][k]-v[i][j][k-1])/(v[i][j][k-1]-v[i][j][k-2]);
                    r_wb=(w[i][j][k]-w[i][j][k-1])/(w[i][j][k-1]-w[i][j][k-2]);

                    a_ub=fmin(1.0,2.0*r_ub);
                    b_ub=fmin(2.0,r_ub);
                    c_ub=fmax(a_ub,b_ub);
                    limiter_ub=fmax(0.0,c_ub);

                    a_vb=fmin(1.0,2.0*r_vb);
                    b_vb=fmin(2.0,r_vb);
                    c_vb=fmax(a_vb,b_vb);
                    limiter_vb=fmax(0.0,c_vb);

                    a_wb=fmin(1.0,2.0*r_wb);
                    b_wb=fmin(2.0,r_wb);
                    c_wb=fmax(a_wb,b_wb);
                    limiter_wb=fmax(0.0,c_wb);

					uba=u[i][j][k]+0.5*limiter_ub*(u[i][j][k-1]-u[i][j][k]);
					vba=v[i][j][k]+0.5*limiter_vb*(v[i][j][k-1]-v[i][j][k]);
					wba=w[i][j][k]+0.5*limiter_wb*(w[i][j][k-1]-w[i][j][k]);
				}
}

void Superbee1()
{

       //TVD face velocities
		if(ue1>0.0)
			{
                r_ue1=(u1[i][j][k]-u1[i-1][j][k])/(u1[i+1][j][k]-u1[i][j][k]);
                r_ve1=(v1[i][j][k]-v1[i-1][j][k])/(v1[i+1][j][k]-v1[i][j][k]);
                r_we1=(w1[i][j][k]-w1[i-1][j][k])/(w1[i+1][j][k]-w1[i][j][k]);

                a_ue=fmin(1.0,2.0*r_ue);
                b_ue=fmin(2.0,r_ue);
                c_ue=fmax(a_ue,b_ue);
                limiter_ue=fmax(0.0,c_ue);

                a_ve=fmin(1.0,2.0*r_ve);
                b_ve=fmin(2.0,r_ve);
                c_ve=fmax(a_ve,b_ve);
                limiter_ve=fmax(0.0,c_ve);

                a_we=fmin(1.0,2.0*r_we);
                b_we=fmin(2.0,r_we);
                c_we=fmax(a_we,b_we);
                limiter_we=fmax(0.0,c_we);

				uea1=u1[i][j][k]+0.5*limiter_ue*(u1[i+1][j][k]-u1[i][j][k]);
				vea1=v1[i][j][k]+0.5*limiter_ve*(v1[i+1][j][k]-v1[i][j][k]);
				wea1=w1[i][j][k]+0.5*limiter_we*(w1[i+1][j][k]-w1[i][j][k]);
			}else
				{
                    r_ue=(u1[i+1][j][k]-u1[i][j][k])/(u1[i][j][k]-u1[i-1][j][k]);
                    r_ve=(v1[i+1][j][k]-v1[i][j][k])/(v1[i][j][k]-v1[i-1][j][k]);
                    r_we=(w1[i+1][j][k]-w1[i][j][k])/(w1[i][j][k]-w1[i-1][j][k]);

                    a_ue=fmin(1.0,2.0*r_ue);
                    b_ue=fmin(2.0,r_ue);
                    c_ue=fmax(a_ue,b_ue);
                    limiter_ue=fmax(0.0,c_ue);

                    a_ve=fmin(1.0,2.0*r_ve);
                    b_ve=fmin(2.0,r_ve);
                    c_ve=fmax(a_ve,b_ve);
                    limiter_ve=fmax(0.0,c_ve);

                    a_we=fmin(1.0,2.0*r_we);
                    b_we=fmin(2.0,r_we);
                    c_we=fmax(a_we,b_we);
                    limiter_we=fmax(0.0,c_we);

					uea1=u1[i+1][j][k]+0.5*limiter_ue*(u1[i][j][k]-u1[i+1][j][k]);
					vea1=v1[i+1][j][k]+0.5*limiter_ve*(v1[i][j][k]-v1[i+1][j][k]);
					wea1=w1[i+1][j][k]+0.5*limiter_we*(w1[i][j][k]-w1[i+1][j][k]);
				}
		if(uw1>0.0)
			{

                r_uw=(u1[i-1][j][k]-u1[i-2][j][k])/(u1[i][j][k]-u1[i-1][j][k]);
                r_vw=(v1[i-1][j][k]-v1[i-2][j][k])/(v1[i][j][k]-v1[i-1][j][k]);
                r_ww=(w1[i-1][j][k]-w1[i-2][j][k])/(w1[i][j][k]-w1[i-1][j][k]);

                a_uw=fmin(1.0,2.0*r_uw);
                b_uw=fmin(2.0,r_uw);
                c_uw=fmax(a_uw,b_uw);
                limiter_uw=fmax(0.0,c_uw);

                a_vw=fmin(1.0,2.0*r_vw);
                b_vw=fmin(2.0,r_vw);
                c_vw=fmax(a_vw,b_vw);
                limiter_vw=fmax(0.0,c_vw);

                a_ww=fmin(1.0,2.0*r_ww);
                b_ww=fmin(2.0,r_ww);
                c_ww=fmax(a_ww,b_ww);
                limiter_ww=fmax(0.0,c_ww);

				uwa1=u1[i-1][j][k]+0.5*limiter_uw*(u1[i][j][k]-u1[i-1][j][k]);
				vwa1=v1[i-1][j][k]+0.5*limiter_vw*(v1[i][j][k]-v1[i-1][j][k]);
				wwa1=w1[i-1][j][k]+0.5*limiter_ww*(w1[i][j][k]-w1[i-1][j][k]);

			}else
				{
                    r_uw=(u1[i][j][k]-u1[i-1][j][k])/(u1[i-1][j][k]-u1[i-2][j][k]);
                    r_vw=(v1[i][j][k]-v1[i-1][j][k])/(v1[i-1][j][k]-v1[i-2][j][k]);
                    r_ww=(w1[i][j][k]-w1[i-1][j][k])/(w1[i-1][j][k]-w1[i-2][j][k]);

                    a_uw=fmin(1.0,2.0*r_uw);
                    b_uw=fmin(2.0,r_uw);
                    c_uw=fmax(a_uw,b_uw);
                    limiter_uw=fmax(0.0,c_uw);

                    a_vw=fmin(1.0,2.0*r_vw);
                    b_vw=fmin(2.0,r_vw);
                    c_vw=fmax(a_vw,b_vw);
                    limiter_vw=fmax(0.0,c_vw);

                    a_ww=fmin(1.0,2.0*r_ww);
                    b_ww=fmin(2.0,r_ww);
                    c_ww=fmax(a_ww,b_ww);
                    limiter_ww=fmax(0.0,c_ww);

					uwa1=u1[i][j][k]+0.5*limiter_uw*(u1[i-1][j][k]-u1[i][j][k]);
					vwa1=v1[i][j][k]+0.5*limiter_vw*(v1[i-1][j][k]-v1[i][j][k]);
                    wwa1=w1[i][j][k]+0.5*limiter_ww*(w1[i-1][j][k]-w1[i][j][k]);
				}
		if(vn1>0.0)
			{
                r_un=(u1[i][j][k]-u1[i][j-1][k])/(u1[i][j+1][k]-u1[i][j][k]);
                r_vn=(v1[i][j][k]-v1[i][j-1][k])/(v1[i][j+1][k]-v1[i][j][k]);
                r_wn=(w1[i][j][k]-w1[i][j-1][k])/(w1[i][j+1][k]-w1[i][j][k]);

                a_un=fmin(1.0,2.0*r_un);
                b_un=fmin(2.0,r_un);
                c_un=fmax(a_un,b_un);
                limiter_un=fmax(0.0,c_un);

                a_vn=fmin(1.0,2.0*r_vn);
                b_vn=fmin(2.0,r_vn);
                c_vn=fmax(a_vn,b_vn);
                limiter_vn=fmax(0.0,c_vn);

                a_wn=fmin(1.0,2.0*r_wn);
                b_wn=fmin(2.0,r_wn);
                c_wn=fmax(a_wn,b_wn);
                limiter_wn=fmax(0.0,c_wn);


				una1=u1[i][j][k]+0.5*limiter_un*(u1[i][j+1][k]-u1[i][j][k]);
				vna1=v1[i][j][k]+0.5*limiter_vn*(v1[i][j+1][k]-v1[i][j][k]);
                wna1=w1[i][j][k]+0.5*limiter_wn*(w1[i][j+1][k]-w1[i][j][k]);
			}else
				{
                    r_un=(u1[i][j+1][k]-u1[i][j][k])/(u1[i][j][k]-u1[i][j-1][k]);
                    r_vn=(v1[i][j+1][k]-v1[i][j][k])/(v1[i][j][k]-v1[i][j-1][k]);
                    r_wn=(w1[i][j+1][k]-w1[i][j][k])/(w1[i][j][k]-w1[i][j-1][k]);

                    a_un=fmin(1.0,2.0*r_un);
                    b_un=fmin(2.0,r_un);
                    c_un=fmax(a_un,b_un);
                    limiter_un=fmax(0.0,c_un);

                    a_vn=fmin(1.0,2.0*r_vn);
                    b_vn=fmin(2.0,r_vn);
                    c_vn=fmax(a_vn,b_vn);
                    limiter_vn=fmax(0.0,c_vn);

                    a_wn=fmin(1.0,2.0*r_wn);
                    b_wn=fmin(2.0,r_wn);
                    c_wn=fmax(a_wn,b_wn);
                    limiter_wn=fmax(0.0,c_wn);

					una1=u1[i][j+1][k]+0.5*limiter_un*(u1[i][j][k]-u1[i][j+1][k]);
					vna1=v1[i][j+1][k]+0.5*limiter_vn*(v1[i][j][k]-v1[i][j+1][k]);
                    wna1=w1[i][j+1][k]+0.5*limiter_wn*(w1[i][j][k]-w1[i][j+1][k]);
				}
		if(vs1>0.0)
			{

                r_us=(u1[i][j-1][k]-u1[i][j-2][k])/(u1[i][j][k]-u1[i][j-1][k]);
                r_vs=(v1[i][j-1][k]-v1[i][j-2][k])/(v1[i][j][k]-v1[i][j-1][k]);
                r_ws=(w1[i][j-1][k]-w1[i][j-2][k])/(w1[i][j][k]-w1[i][j-1][k]);

                a_us=fmin(1.0,2.0*r_us);
                b_us=fmin(2.0,r_us);
                c_us=fmax(a_un,b_us);
                limiter_us=fmax(0.0,c_us);

                a_vs=fmin(1.0,2.0*r_vs);
                b_vs=fmin(2.0,r_vs);
                c_vs=fmax(a_vs,b_vs);
                limiter_vs=fmax(0.0,c_vs);

                a_ws=fmin(1.0,2.0*r_ws);
                b_ws=fmin(2.0,r_ws);
                c_ws=fmax(a_ws,b_ws);
                limiter_ws=fmax(0.0,c_ws);

				usa1=u1[i][j-1][k]+0.5*limiter_us*(u1[i][j][k]-u1[i][j-1][k]);
				vsa1=v1[i][j-1][k]+0.5*limiter_vs*(v1[i][j][k]-v1[i][j-1][k]);
				wsa1=w1[i][j-1][k]+0.5*limiter_ws*(w1[i][j][k]-w1[i][j-1][k]);
			}else
				{
                    r_us=(u1[i][j][k]-u1[i][j-1][k])/(u1[i][j-1][k]-u1[i][j-2][k]);
                    r_vs=(v1[i][j][k]-v1[i][j-1][k])/(v1[i][j-1][k]-v1[i][j-2][k]);
                    r_ws=(w1[i][j][k]-w1[i][j-1][k])/(w1[i][j-1][k]-w1[i][j-2][k]);

                    a_us=fmin(1.0,2.0*r_us);
                    b_us=fmin(2.0,r_us);
                    c_us=fmax(a_un,b_us);
                    limiter_us=fmax(0.0,c_us);

                    a_vs=fmin(1.0,2.0*r_vs);
                    b_vs=fmin(2.0,r_vs);
                    c_vs=fmax(a_vs,b_vs);
                    limiter_vs=fmax(0.0,c_vs);

                    a_ws=fmin(1.0,2.0*r_ws);
                    b_ws=fmin(2.0,r_ws);
                    c_ws=fmax(a_ws,b_ws);
                    limiter_ws=fmax(0.0,c_ws);

					usa1=u1[i][j][k]+0.5*limiter_us*(u1[i][j-1][k]-u1[i][j][k]);
					vsa1=v1[i][j][k]+0.5*limiter_vs*(v1[i][j-1][k]-v1[i][j][k]);
                    wsa1=w1[i][j][k]+0.5*limiter_ws*(w1[i][j-1][k]-w1[i][j][k]);
				}
		if(wt1>0.0)
			{
                r_ut=(u1[i][j][k]-u1[i][j][k-1])/(u1[i][j][k+1]-u1[i][j][k]);
                r_vt=(v1[i][j][k]-v1[i][j][k-1])/(v1[i][j][k+1]-v1[i][j][k]);
                r_wt=(w1[i][j][k]-w1[i][j][k-1])/(w1[i][j][k+1]-w1[i][j][k]);

                a_ut=fmin(1.0,2.0*r_ut);
                b_ut=fmin(2.0,r_ut);
                c_ut=fmax(a_ut,b_ut);
                limiter_ut=fmax(0.0,c_ut);

                a_vt=fmin(1.0,2.0*r_vt);
                b_vt=fmin(2.0,r_vt);
                c_vt=fmax(a_vt,b_vt);
                limiter_vt=fmax(0.0,c_vt);

                a_wt=fmin(1.0,2.0*r_wt);
                b_wt=fmin(2.0,r_wt);
                c_wt=fmax(a_wt,b_wt);
                limiter_wt=fmax(0.0,c_wt);


				uta1=u1[i][j][k]+0.5*limiter_ut*(u1[i][j][k+1]-u1[i][j][k]);
				vta1=v1[i][j][k]+0.5*limiter_vt*(v1[i][j][k+1]-v1[i][j][k]);
                wta1=w1[i][j][k]+0.5*limiter_wt*(w1[i][j][k+1]-w1[i][j][k]);

			}else
				{

                    r_ut=(u1[i][j][k+1]-u1[i][j][k])/(u1[i][j][k]-u1[i][j][k-1]);
                    r_vt=(v1[i][j][k+1]-v1[i][j][k])/(v1[i][j][k]-v1[i][j][k-1]);
                    r_wt=(w1[i][j][k+1]-w1[i][j][k])/(w1[i][j][k]-w1[i][j][k-1]);


                    a_ut=fmin(1.0,2.0*r_ut);
                    b_ut=fmin(2.0,r_ut);
                    c_ut=fmax(a_ut,b_ut);
                    limiter_ut=fmax(0.0,c_ut);

                    a_vt=fmin(1.0,2.0*r_vt);
                    b_vt=fmin(2.0,r_vt);
                    c_vt=fmax(a_vt,b_vt);
                    limiter_vt=fmax(0.0,c_vt);

                    a_wt=fmin(1.0,2.0*r_wt);
                    b_wt=fmin(2.0,r_wt);
                    c_wt=fmax(a_wt,b_wt);
                    limiter_wt=fmax(0.0,c_wt);

					uta1=u1[i][j][k+1]+0.5*limiter_ut*(u1[i][j][k]-u1[i][j][k+1]);
					vta1=v1[i][j][k+1]+0.5*limiter_vt*(v1[i][j][k]-v1[i][j][k+1]);
                    wta1=w1[i][j][k+1]+0.5*limiter_wt*(w1[i][j][k]-w1[i][j][k+1]);
				}
		if(wb1>0.0)
			{

                r_ub=(u1[i][j][k-1]-u1[i][j][k-2])/(u1[i][j][k]-u1[i][j][k-1]);
                r_vb=(v1[i][j][k-1]-v1[i][j][k-2])/(v1[i][j][k]-v1[i][j][k-1]);
                r_wb=(w1[i][j][k-1]-w1[i][j][k-2])/(w1[i][j][k]-w1[i][j][k-1]);


                a_ub=fmin(1.0,2.0*r_ub);
                b_ub=fmin(2.0,r_ub);
                c_ub=fmax(a_ub,b_ub);
                limiter_ub=fmax(0.0,c_ub);

                a_vb=fmin(1.0,2.0*r_vb);
                b_vb=fmin(2.0,r_vb);
                c_vb=fmax(a_vb,b_vb);
                limiter_vb=fmax(0.0,c_vb);

                a_wb=fmin(1.0,2.0*r_wb);
                b_wb=fmin(2.0,r_wb);
                c_wb=fmax(a_wb,b_wb);
                limiter_wb=fmax(0.0,c_wb);

				uba1=u1[i][j][k-1]+0.5*limiter_ub*(u1[i][j][k]-u1[i][j][k-1]);
				vba1=v1[i][j][k-1]+0.5*limiter_vb*(v1[i][j][k]-v1[i][j][k-1]);
				wba1=w1[i][j][k-1]+0.5*limiter_wb*(w1[i][j][k]-w1[i][j][k-1]);
			}else
				{
                    r_ub=(u1[i][j][k]-u1[i][j][k-1])/(u1[i][j][k-1]-u1[i][j][k-2]);
                    r_vb=(v1[i][j][k]-v1[i][j][k-1])/(v1[i][j][k-1]-v1[i][j][k-2]);
                    r_wb=(w1[i][j][k]-w1[i][j][k-1])/(w1[i][j][k-1]-w1[i][j][k-2]);

                    a_ub=fmin(1.0,2.0*r_ub);
                    b_ub=fmin(2.0,r_ub);
                    c_ub=fmax(a_ub,b_ub);
                    limiter_ub=fmax(0.0,c_ub);

                    a_vb=fmin(1.0,2.0*r_vb);
                    b_vb=fmin(2.0,r_vb);
                    c_vb=fmax(a_vb,b_vb);
                    limiter_vb=fmax(0.0,c_vb);

                    a_wb=fmin(1.0,2.0*r_wb);
                    b_wb=fmin(2.0,r_wb);
                    c_wb=fmax(a_wb,b_wb);
                    limiter_wb=fmax(0.0,c_wb);

					uba1=u1[i][j][k]+0.5*limiter_ub*(u1[i][j][k-1]-u1[i][j][k]);
					vba1=v1[i][j][k]+0.5*limiter_vb*(v1[i][j][k-1]-v1[i][j][k]);
					wba1=w1[i][j][k]+0.5*limiter_wb*(w1[i][j][k-1]-w1[i][j][k]);
				}
}
