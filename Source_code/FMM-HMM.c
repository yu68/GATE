/* ****** VERSION HISTORY ******

06/11/2011

****** END VERSION HISTORY ****** */

#define NRANSI
#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <time.h>
#include "nrutil.h"
#include "FMM-HMM.h"

#define LARGENUM 1e50
#define SMALLNUM 1e-50


// ********** random number generator **********

// The relation between a 3d r matrix to the relative vector
// Suppose x is n*m*v dimension matrix, y = as.vector(x), then x[i,j,k]=y[i+(j-1)*n+(k-1)*(n*m)]
// Caution v is not included!!!
double emission(double *observation, double *lambda, long nstate, long nmarker, long ntime, long nregion, long ncluster, 
long k, long w, long h, long t)
{
 double p = 1.0;
 long m;
 for(m=1;m<=nmarker;m++)
 {
  p *= dpois(observation[w+(t-1)*nregion+(m-1)*(nregion*ntime)],lambda[k+h*ncluster+(m-1)*(ncluster*nstate)],0);
 }
 
 return p;
}


void forward(double *f, double *b, double *lambda, double *observation, long nstate, long nmarker, long ntime, long nregion, long ncluster)
{
 long w, k, t;
 for(w = 1; w<= nregion*ncluster*2*ntime; w++)
    {
	 f[w] = 1;
	}
 for(w=1; w<=nregion; w++)
    {
	  for(k=1; k<=ncluster; k++)
	     {
		    f[w+(k-1)*nregion] = emission(observation, lambda, nstate, nmarker, ntime, nregion, ncluster, k, w, 0, 1);
			f[w+(k-1)*nregion+nregion*ncluster] = 0;
			for(t=2; t<=ntime; t++)
			   {
			    f[w+(k-1)*nregion+(t-1)*nregion*ncluster*2] = emission(observation, lambda, nstate, nmarker, ntime, nregion, ncluster, k, w, 0, t)*
				(f[w+(k-1)*nregion+(t-2)*nregion*ncluster*2]*b[k]+f[w+(k-1)*nregion+nregion*ncluster+(t-2)*nregion*ncluster*2]*b[k+ncluster]);
				f[w+(k-1)*nregion+nregion*ncluster+(t-1)*nregion*ncluster*2] = emission(observation, lambda, nstate, nmarker, ntime, nregion, ncluster, k, w, 1, t)
				*(f[w+(k-1)*nregion+(t-2)*nregion*ncluster*2]*b[k+ncluster*2]+f[w+(k-1)*nregion+nregion*ncluster+(t-2)*nregion*ncluster*2]*b[k+ncluster*3]);
			   }
		 }
	} 
}


void backward(double *back, double *b, double *lambda, double *observation, long nstate, long nmarker, long ntime, long nregion, long ncluster)
{
 long w, k, t;
 for(w = 1; w<= nregion*ncluster*2*ntime; w++)
    {
	 back[w] = 1;
	}
 for(w = 1; w <= nregion; w++)
    {
	 for(k = 1; k <= ncluster; k++)
	    {
		 back[w+(k-1)*nregion+nregion*ncluster*2*(ntime-1)] = 1;
		 back[w+(k-1)*nregion+nregion*ncluster+nregion*ncluster*2*(ntime-1)] = 1;
		 for(t = ntime-1; t >= 1; t--)
		    {
			 back[w+(k-1)*nregion+nregion*ncluster*2*(t-1)] = back[w+(k-1)*nregion+nregion*ncluster*2*t]*b[k]*emission(observation, lambda, nstate, nmarker, ntime, nregion, ncluster, k, w, 0, t+1)+
			 back[w+(k-1)*nregion+nregion*ncluster+nregion*ncluster*2*t]*b[k+ncluster*2]*emission(observation, lambda, nstate, nmarker, ntime, nregion, ncluster, k, w, 1, t+1);
			 back[w+(k-1)*nregion+nregion*ncluster+nregion*ncluster*2*(t-1)] = back[w+(k-1)*nregion+nregion*ncluster*2*t]*b[k+ncluster]*emission(observation, lambda, nstate, nmarker, ntime, nregion, ncluster, k, w, 0, t+1)+
			 back[w+(k-1)*nregion+nregion*ncluster+nregion*ncluster*2*t]*b[k+ncluster+ncluster*2]*emission(observation, lambda, nstate, nmarker, ntime, nregion, ncluster, k, w, 1, t+1);
			}
		}
	}
 

}


double log_likelihood(long nregion, long ncluster, long ntime, long *cluster, double *pie, double *forwar)
{
 double P = 0.0;
 double Q = 0.0;
 long w, k;
 for(w = 1; w<=nregion; w++)
 {
	 P=0.0;
	 for (k=1; k <= ncluster; k++)
	 {
	     P += (forwar[w+(k-1)*nregion+2*nregion*ncluster*(ntime-1)]+forwar[w+(k-1)*nregion+nregion*ncluster+2*nregion*ncluster*(ntime-1)])*pie[k];
	
	 }
	 Q +=log(P);
 }
 return Q;
}


//Maximum function
double f_max(double *a, long a_length, long *index)
{
 long i, max = 1;
 double x = a[1];
 for(i=2;i<=a_length;i++)
 {
    if(x<a[i]) {
	  x = a[i];
	  max = i;
	}
 }
 if(index != 0) {
	(*index) = max;
 }
 return x;
}

// R which function
long f_which(double *a, long a_length, double b)
{
 double r = 1.0e-100;
 long i,n;
 for(i=1; i<=a_length; i++)
 {
    if(fabs(a[i]-b)<r)
	{
	 n = i;
	 break;
	}
	
 }
 return n;
}


double *rowSums(double **probability, long nregion, long ncluster)
{
 long i, j;
 double *a = dvector(1, nregion);
 
 for(i=1; i<=nregion; i++)
    {
	 a[i] = 0;
	 for(j=1; j<=ncluster; j++)
	    {
		 a[i]+=probability[i][j];
		}
	}
 return a;
}

void matrixtovector(double **matrixf, double *vectorf, long nregion, long ncluster)
{
 long i, j;
 for(i=1; i<=nregion; i++)
    {
	 for(j=1; j<=ncluster; j++)
	    {
		 matrixf[i][j] = matrixf[i][j]/vectorf[i];
		}
	}
}

void colMeans(double **probability, double *pie, long nregion, long ncluster)
{
 long i, j;
 for(j=1; j<=ncluster; j++)
    {
	 pie[j] = 0;
	 for(i=1; i<=nregion; i++)
	    {
		 pie[j] += probability[i][j];
		}
	 pie[j] = pie[j]/nregion;
	}
}

// probability is nregion*ncluster matrix
// probability should be assigned memory space beforehand.
void cluster_matrix(long nregion, long ncluster, long ntime, long *cluster, double *pie, double *forwar, double **probability)
{
 long w, ki;
 long i, j;
 double *prob_sum;
 double a;
 for(i=1; i<=nregion; i++)
     {
	    for(j=1; j<=ncluster; j++)
		   probability[i][j] = 0;
	 }

 
 for(w=1; w<=nregion; w++)
    {
	 for(ki=1; ki<=ncluster; ki++)
	    {
		 probability[w][ki] = LARGENUM*(forwar[w+(ki-1)*nregion+2*nregion*ncluster*(ntime-1)]+forwar[w+(ki-1)*nregion+nregion*ncluster+2*nregion*ncluster*(ntime-1)])*pie[ki];
		}
	 a = f_max(probability[w], ncluster, cluster+w);
     
     	 
	} 
 
 prob_sum = rowSums(probability, nregion, ncluster);
 matrixtovector(probability, prob_sum, nregion, ncluster);
 colMeans(probability, pie, nregion, ncluster);
 free_dvector(prob_sum, 1, nregion);
 
}



// Update lambda and b.
// Add Hidden state H: nregion * ntime;
void HMM_Recursion(double *b, double *observation, double *lambda, long *cluster, double *forwar, double *backwar, long nstate, 
long ncluster, long nmarker, long ntime, long nregion, long *H)
{
 double *new_b = dvector(1, ncluster*2*2);
 double *new_lambda = dvector(1, ncluster*2*nmarker);
 double *denomi_lambda = dvector(1, ncluster*2*nmarker);
 long i, j, w, k, t, m;
 long b_l = ncluster*2*2;
 long lam_l = ncluster*2*nmarker;
 double p_obser_w, temp, b_sum;
 for(i=1; i<=b_l; i++)
   {
    new_b[i] = 0;
   }
 for(i=1; i<=lam_l; i++)
   {
    new_lambda[i] = 0;
	denomi_lambda[i] = 0;
   }
 for(w=1; w<=nregion; w++)
    {
	 k = cluster[w];
	 p_obser_w = forwar[w+(k-1)*nregion+2*nregion*ncluster*2] + forwar[w+(k-1)*nregion+nregion*ncluster+2*nregion*ncluster*2];
	 for(i=1; i<=2; i++)
	    {
		 for(j=1; j<=2; j++)
		    {
			 for(t=1; t<=ntime-1; t++)
			    {
				 temp = LARGENUM*forwar[w+(k-1)*nregion+(i-1)*nregion*ncluster+(t-1)*nregion*ncluster*2]*b[k+(i-1)*ncluster+(j-1)*ncluster*2]*backwar[w+(k-1)*nregion+(j-1)*nregion*ncluster+t*nregion*ncluster*2]*
				        emission(observation, lambda, nstate, nmarker, ntime, nregion, ncluster,k,w,j-1,t+1)/p_obser_w;
				 new_b[k+(i-1)*ncluster+(j-1)*2*ncluster] += temp;
				}
			}
		 for(t=1; t<=ntime; t++)
		    {
			 H[w+(t-1)*nregion] = ((forwar[w+(k-1)*nregion+(t-1)*nregion*ncluster*2]*backwar[w+(k-1)*nregion+(t-1)*nregion*ncluster*2])<
			                       (forwar[w+(k-1)*nregion+nregion*ncluster+(t-1)*nregion*ncluster*2]*backwar[w+(k-1)*nregion+nregion*ncluster+(t-1)*nregion*ncluster*2]) ? 1:0);
			 for(m=1; m<=nmarker; m++)
			    {
				 new_lambda[k+(i-1)*ncluster+(m-1)*ncluster*2] += abs(i-2+H[w+(t-1)*nregion])*observation[w+(t-1)*nregion+(m-1)*nregion*ntime];
				 denomi_lambda[k+(i-1)*ncluster+(m-1)*ncluster*2] += abs(i-2+H[w+(t-1)*nregion]);
				}
			}
		}
	}
	
 for(k=1; k<=ncluster; k++)
    {
	 for(i=1; i<=2; i++)
	    {
		 b_sum = new_b[k+(i-1)*ncluster] + new_b[k+(i-1)*ncluster+ncluster*2];
		 b[k+(i-1)*ncluster] = new_b[k+(i-1)*ncluster]/b_sum;
		 b[k+(i-1)*ncluster+ncluster*2] = new_b[k+(i-1)*ncluster+ncluster*2]/b_sum;
		}
	}
 for(k=1; k<=ncluster; k++)
    {
	 for(i=1; i<=2; i++)
	    {
		 for(m=1; m<=nmarker; m++)
		    {
			 lambda[k+(i-1)*ncluster+(m-1)*ncluster*2] = new_lambda[k+(i-1)*ncluster+(m-1)*ncluster*2]/denomi_lambda[k+(i-1)*ncluster+(m-1)*ncluster*2];
			}
		}
	}
 
 
 free_dvector(new_b, 1, ncluster*2*2);
 free_dvector(new_lambda, 1, ncluster*2*nmarker);
 free_dvector(denomi_lambda, 1, ncluster*2*nmarker);
 
}

double Edistance(double *a, double *b, long length)
{
 double c=0;
 long i;
 for(i=1; i<=length; i++)
 {
    c += pow(a[i] - b[i],2);
 }
 
 return sqrt(c);
}

// Compute distance of two matrix;
double mEdistance(double **a, double **b, long nrow, long ncol)
{
 long i,j;
 double d=0;
 for(i=1; i<=nrow; i++)
 {
    for(j=1; j<=ncol; j++)
	   {
	    d += pow(a[i][j]-b[i][j],2);
	   }
 }
 return sqrt(d);
}

//Copy b to a.
void copy_vector(double *a, double *b, long length)
{
 long i;
 for(i=1;i<=length;i++)
    {
	 a[i] = b[i];
	}
}

// Copy matrix b to a.
void copy_matrix(double **a, double **b, long nrow, long ncol)
{
 long i, j;
 for(i=1;i<=nrow;i++)
 {
    for(j=1; j<=ncol; j++)
	{
	   a[i][j] = b[i][j];
	}
 }

}




// Main function
void f_HMM(double *observation, double *b, double *lambda, int32_t *cluster, double *pie, double *forwar, double *backwar, 
double *probability, double *diff_r, double *log_likeli_r, double *pie_r, int32_t *nstate, int32_t *ncluster, int32_t *nmarker, int32_t *ntime, int32_t *nregion, 
int32_t *maxiteration, int32_t *stop, int32_t *nstep, double *ndistance, int32_t *H)
{
 long nstate_c = *nstate, ncluster_c = *ncluster, nmarker_c = *nmarker, ntime_c = *ntime, nregion_c = *nregion, maxiteration_c = *maxiteration,
      nstep_c = *nstep;
 double ndistance_c = *ndistance;
 long i, j, k, l;
 double *observation_c = dvector(1, nregion_c*ntime_c*nmarker_c);
 double *b_c = dvector(1, ncluster_c*4);
 double *lambda_c = dvector(1, ncluster_c*2*nmarker_c);
 long *cluster_c = lvector(1, nregion_c);
 double *pie_c = dvector(1, ncluster_c);
 double *forwar_c = dvector(1, nregion_c*ncluster_c*2*ntime_c);
 double *backwar_c = dvector(1, nregion_c*ncluster_c*2*ntime_c);
 double **probability_c = dmatrix(1, nregion_c, 1, ncluster_c);
 double **probability_temp = dmatrix(1, nregion_c, 1, ncluster_c);
 long n = 0, *H_c = lvector(1, nregion_c*ntime_c);
 double b_temp, d, d_p, d_lambda;
 double *B = dvector(1, ncluster_c*4), *Lambda_temp=dvector(1, ncluster_c*2*nmarker_c);
 
 // Initialization
 for(i=1; i<=nregion_c*ntime_c; i++)
    {
	 H_c[i] = H[i-1];
	}
 
 for(i=1; i<= (nregion_c*ntime_c*nmarker_c); i++)
    {
	 observation_c[i] = observation[i-1];
	}
 for(i=1; i<= (ncluster_c*4); i++)
    {
	 b_c[i] = b[i-1];
	}
 
  for(i=1; i<= (ncluster_c*2*nmarker_c); i++)
    {
	 lambda_c[i] = lambda[i-1];
	}
  for(i=1; i<= nregion_c; i++)
    {
	 cluster_c[i] = cluster[i-1];
	}
  for(i=1; i<= (nregion_c*ncluster_c*2*ntime_c); i++)
    {
	 forwar_c[i] = forwar[i-1];
	 backwar_c[i] = backwar[i-1];
	}

  for(i=1; i<= nregion_c; i++)
    {
        for(j=1; j<=ncluster_c; j++)
		   {
		    probability_c[i][j] = probability[i+(j-1)*nregion_c];
			probability_temp[i][j] = probability_c[i][j];
		   }
    }
  for(i=1; i<= ncluster_c; i++)
    {
	 pie_c[i] = pie[i-1];
	}
 
 forward(forwar_c, b_c, lambda_c, observation_c, nstate_c, nmarker_c, ntime_c, nregion_c, ncluster_c);
 backward(backwar_c, b_c, lambda_c, observation_c, nstate_c, nmarker_c, ntime_c, nregion_c, ncluster_c);
 cluster_matrix(nregion_c, ncluster_c, ntime_c, cluster_c, pie_c, forwar_c, probability_c);

 for(i=1; i<=maxiteration_c; i++)
    {
	 
	 
	 for(k=1; k<=ncluster_c*4; k++)
      {
	    B[k] = b_c[k];
	  }
	for(k=1; k<= (ncluster_c*2*nmarker_c); k++)
    {
	 Lambda_temp[k] = lambda_c[k];
	}

	    for(k = 1; k<=ncluster_c; k++)
	    {
		 for(j = 1; j<=2; j++)
		    {
			 b_temp = b_c[k+(j-1)*ncluster_c*2];
			 if((b_temp<0.001)||(b_temp>0.999))
			   {
			    Rprintf("change b to 0.5\n");
				b_c[k+(j-1)*ncluster_c*2] = 0.5;
			   }
			}
		 
		}
		
	   
	 forward(forwar_c, b_c, lambda_c, observation_c, nstate_c, nmarker_c, ntime_c, nregion_c, ncluster_c);
	 backward(backwar_c, b_c, lambda_c, observation_c, nstate_c, nmarker_c, ntime_c, nregion_c, ncluster_c);
	 HMM_Recursion(b_c, observation_c, lambda_c, cluster_c, forwar_c, backwar_c, nstate_c, ncluster_c, nmarker_c, ntime_c, nregion_c, H_c);


	   
	    for(j=1; j<=nregion_c*ntime_c; j++)
        {
	     H[j-1] = H_c[j];

	    }

	 
	 for(l=1; l<=ncluster_c; l++)
	    {
		 pie_r[(i-1)*ncluster_c+l-1] = pie_c[l];
		}
	 d_lambda = Edistance(Lambda_temp, lambda_c, (ncluster_c*2*nmarker_c)); 
	 d = Edistance(b_c, B, ncluster_c*4);
	 diff_r[i-1] = d;
	 log_likeli_r[i-1]=log_likelihood(nregion_c, ncluster_c, ntime_c, cluster_c, pie_c, forwar_c);
	 n++;

	 Rprintf("%d, %lf, %lf, %d \n", i, d, d_lambda, n);
	 for(j=1; j<=ncluster_c; j++)
	 {
	    Rprintf("%lf \n", b_c[j+2*ncluster_c]);
		Rprintf("%lf \n", pie_c[j]);
	 }
	 if((n==nstep_c)||(d<ndistance_c))
	   {
	    forward(forwar_c, b_c, lambda_c, observation_c, nstate_c, nmarker_c, ntime_c, nregion_c, ncluster_c);
		cluster_matrix(nregion_c, ncluster_c, ntime_c, cluster_c, pie_c, forwar_c, probability_c);
	    d_p = mEdistance(probability_c, probability_temp, nregion_c, ncluster_c);
		copy_matrix(probability_temp, probability_c, nregion_c, ncluster_c);
		Rprintf("difference of p is: %lf \n",d_p);
		if (d_p<0.0001)
		  {break;}
		n = 0;
	   }
	  
	}
 *stop=i-1;



 for(j = 1; j<=ncluster_c*2*nmarker_c; j++)
    {
	 lambda[j-1] = lambda_c[j];
	}
 
 for(j = 1; j<=ncluster_c*4; j++)
    {
	 b[j-1] = b_c[j];
	}
for(i=1; i<= nregion_c; i++)
    {
        for(j=1; j<=ncluster_c; j++)
		   {
		     probability[i+(j-1)*nregion_c] = probability_c[i][j];
		   }
    }

for(i=1; i<= nregion_c; i++)
   {
	 cluster[i-1] = cluster_c[i];
   }
 
 free_dvector(observation_c, 1, nregion_c*ntime_c*nmarker_c);
 free_dvector(b_c, 1, ncluster_c*4);
 free_dvector(pie_c, 1, ncluster_c);
 free_dvector(lambda_c, 1, ncluster_c*2*nmarker_c);
 free_lvector(cluster_c, 1, nregion_c);
 free_dvector(forwar_c, 1, nregion_c*ncluster_c*2*ntime_c);
 free_dvector(backwar_c, 1, nregion_c*ncluster_c*2*ntime_c);
 free_dmatrix(probability_c, 1, nregion_c, 1, ncluster_c);
 free_dmatrix(probability_temp, 1, nregion_c, 1, ncluster_c);
 free_dvector(B, 1, ncluster_c*4);
 free_lvector(H_c, 1, nregion_c*ntime_c);
}





#undef NRANSI
