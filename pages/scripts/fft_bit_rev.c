#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>

#define N 2048
#define logN 11
#define T 1.

int rev(int n){
int i,nbin[logN];
for(i=0;i<logN;i++){
nbin[i] = n%2;n/=2;
}
int sum=0;
for(i=0;i<logN;i++){
sum+= (1<<i)*nbin[logN-i-1];
}
return sum;
}
void bit_rev(complex double* x){
complex double temp;

for(int i=0,j=0;i<N-1;i++){
	j = rev(i);
if(i<j){temp = x[i]; x[i] = x[j]; x[j] = temp;}
}
return;
}
void fft_bit_rev(complex double* x){
	bit_rev(x);
	int Ncurr;
	complex double W,Wtot,u;
	for(int i=1,j;i<logN+1;i++){
		Ncurr = (1 << i);
		W = cexp(-2.*M_PI*I/Ncurr);
		for(j=0;j<N;j+=Ncurr){
			Wtot = 1.;
			for(int k=0;k<Ncurr/2;k++){
			u = Wtot*x[j+k+Ncurr/2];
			x[j+k+Ncurr/2] = x[j+k] - u;
			x[j+k] = x[j+k] + u ;
			Wtot*=W;
			}
		}
	}
return;
}

void fft(complex double* x, int dim){
	complex double *odd, *even, temp;
	int k;
	if(dim == 1) return;
	else {
		 odd = ( complex double*) malloc( dim/2 * sizeof(complex double));
		 even = (complex double*) malloc( dim/2 * sizeof(complex double));
		 for(k=0;k<dim/2;k++){
			 odd[k] = x[2*k+1];
			 even[k] = x[2*k];
		 }
		 fft(odd,dim/2);fft(even,dim/2);
		 for(k=0;k<dim/2;k++){
			 temp = cexp(-2.*M_PI*I/dim*k);
			 x[k]=even[k] + temp*odd[k];

			 x[k+dim/2]=even[k] - temp*odd[k];
		 }
	}
free(odd); free(even);	
	return;
}

double f(double t){
	return exp(-t*t/T/T*100)*sin(2*M_PI*t*20);

}
double fill(double complex *x_t, double* t_t){
double dt = 2*T/N;	
int i;
for (i=0;i<N;i++){
x_t[i] = f(-T + dt*i); 
t_t[i] = -T + dt*i;
}
return dt;
}

int main(){
double complex x_t[N];
double t_t[N];
double complex X_n[N];
double f_n[N];
int n,k;
double complex sum;
double dt=fill(x_t,t_t);
double fs =  1./dt;
/*
for(n=0;n<N;n++){
	sum = 0;
	for(k=0;k<N;k++){
		sum+=x_t[k]*cexp(- 1.*2*M_PI*I*k*(-N/2+n)/N);
	}
X_n[n]=dt*sum; 
}
*/
for(n=0;n<N;n++){
X_n[n] = x_t[n];
}
fft_bit_rev(X_n);
for(n=0;n<N;n++){
 f_n[n] = (-N/2+n)*fs/N;
 if(n<N/2) f_n[n] = n*fs/N;
 if(n>N/2) f_n[n] = (-N + n)*fs/N;
if(f_n[n] < 30. && f_n[n] > -30.){
printf("%lf %lf\n",f_n[n],cabs(X_n[n]));}
}
	return 0;
}
