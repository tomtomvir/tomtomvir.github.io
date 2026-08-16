#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>

#define N 2048
#define T 1.


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
fft(X_n,N);
for(n=0;n<N;n++){
 f_n[n] = (-N/2+n)*fs/N;
 if(n<N/2) f_n[n] = n*fs/N;
 if(n>N/2) f_n[n] = (-N + n)*fs/N;
if(f_n[n] < 30. && f_n[n] > -30.){
printf("%lf %lf\n",f_n[n],cabs(X_n[n]));}
}
	return 0;
}
