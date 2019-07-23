#include <math.h>
#include <mex.h>

extern "C" {
    int AS_Compare(const void*,const void*);
    double AccurateSum(int, double *);
    double Newellf(double, double, double);
	double NewellF2(double, double, double);
    void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
}

int AS_Compare(const void* px,const void* py)
{
  double x=fabs(*((const double *)px));
  double y=fabs(*((const double *)py));
  if(x<y) return 1;
  if(x>y) return -1;
  return 0;
}

double AccurateSum(int n, double *arr) 
{
	qsort(arr,n,sizeof(double),AS_Compare);
	
	double sum,corr,y,u,t,v,z,x,tmp;
	sum=arr[0]; corr=0;
	int i;
	for(i=1;i<n;i++) {	
		x=arr[i];
		y=corr+x;
		tmp=y-corr;
		u=x-tmp;
		t=y+sum;
		tmp=t-sum;
		v=y-tmp;
		z=u+v;
		sum=t+z;
		tmp=sum-t;
		corr=z-tmp;
	}
	return sum;
}

double Newellf(double x, double y, double z)
{
	x = fabs(x); double xsq = x * x;
	y = fabs(y); double ysq = y * y;
	z = fabs(z); double zsq = z * z;
	
	double R = xsq + ysq + zsq;
	if (R <= 0.0) return 0.0;
	else R = sqrt(R);
    
    double piece[8];
	int piececount = 0;

  	if( z>0. ) { 
		double temp1, temp2, temp3;
	    piece[piececount++] = 2*(2*xsq-ysq-zsq)*R;
		if( (temp1 = x*y*z) > 0. ) 
			piece[piececount++] = -12*temp1*atan2(y*z,x*R);
		if( y > 0. && ( temp2 = xsq + zsq ) > 0. ) {
			double dummy = log(((y+R)*(y+R))/temp2);
            piece[piececount++] = 3*y*zsq*dummy;
			piece[piececount++] = -3*y*xsq*dummy;
		}
		if( ( temp3 = xsq + ysq ) > 0. ) {
			double dummy = log(((z+R)*(z+R))/temp3);
			piece[piececount++] = 3*z*ysq*dummy;
			piece[piececount++] = -3*z*xsq*dummy;
		}
	} 
	else {
    	if( x == y ) {
			const double K = 2 * sqrt(2.0f) -6 * log(1+sqrt(2.0f));
      		piece[piececount++] = K*xsq*x;
		} 
		else {
			piece[piececount++] = 2*(2*xsq-ysq)*R;
			if(y>0. && x>0.)
				piece[piececount++] = -6*y*xsq*log((y+R)/x);
    	}
	}	

    return AccurateSum(piececount,piece)/12.;
}


double NewellF2(double x, double y, double z)
{
    double F = Newellf(x,y,z) - Newellf(x,0,z) - Newellf(x,y,0) + Newellf(x,0,0);
    return F;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 10) {
		mexPrintf("Improper Syntax\n");
		return;
	}
	
    double *A = mxGetPr(prhs[1]);
    double *B = mxGetPr(prhs[2]);
    double *C = mxGetPr(prhs[3]);
    int Nx = (int)*A;
    int Ny = (int)*B;
    int Nz = (int)*C;
    A = mxGetPr(prhs[4]);
    B = mxGetPr(prhs[5]);
    C = mxGetPr(prhs[6]);
    double dx = *A;
    double dy = *B;
    double dz = *C;
    A = mxGetPr(prhs[7]);
    B = mxGetPr(prhs[8]);
    C = mxGetPr(prhs[9]);
    double Xo = *A;
    double Yo = *B;
    double Zo = *C;
    
    mwSize M = mxGetNumberOfDimensions(prhs[0]);
    mwSize S = mxGetM(prhs[0]);
    const mwSize *D = mxGetDimensions(prhs[0]);
	int N    = mxGetNumberOfElements(prhs[0]);
	
    plhs[0] = (mxArray *)mxCreateNumericArray(M,D,mxGetClassID(prhs[0]),(mxComplexity)mxIsComplex(prhs[0]));
	double *out = mxGetPr(plhs[0]);

    int i,j,k;
    int k_pos, j_pos, i_pos;
    int k_inv, j_inv, i_inv;
    double k_mul, j_mul, i_mul;
    double k_min, j_min, i_min;
    
    for (k = 0; k < Nz/2; k++) {
        k_pos = k*Nx*Ny;
        k_inv = (Nz-k-1)*Nx*Ny;
        k_mul = Zo+k*dz;;
        k_min = Zo-(k+1)*dz;
        for (j = 0; j < Ny/2; j++) {
            j_pos = j*Nx;
            j_inv = (Ny-j-1)*Nx;
            j_mul = Yo+j*dy;
            j_min = Yo-(j+1)*dy;
            for (i = 0; i < Nx/2; i++) {
                i_pos = i;
                i_inv = Nx-i-1;
                i_mul = Xo+i*dx;
                i_min = Xo-(i+1)*dx;
                out[i_pos+j_pos+k_pos] = NewellF2(k_mul,i_mul,j_mul);
                out[i_inv+j_pos+k_pos] = NewellF2(k_mul,i_min,j_mul);
                out[i_pos+j_inv+k_pos] = NewellF2(k_mul,i_mul,j_min);
                out[i_inv+j_inv+k_pos] = NewellF2(k_mul,i_min,j_min);
                out[i_pos+j_pos+k_inv] = NewellF2(k_min,i_mul,j_mul);
                out[i_inv+j_pos+k_inv] = NewellF2(k_min,i_min,j_mul);
                out[i_pos+j_inv+k_inv] = NewellF2(k_min,i_mul,j_min);
                out[i_inv+j_inv+k_inv] = NewellF2(k_min,i_min,j_min);
            }
        }
    }

    return;
}
