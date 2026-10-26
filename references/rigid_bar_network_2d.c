#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


#define	N						1600
#define	DOUBLE_N					((double)N)
#define	DIM							2
#define	DELTA(X,Y)			((X)==(Y)?1.0:0.0)
#define	X_COMP					0
#define 	Y_COMP 					1
#define	I_INDEX					0
#define	J_INDEX					1
#define	CUTOFF_SQRD					(1.00*1.00)
#define	CGSOLVER_TOL				1e-10

/* constant constants... :P*******************************************************/
#define	LENGTH				1.0
#define	HALF_LENGTH			0.5
#define	HUGE_INT				10000000
#define	UNI		((double)rand()/((double)RAND_MAX + 1.0))
#define	TWO_PI					6.28318530717958


/*globals*/
double L; 			/*	length of <square> box 	*/
double invL;			/* inverse length of box	*/
double V;			/*	volume			*/
double u;			/*	potential energy 		*/
double kinetic;		/*	kinetic energy		*/
double stress;
double leDisplacement;
double strain;
int contacts;

double *rx;		/*	x component of position	*/
double *ry;		/*	y component of position	*/

double *auxVec;

/************************ new bondnetwork nebz ******************************************************/
int *lookupBond;
double *nij;
double *restLength;
/*******************************************************************************************************/

int readNetwork(char *inFileName){
	int i,j,l,foo;
	double dummy,dx,dy,r;
	FILE *inFile;
	
	rx = (double *)malloc(sizeof(double)*N);
	ry = (double *)malloc(sizeof(double)*N);
	auxVec = (double *)malloc(sizeof(double)*N*DIM);
	
	inFile = fopen(inFileName,"rb");
	if ( inFile ){
		foo = fscanf(inFile,"%lf %lf %lf",&L, &strain, &dummy);
		invL = 1.0/L; V = L*L;
		leDisplacement = fmod(strain,1.0);
		if (leDisplacement < -0.5)
			leDisplacement += 1.0;
		else if (leDisplacement > 0.5)
			leDisplacement -= 1.0;
		
		for (i=0;i<N;i++)
			foo = fscanf(inFile,"%lf %lf %lf", &(rx[i]), &(ry[i]), &dummy);
		
		foo = fscanf(inFile,"%d %lf %lf", &contacts, &dummy, &dummy);
		
		printf("going to read %d contacts\n\n", contacts); 
		
		lookupBond = (int *)malloc(sizeof(int)*2*contacts);				// allocate lookup-bond table memory
		nij = (double *)malloc(sizeof(double)*contacts*DIM);
		restLength = (double *)malloc(sizeof(double)*contacts);
		
		for (l=0;l<contacts;l++){
			foo = fscanf(inFile,"%d %d %lf", &(lookupBond[2*l + I_INDEX]), &(lookupBond[2*l + J_INDEX]), &(restLength[l]) );
			i = lookupBond[2*l + I_INDEX];
			j = lookupBond[2*l + J_INDEX];
			dx = rx[j] - rx[i];
			dy = ry[j] - ry[i];
			// mess due to LEES EDWARDS periodic boundary conditions
			if ( dy >= HALF_LENGTH ){
				dy -= LENGTH;
				dx -= leDisplacement;
			}
			else if ( dy < -HALF_LENGTH ){
				dy += LENGTH;
				dx += leDisplacement;
			}
			if ( dx >= HALF_LENGTH )
				dx -= LENGTH;
			else if ( dx < -HALF_LENGTH )
				dx += LENGTH;
			// end of mess
			r = sqrt(dx*dx + dy*dy); //not real distance but doesn't matter since we want nij
			nij[DIM*l + X_COMP] = dx/r;
			nij[DIM*l + Y_COMP] = dy/r;
		}
		fclose(inFile);
		return 1;
	}
	printf("couldn't read %s\n",inFileName);
	return 0;
}

void saveNetwork(char *networkFileName){
	int i,j,l;
	FILE *networkFile;
	networkFile = fopen(networkFileName,"wb");
	fprintf(networkFile,"%.15g\t%.15g\t-1\n",L,strain);
	for (i=0;i<N;i++)
		fprintf(networkFile,"%.15g\t%.15g\t-1\n",rx[i],ry[i]);
	
	fprintf(networkFile,"%d\t-1\t-1\n",contacts);
	for (l=0;l<contacts;l++){
		i = lookupBond[2*l+I_INDEX];
		j = lookupBond[2*l+J_INDEX];
		fprintf(networkFile,"%d\t%d\t%.15g\n",i,j,restLength[l]);
	}
	fclose(networkFile);
	return;
}


void free_everything(){
	free(rx); free(ry);
	free(lookupBond);
	free(restLength);
	free(nij); free(auxVec);
	return;
}

void imposeStrainIncrement(double strainIncrement){
	int i;
	double temp;
	
	strain += strainIncrement;
	leDisplacement = fmod(strain,1.0);
	if (leDisplacement < -0.5)
		leDisplacement += 1.0;
	else if (leDisplacement > 0.5)
		leDisplacement -= 1.0;
	for (i=0;i<N;i++){
		temp = rx[i] + strainIncrement*ry[i]; 
		
		if (temp >= LENGTH)
			rx[i] = temp - LENGTH;
		else if (temp < 0.0)
			rx[i] = temp + LENGTH;
		else
			rx[i] = temp;
	}
	return;
}

void esOnVector(double *vec, double *res){
	int l,i,j;
	
	for (l=0;l<contacts;l++){
		i = lookupBond[2*l+I_INDEX];
		j = lookupBond[2*l+J_INDEX];
		res[l] =  	nij[DIM*l+X_COMP]*( vec[DIM*j + X_COMP] - vec[DIM*i + X_COMP] ) + 
				nij[DIM*l+Y_COMP]*( vec[DIM*j + Y_COMP] - vec[DIM*i + Y_COMP] );
	}
	return;
}

void esTransposeOnVector(double *vec, double *res){
	int i,j,l;

	for (i=0;i<N*DIM;i++)
		res[i] = 0.0;
	
	for (l=0;l<contacts;l++){
		i = lookupBond[2*l+I_INDEX];
		j = lookupBond[2*l+J_INDEX];
		res[DIM*i + X_COMP] -= vec[l]*nij[DIM*l + X_COMP]; //mind the minus
		res[DIM*j + X_COMP] += vec[l]*nij[DIM*l + X_COMP]; //mind the plus
		res[DIM*i + Y_COMP] -= vec[l]*nij[DIM*l + Y_COMP]; //mind the minus
		res[DIM*j + Y_COMP] += vec[l]*nij[DIM*l + Y_COMP]; //mind the plus
	}
	return;
}

//here <vec> and <res> are both in contact space
void enOnVector(double *vec, double *res){
	esTransposeOnVector(vec,auxVec); //auxVec is now global...
	esOnVector(auxVec,res);
	return;
}

//contact space contraction
double cdot(double *v1, double *v2){
	int l;
	double res = 0.0;
	for (l=0;l<contacts;l++)
		res += v1[l]*v2[l];
	return res;
}	

//solves contact-space problems
int cgSolver( double *x, double *b){
	double *r, *w, *z;
	double a,B,maxNorm,scale,sum;
	int i,j,k,t_0,l,q;
		
	t_0 = (int)time(0);
	
	scale = sqrt( cdot(b,b)/( (double)contacts ) );
	
	r = (double *)malloc(contacts*sizeof(double));
	w = (double *)malloc(contacts*sizeof(double));
	z = (double *)malloc(contacts*sizeof(double));
	
	enOnVector(x,r);
	for (i=0;i<contacts;i++){
		r[i] = b[i] - r[i];
		w[i] = -r[i];
	}
	
	enOnVector(w,z);
	a = cdot(r,w)/cdot(w,z);
	for (i=0;i<contacts;i++)
		x[i] += a*w[i];
	
	B = 0.0;
	for (k=0; k<HUGE_INT; k++){
		maxNorm = 0.0;
		for (i=0;i<contacts;i++){
			r[i] = r[i] - a*z[i];
			if (r[i]*r[i] > maxNorm)
				maxNorm = r[i]*r[i];
		}
		maxNorm = sqrt(maxNorm);
		fprintf(stderr, "\r");
		fprintf(stderr, "maxNorm/scale: %.10g",maxNorm/scale);
		
		if ( maxNorm/scale < CGSOLVER_TOL )
			break;
		//~ if ( maxNorm < CGSOLVER_TOL )
			//~ break;
		
		B = cdot(r,z)/cdot(w,z);
		for (i=0;i<contacts;i++)
			w[i] = -r[i] + B*w[i];
		
		enOnVector(w,z);
		a = cdot(r,w)/cdot(w,z);
		
		for (i=0;i<contacts;i++)
			x[i] += a*w[i];
	}
	//fprintf(stderr, "\n");
	free(r); free(w); free(z);
	//printf("took %d seconds to linSolve\n",(int)(time(0) - t_0));
	if (k==HUGE_INT){
		printf("couldn't solve linear equation\n");
		return 0;
	}
	else
		return 1;
}



int main(int n,char **inputStrings){
	int q;
	double coordination;
        char networkFileName[1024];
        
	srand( time(0) );
	
	coordination = 3.5;
	q = 0;
	sprintf(networkFileName,"floppy_network_N%d_z%.3f_s%.5d.dat",N,coordination,q);
	if ( readNetwork(networkFileName) ){
		printf("contacts %d\n",contacts);
		
		free_everything();
	}
	else
		printf("%s not about\n",networkFileName);
	
	return 0;
}
