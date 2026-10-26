#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>


#define	N						1600
#define	DOUBLE_N					((double)N)
#define	LARGE						0.7
#define	SMALL						0.5
#define	DIM							2
#define	TIME_STEP						0.02
#define	HALF_TIME_STEP				0.5*TIME_STEP
#define	DELTA(X,Y)			((X)==(Y)?1.0:0.0)
#define	MAX_NEBZ						32
#define	MAX_CELLS					(N)
#define	X_COMP					0
#define 	Y_COMP 					1
#define	CUTOFF_SQRD					(1.00*1.00)
#define	TAU_T						1.0 //for thermostating
#define	DR							0.20

/*fire shit */
#define	MAX_ITERATIONS		1000000
#define	RATIO_TOL					1e-7
#define	TOL							1e-10
#define	N_MIN						5
#define	F_INC						1.1
#define	F_DEC						0.5
#define	F_ALPHA						0.99
#define	ALPHA_START					0.1
#define	DELTA_T_MAX					(10.0*TIME_STEP)

/* constant constants... :P*******************************************************/
#define	LENGTH				1.0
#define	HALF_LENGTH			0.5

#define	UNI		((double)rand()/((double)RAND_MAX + 1.0))
#define	TWO_PI					6.28318530717958

void updateNebzLists();
void initializeSystem();
void calculateForces();

/*globals*/
double L; 			/*	length of <square> box 	*/
double packingFraction;	/*	 packing fraction		*/
double invL;			/* inverse length of box	*/
double V;			/*	volume			*/
double Vp;			/*	volume occupied by particles			*/
double u;			/*	potential energy 		*/
double kinetic;		/*	kinetic energy		*/
double stress;
double strain;
double pressure;
double T;
int contacts;
double typicalGrad;
double typicalContactForce;
double coordination;


double *rx, *ry;
double *fx, *fy;
double *px, *py;
double *size;


/******************************* for nebz list & cell subdivision *******************************/
int **nebz;
int *numOfNebz;
int *numInCell;
int **cells;
double cellCutOff,maxD,listCutOffSqrd;
int nebListCounter;
const int shift[4][2] = {{-1,0},{-1,1},{0,1},{1,1}}; //{yShift,xShift};
/************************************************************************************/

/*******************************for network*******************************************/
int **networkNebz; 
int *numOfNetworkNebz;
double **pairwiseDistance;
double **angle;
/************************************************************************************/

/******************************* for fire **********************************************/
double timeStep;
double alpha;
int numOfStepsSinceNegativePower;
/************************************************************************************/


void fixDrift(){
	int i;
	double Px,Py;
	Px = 0.0; Py = 0.0;
	for (i=0;i<N;i++){
		Px += px[i];
		Py += py[i];
	}
	Px = Px/DOUBLE_N; Py = Py/DOUBLE_N; 
	for (i=0;i<N;i++){
		px[i] -= Px;
		py[i] -= Py;
	}
	return;
}

void initializeGrid(double _volumeFraction){
	int i,sqrtN,j;
	int *type;
	double space;
	
	type = (int *)malloc(sizeof(int)*N);
	for (i=0;i<N;i++)
		type[i] = 0;
	
	sqrtN=1;
	while (sqrtN*sqrtN<N)
		sqrtN++;
	space = 1.0/(double)sqrtN;
	for (i=0; i<N; i++){
		rx[i] = space*(double)(i%sqrtN);
		ry[i] = space*(double)(i/sqrtN);
		px[i] = 0.0; py[i] = 0.0; type[i] = 0;
	}
	for (i=0; i<N/2; i++){
		j = (int)(UNI*DOUBLE_N);
		while ( type[j] )
			j = (int)(UNI*DOUBLE_N);
		type[j] = 1;
	}
	
	Vp = 0.0;
	for (i=0;i<N;i++){
		if ( type[i] ) //large
			size[i] = LARGE;
		else
			size[i] = SMALL;
		Vp += 0.5*TWO_PI*size[i]*size[i] ;
	}
	V = Vp/_volumeFraction;
	packingFraction = _volumeFraction;
	L = sqrt(V); invL = 1.0/L;
	initializeSystem();
	free(type);
	return;
}

void initializeSystem(){
	int i;
	double maxSize = 0.0;
	for (i=0;i<N;i++)
		if (size[i] > maxSize)
			maxSize = size[i];
	cellCutOff = 2.0*maxSize*sqrt(CUTOFF_SQRD) + DR;
	listCutOffSqrd = cellCutOff*cellCutOff;
	Vp = 0.0;
	for (i=0;i<N;i++)
		Vp += 0.5*TWO_PI*size[i]*size[i] ;
	packingFraction = Vp/V;
	updateNebzLists();
	calculateForces();
	return;
}

void saveSnapShot(char *fileName){
	int i;
	FILE *file;
	file = fopen(fileName,"wb");
	fprintf(file,"%.15g\t0.0\t%10g\n",L,pressure);
	for (i=0; i<N; i++)
		fprintf(file,"%.15g\t%.15g\t%.15g\n",rx[i],ry[i],size[i]);
	fclose(file);
	return;
}

int readSnapShot(char *fileName){
	int i,k,foo;
	double dummy;
	FILE *inFile;
	inFile = fopen(fileName, "rb");
	if ( inFile ){
		foo = fscanf(inFile, "%lf",&L);
		V = L*L; invL = 1.0/L;
		foo = fscanf(inFile, "%lf",&dummy);
		foo = fscanf(inFile, "%lf",&dummy);
		for (i=0;i<N;i++)
			foo = fscanf(inFile, "%lf %lf %lf",&(rx[i]),&(ry[i]),&(size[i]));
		fclose(inFile);
		strain = 0.0;
		initializeSystem();
		return 1;
	}
	return 0;
}

void updateNebzLists(){
	int i,j,current,mSqrd,a,b,k,l,numHere,numThere,target,w,m;
	int ii,jj;
	double dx,dy,r2;
	double includeSqrd,invCellSize;
	
	maxD = 0.0;
	nebListCounter = 0;
	
	m =(int)(L/cellCutOff);
	mSqrd = m*m; 		//the length of cells and numInCell...
	invCellSize = (double)m;
	
	//set number in each cell to zero
	for (i=0; i<mSqrd; i++)
		numInCell[i] = 0;
	
	/*********cells sorting********/
	for (i=0;i<N;i++){
		current = m*(int)(invCellSize*ry[i]) + (int)(invCellSize*rx[i]);
		cells[current][numInCell[current]] = i;
		numInCell[current]++;
	}

	//set all neb lists to zero
	for (i=0; i<N; i++)
		numOfNebz[i] = 0;
	
	for (i=0;i<m;i++){
		for (j=0;j<m;j++){
			current = m*i + j;
			numHere = numInCell[current];
			//first check interactions within a cell
			for (k=0;k<numHere-1;k++){
				for (l=k+1;l<numHere;l++){
					ii = cells[current][k];
					jj = cells[current][l];
					dx = rx[jj] - rx[ii];
					dy = ry[jj] - ry[ii];
					r2 = L*L*( dx*dx + dy*dy );
					if (r2 < listCutOffSqrd){
						nebz[ii][numOfNebz[ii]] = jj;
						nebz[jj][numOfNebz[jj]] = ii;
						numOfNebz[ii]++;
						numOfNebz[jj]++;
					}
				}
			}
			
			//now there are 4 nebz that need to be checked.
			for (w=0;w<4;w++){
				a = i + shift[w][0];
				if (a<0) a = m-1;
				else if (a==m) a = 0;
				b = j + shift[w][1];
				if (b==m) b = 0;
				target = m*a + b;
				numThere = numInCell[target];
				for (k=0;k<numHere;k++){
					for (l=0;l<numThere;l++){
						ii = cells[current][k];
						jj = cells[target][l];
						dx = rx[jj] - rx[ii];
						dy = ry[jj] - ry[ii];
						// mess due to periodic boundary conditions 
						if ( dx >= HALF_LENGTH )
							dx -= LENGTH;
						else if ( dx < -HALF_LENGTH )
							dx += LENGTH;
						if ( dy >= HALF_LENGTH )
							dy -= LENGTH;
						else if ( dy < -HALF_LENGTH )
							dy += LENGTH;
						// end of mess 
						r2 = L*L*( dx*dx + dy*dy );
						if (r2 < listCutOffSqrd){
							nebz[ii][numOfNebz[ii]] = jj;
							nebz[jj][numOfNebz[jj]] = ii;
							numOfNebz[ii]++;
							numOfNebz[jj]++;
						}
					}
				}
			}
		}
	}
	return;
}

void calculateForces(){
	int i,j,m,k;
	double g,dx,dy,rxi,ryi,r2;
	double r,sizeI,sizeIJ,temp;
	
	//set to zero.
	u = 0.0; pressure = 0.0; contacts = 0; typicalContactForce = 0.0;
	for (i=0;i<N;i++){
		fx[i] = 0.0; fy[i] = 0.0;
	}
	
	for (i=0;i<N-1;i++){
		m = numOfNebz[i];
		for (k=0;k<m;k++){
			j = nebz[i][k];
			if (j > i){
				dx = rx[j] - rx[i];
				dy = ry[j] - ry[i];
				// mess due to periodic boundary conditions 
				if ( dy >= HALF_LENGTH )
					dy -= LENGTH;
				else if ( dy < -HALF_LENGTH )
					dy += LENGTH;
				if ( dx >= HALF_LENGTH )
					dx -= LENGTH;
				else if ( dx < -HALF_LENGTH )
					dx += LENGTH;
				// end of mess 
				sizeIJ = size[i] + size[j];
				//invSizeSqrd = 1.0/(sizeIJ*sizeIJ);
				r2 = L*L*( dx*dx + dy*dy );
				if (r2 < sizeIJ*sizeIJ){
					contacts++;
					r = sqrt(r2);
					g = (sizeIJ - r)/r; //this must be positive
					temp = L*dx*g;
					fx[j] += temp;
					fx[i] -= temp;
					temp = L*dy*g;
					fy[j] += temp;
					fy[i] -= temp;
					stress += -temp*L*dx;
					pressure += g*r2;
					u += 0.5*(sizeIJ - r)*(sizeIJ - r);
					typicalContactForce += g*r; //purely repulsive!
				}
			}
		}
	}
	stress = stress/V;
	pressure = pressure/(V*(double)DIM);
	for (typicalGrad = 0.0,i=0;i<N;i++)
		typicalGrad += fx[i]*fx[i] + fy[i]*fy[i];
	typicalGrad = sqrt(typicalGrad/DOUBLE_N);
	typicalContactForce = sqrt(typicalContactForce/DOUBLE_N);
	return;
}

//runs, with thermostating.
void advanceTime(){
	int i;
	double temp,dx,dy,dz,inst_T,xi_T;
	
	for (i=0; i<N; i++){
		px[i] += HALF_TIME_STEP*fx[i];
		py[i] += HALF_TIME_STEP*fy[i];
		
		dy = TIME_STEP*py[i]; //real space displacement
		temp = ry[i]  + dy*invL;
		if (temp >= LENGTH)
			ry[i] = temp - LENGTH;
		else if (temp < 0.0)
			ry[i] = temp + LENGTH;
		else
			ry[i] = temp;
		
		dx = TIME_STEP*px[i];
		temp = rx[i] + dx*invL;
		if (temp >= LENGTH)
			rx[i] = temp - LENGTH;
		else if (temp < 0.0)
			rx[i] = temp + LENGTH;
		else
			rx[i] = temp;
		
		temp = dx*dx + dy*dy;
		if (temp > maxD)
			maxD = temp;
	}

	nebListCounter++;
	if ( 2.0*((double)nebListCounter)*sqrt(maxD) > DR )
		updateNebzLists();
	
	calculateForces();
	kinetic = 0.0;
	for (i=0; i<N; i++){
		px[i] += HALF_TIME_STEP*fx[i];
		py[i] += HALF_TIME_STEP*fy[i];
		kinetic += px[i]*px[i] + py[i]*py[i];
	}
	kinetic = 0.5*kinetic;
	
	/************************** thermostating part *************************************/
	inst_T = 2.0*kinetic/((double)DIM*DOUBLE_N);
	xi_T = sqrt(1.0 + TIME_STEP*(T/inst_T - 1.0)/TAU_T);
	for (i=0;i<N;i++){
		px[i] *= xi_T;
		py[i] *= xi_T;
	}
	/********************************************************************************/
	return;
}

void run(double duration){
	int k,steps;
	double t;
	//~ FILE *conFile = fopen("conserve.dat","wb");
	
	steps = (int)(duration/TIME_STEP);
	for (k=0;k<steps;k++){
		advanceTime();
		//if ( !(k%128) ){
			//~ fprintf(conFile,"%.10g\t%.10g\t%.10g\n",(double)k*TIME_STEP,kinetic,u);
			//printf("%.10g\t%.10g\t%.10g\n",(double)k*TIME_STEP,kinetic,u);
		//}
		if ( !(k%4096) )
			fixDrift();
	}
	//~ fclose(conFile);
	return;
}

void fireAdvanceTime(){
	int i;
	double temp,dx,dy;
	double power,mod_force,inv_mod_force,mod_velocity;
	
	for (i=0; i<N; i++){
		px[i] += 0.5*timeStep*fx[i];
		py[i] += 0.5*timeStep*fy[i];
		
		dy = timeStep*py[i]; //real space displacement
		temp = ry[i]  + dy*invL;
		if (temp >= LENGTH)
			ry[i] = temp - LENGTH;
		else if (temp < 0.0)
			ry[i] = temp + LENGTH;
		else
			ry[i] = temp;
		
		dx = timeStep*px[i]; //real space displacement
		temp = rx[i] + dx*invL;
		if (temp >= LENGTH)
			rx[i] = temp - LENGTH;
		else if (temp < 0.0)
			rx[i] = temp + LENGTH;
		else
			rx[i] = temp;
		
		temp = dx*dx + dy*dy;
		if (temp > maxD)
			maxD = temp;
	}

	nebListCounter++;
	if ( 2.0*((double)nebListCounter)*sqrt(maxD) > DR )
		updateNebzLists();
	
	calculateForces();
	kinetic = 0.0;
	for (i=0; i<N; i++){
		px[i] += 0.5*timeStep*fx[i];
		py[i] += 0.5*timeStep*fy[i];
		kinetic += px[i]*px[i] + py[i]*py[i];
	}
	kinetic = 0.5*kinetic;
	
	//calculate power
	power = 0.0;
	mod_force = 0.0;
	mod_velocity = 0.0;
	for (i=0;i<N;i++){
		power += fx[i]*px[i] + fy[i]*py[i];
		mod_force += fx[i]*fx[i] + fy[i]*fy[i];
		mod_velocity += px[i]*px[i] + py[i]*py[i];
	}
	inv_mod_force = 1.0/sqrt(mod_force);
	mod_velocity = sqrt(mod_velocity);
	
	for (i=0;i<N;i++){
		px[i] = (1.0 - alpha)*px[i] + alpha*mod_velocity*fx[i]*inv_mod_force;
		py[i] = (1.0 - alpha)*py[i] + alpha*mod_velocity*fy[i]*inv_mod_force;
	}
	
	if ( power > 0.0 ){
		if ( numOfStepsSinceNegativePower > N_MIN ){
			if ( F_INC*timeStep < DELTA_T_MAX )
				timeStep *= F_INC;
			else
				timeStep = DELTA_T_MAX;
			alpha *= F_ALPHA;
		}
		numOfStepsSinceNegativePower++;
	}
	else{
		numOfStepsSinceNegativePower = 0;
		timeStep *= F_DEC;
		for (i=0;i<N;i++){ //freeze the system
			px[i] = 0.0;
			py[i] = 0.0;
		}
		alpha = ALPHA_START;
	}
	return;
}


int fire(){
	int itr,t0;
	
	t0 = (int)time(0);
	
	numOfStepsSinceNegativePower = 0;
	alpha = ALPHA_START;
	calculateForces();
	timeStep = TIME_STEP;
	
	itr = 0;
	while ( itr < MAX_ITERATIONS && typicalContactForce > TOL && typicalGrad/typicalContactForce > RATIO_TOL ){
		fireAdvanceTime();
		itr++;
	}
	printf("iterations %d, pressure %.10g, ratio %.10g\n",itr,pressure,typicalGrad/typicalContactForce);
	printf("\n\nfire took %d seconds\n\n",(int)time(0) - t0);
	if ( itr == MAX_ITERATIONS ){
		printf("failed to minimize\n\n");
		return 0;
	}
	return 1;
}

/*************************************network stuff ***********************************************/

int countContacts(){
	int i,sum = 0;
	for (i=0;i<N;i++)
		sum += numOfNetworkNebz[i];
	return sum;
}

void makeInitialNetwork(){
	int i,j,k,m,l;
	double dx,dy,r2,sizeIJ,nx,temp;
	
	for (i=0;i<N;i++)
		numOfNetworkNebz[i] = 0;
	
	for (i=0;i<N-1;i++){
		m = numOfNebz[i];
		for (k=0;k<m;k++){
			j = nebz[i][k];
			if (j > i){
				dx = rx[j] - rx[i];
				dy = ry[j] - ry[i];
				// mess due to periodic boundary conditions 
				if ( dy >= HALF_LENGTH )
					dy -= LENGTH;
				else if ( dy < -HALF_LENGTH )
					dy += LENGTH;
				if ( dx >= HALF_LENGTH )
					dx -= LENGTH;
				else if ( dx < -HALF_LENGTH )
					dx += LENGTH;
				// end of mess 
				sizeIJ = size[i] + size[j];
				r2 = L*L*( dx*dx + dy*dy);
				if (r2 < sizeIJ*sizeIJ){
					networkNebz[i][numOfNetworkNebz[i]] = j;
					networkNebz[j][numOfNetworkNebz[j]] = i;
					
					nx = dx/sqrt(dx*dx + dy*dy);
					if ( dy > 0.0 ){
						angle[i][numOfNetworkNebz[i]] = acos(nx);
						angle[j][numOfNetworkNebz[j]] = 0.5*TWO_PI + acos(nx);
					}
					else{
						angle[i][numOfNetworkNebz[i]] = TWO_PI - acos(nx);
						angle[j][numOfNetworkNebz[j]] = 0.5*TWO_PI - acos(nx);
					}
						
					numOfNetworkNebz[i]++;
					numOfNetworkNebz[j]++;
				}
			}
		}
	}
	
	//sort angles ascending order
	for (i=0;i<N;i++){
		m = numOfNetworkNebz[i];
		for (k=0;k<m-1;k++){
			for (l=k+1;l<m;l++){
				if (angle[i][k]>angle[i][l]){
					//switch them around
					temp = angle[i][k];
					angle[i][k] = angle[i][l];
					angle[i][l] = temp;
					j = networkNebz[i][k];
					networkNebz[i][k] = networkNebz[i][l];
					networkNebz[i][l] = j;
				}
			}
		}
	}
	
	//choose a random node and output angles
	//~ i=(int)(UNI*DOUBLE_N);
	//~ printf("angles of %d: ",i);
	//~ m = numOfNetworkNebz[i];
	//~ for (k=0;k<m;k++)
		//~ printf("%g, ",angle[i][k]);
	//~ printf("\n");
	
	//here angles are all sorted
	return;
}

//~ void removeBond(){
	//~ int i,maxSum,minDifference;
	//~ int m,k,j,iBond,jBond;
	//~ int maxNebCoordination, nebWithMaxCoordination;
	
	//~ //first find the maxSum pair
	//~ maxSum = 0; 
	//~ for (i=0;i<N;i++){
		//~ //find neighbor with largest coordination
		//~ maxNebCoordination = 0; nebWithMaxCoordination = 0;
		//~ m = numOfNetworkNebz[i];
		//~ for (k=0;k<m;k++){
			//~ j = networkNebz[i][k];
			//~ if (j > i){ //go over every pair only once... 
				//~ if ( numOfNetworkNebz[j] > maxNebCoordination ){
					//~ maxNebCoordination = numOfNetworkNebz[j];
					//~ nebWithMaxCoordination = k;
				//~ }
			//~ }
		//~ }
		//~ if ( m + maxNebCoordination >  maxSum )
			//~ maxSum = m + maxNebCoordination;
	//~ }
	
	//~ //now find that pair that has maxSum, and the smallest difference.
	//~ minDifference = 1000000;
	//~ for (i=0;i<N;i++){
		//~ //find neighbor with largest coordination
		//~ m = numOfNetworkNebz[i];
		//~ for (k=0;k<m;k++){
			//~ j = networkNebz[i][k];
			//~ if (j > i){ //go over every pair only once... 
				//~ if ( numOfNetworkNebz[j] + m == maxSum ){
					//~ if ( abs(numOfNetworkNebz[j] - m) < minDifference ){
						//~ minDifference = abs(numOfNetworkNebz[j] - m);
						//~ iBond = i;
						//~ jBond = j;
					//~ }
				//~ }
			//~ }
		//~ }
	//~ }
	
	//~ //now remove bond
	//~ m = numOfNetworkNebz[iBond];
	//~ for (k=0;k<m;k++){
		//~ j = networkNebz[iBond][k];
		//~ if ( j == jBond ){
			//~ if ( k < m-1 ) //not the last in the list, replace with last
				//~ networkNebz[iBond][k] = networkNebz[iBond][m-1];
			//~ numOfNetworkNebz[iBond]--; //removed last
			//~ break;
		//~ }
	//~ }
	//~ m = numOfNetworkNebz[jBond];
	//~ for (k=0;k<m;k++){
		//~ i = networkNebz[jBond][k];
		//~ if ( i == iBond ){
			//~ if ( k < m-1 ) //not the last in the list, replace with last
				//~ networkNebz[jBond][k] = networkNebz[jBond][m-1];
			//~ numOfNetworkNebz[jBond]--; //removed last
			//~ break;
		//~ }
	//~ }
	//~ return;
//~ }

double remainingAngles(int i, int j){
	double remainingI,remainingJ;
	int k,l,m;
	
	//first find i in j's list
	m = numOfNetworkNebz[j]; k = 0;
	while (k<m){
		if ( networkNebz[j][k] == i )
			break;
		k++;
	}
	//k is the location of i in j's list
	//end cases:
	if ( k==0 )
		remainingJ = TWO_PI - angle[j][m-1] + angle[j][1];
	else if ( k==m-1 )
		remainingJ = TWO_PI - angle[j][m-2] + angle[j][0];
	else
		remainingJ = angle[j][k+1] - angle[j][k-1];
	
	//now find j in i's list
	m = numOfNetworkNebz[i]; k = 0;
	while (k<m){
		if ( networkNebz[i][k] == j )
			break;
		k++;
	}
	//k is the location of j in i's list
	//end cases:
	if ( k==0 )
		remainingI = TWO_PI - angle[i][m-1] + angle[i][1];
	else if ( k==m-1 )
		remainingI = TWO_PI - angle[i][m-2] + angle[i][0];
	else
		remainingI = angle[i][k+1] - angle[i][k-1];
	
	//return the larger of the two.
	if ( remainingI > remainingJ )
		return remainingI;
	else
		return remainingJ;
}

void removeBond(){
	int i,maxSum,minDifference;
	int m,k,j,l,iBond,jBond;
	double minRemainingAngles,temp;
	
	//first find the maxSum pair
	maxSum = 0; 
	for (i=0;i<N;i++){
		m = numOfNetworkNebz[i];
		for (k=0;k<m;k++){
			j = networkNebz[i][k];
			if (j>i){
				if ( m + numOfNetworkNebz[j] > maxSum )
					maxSum = m + numOfNetworkNebz[j];
			}
		}
	}
	
	//now find that pair that has maxSum, and the smallest difference.
	minDifference = 1000000;
	for (i=0;i<N;i++){
		//find neighbor with largest coordination
		m = numOfNetworkNebz[i];
		for (k=0;k<m;k++){
			j = networkNebz[i][k];
			if (j > i){ //go over every pair only once... 
				if ( numOfNetworkNebz[j] + m == maxSum ){
					if ( abs(numOfNetworkNebz[j] - m) < minDifference )
						minDifference = abs(numOfNetworkNebz[j] - m);
				}
			}
		}
	}
	
	minRemainingAngles = 1e20;
	for (i=0;i<N;i++){
		//find neighbor with largest coordination
		m = numOfNetworkNebz[i];
		for (k=0;k<m;k++){
			j = networkNebz[i][k];
			if (j > i){ //go over every pair only once... 
				if ( ( numOfNetworkNebz[j] + m == maxSum ) && ( abs(numOfNetworkNebz[j] - m) == minDifference ) ){
					temp = remainingAngles(i,j);
					if ( temp < minRemainingAngles ){
						minRemainingAngles = temp;
						iBond = i;
						jBond = j;
					}
				}
			}
		}
	}
	
	m = numOfNetworkNebz[iBond]; l = 0;
	while ( l < m ){
		if ( networkNebz[iBond][l] == jBond )
			break;
		l++;
	}
	//remove the k'th, and shift the rest back one place
	for (k=l;k<m-1;k++){
		networkNebz[iBond][k] = networkNebz[iBond][k+1];
		angle[iBond][k] = angle[iBond][k+1];
	}
	numOfNetworkNebz[iBond]--;
	
	//now remove bond
	m = numOfNetworkNebz[jBond]; l = 0;
	while ( l < m ){
		if ( networkNebz[jBond][l] == iBond )
			break;
		l++;
	}
	//remove the k'th, and shift the rest back one place
	for (k=l;k<m-1;k++){
		networkNebz[jBond][k] = networkNebz[jBond][k+1];
		angle[jBond][k] = angle[jBond][k+1];
	}
	numOfNetworkNebz[jBond]--;
	return;
}

void makeNetwork(double z){
	do{		
		removeBond();
		coordination = (double)countContacts()/DOUBLE_N;
		//printf("have %g, target %g\n",coordination,z);
	}while ( coordination > z ); // && fabs(coordination - z) > 1.0/DOUBLE_N );
	return;
}

void allocate_memory(){
	int i;
	
	rx = (double *)malloc(sizeof(double)*N);
	ry = (double *)malloc(sizeof(double)*N);
	px = (double *)malloc(sizeof(double)*N);
	py = (double *)malloc(sizeof(double)*N);
	fx = (double *)malloc(sizeof(double)*N);
	fy = (double *)malloc(sizeof(double)*N);
	size = (double *)malloc(sizeof(double)*N);
	
	numOfNebz = (int *)malloc(sizeof(int)*N);
	nebz = (int **)malloc(sizeof(int *)*N);
	for (i=0;i<N;i++)
		nebz[i] = (int *)malloc(sizeof(int)*MAX_NEBZ);
	
	numInCell = (int *)malloc(sizeof(int)*MAX_CELLS);
	cells = (int **)malloc(sizeof(int *)*MAX_CELLS);
	for (i=0;i<MAX_CELLS;i++)
		cells[i] = (int *)malloc(sizeof(int *)*MAX_NEBZ);
	
	numOfNetworkNebz = (int *)malloc(sizeof(int)*N);
	pairwiseDistance = (double **)malloc(sizeof(double *)*N);
	networkNebz = (int **)malloc(sizeof(int *)*N);
	angle = (double **)malloc(sizeof(double *)*N);
	for (i=0;i<N;i++){
		networkNebz[i] = (int *)malloc(sizeof(int)*MAX_NEBZ);
		pairwiseDistance[i] = (double *)malloc(sizeof(double)*MAX_NEBZ);
		angle[i] = (double *)malloc(sizeof(double)*MAX_NEBZ);
	}	
	return;
}

void free_everything(){
	int i;
	free(rx); free(ry);
	free(px); free(py);
	free(fx); free(fy);
	free(size);
	
	for (i=0;i<N;i++)
		free(nebz[i]);
	free(nebz);
	free(numOfNebz);
	
	for (i=0;i<MAX_CELLS;i++)
		free(cells[i]);
	free(cells);
	free(numInCell);
	
	for (i=0;i<N;i++){
		free(networkNebz[i]);
		free(angle[i]);
		free(pairwiseDistance[i]);
	}
	free(networkNebz);
	free(numOfNetworkNebz);
	free(angle);
	return;
}

void saveNetworkWithLengths(char *networkFileName){
	FILE *networkFile;
	double dx,dy,r;
	int i,j,m,k,counter;
	
	counter = 0;	
	L = sqrt(DOUBLE_N);
	
	for (i=0;i<N;i++){
		m = numOfNetworkNebz[i];
		for (k=0;k<m;k++){
			j = networkNebz[i][k];
			dx = rx[j] - rx[i];
			dy = ry[j] - ry[i];
			// mess due to periodic boundary conditions 
			if ( dy >= HALF_LENGTH )
				dy -= LENGTH;
			else if ( dy < -HALF_LENGTH )
				dy += LENGTH;
			if ( dx >= HALF_LENGTH )
				dx -= LENGTH;
			else if ( dx < -HALF_LENGTH )
				dx += LENGTH;
			// end of mess 
			r = sqrt(dx*dx + dy*dy);
			pairwiseDistance[i][k] = r;
			if ( i < j )
				counter++;
		}
	}
	
	printf("coordination:%.10g\n",(double)countContacts()/DOUBLE_N );
	networkFile = fopen(networkFileName,"wb");
	fprintf(networkFile,"%.15g\t0.0\t-1\n",L);
	for (i=0;i<N;i++)
		fprintf(networkFile,"%.15g\t%.15g\t-1\n",rx[i],ry[i]);
	
	fprintf(networkFile,"%d\t-1\t-1\n",counter);
	for (i=0;i<N;i++){
		for (k=0;k<numOfNetworkNebz[i];k++){
			if ( networkNebz[i][k] > i ){
				fprintf(networkFile,"%d\t%d\t%.15g\n",i,networkNebz[i][k],L*pairwiseDistance[i][k]);
				counter++;
			}
		}
	}
	fclose(networkFile);
	return;
}


void saveNetworkWithLengthsAndAngles(char *networkFileName){
	FILE *networkFile;
	double dx,dy,r,theta_ijk;
	int i,j,m,k,counter;
	double temp_sum=0.0;
	
	counter = 0;
	L = sqrt(DOUBLE_N);
	
	for (i=0;i<N;i++){
		m = numOfNetworkNebz[i];
		for (k=0;k<m;k++){
			j = networkNebz[i][k];
			dx = rx[j] - rx[i];
			dy = ry[j] - ry[i];
			// mess due to periodic boundary conditions 
			if ( dy >= HALF_LENGTH )
				dy -= LENGTH;
			else if ( dy < -HALF_LENGTH )
				dy += LENGTH;
			if ( dx >= HALF_LENGTH )
				dx -= LENGTH;
			else if ( dx < -HALF_LENGTH )
				dx += LENGTH;
			// end of mess 
			r = L*sqrt(dx*dx + dy*dy);
			pairwiseDistance[i][k] = r;
			if ( i < j )
				counter++;
		}
	}
	printf("found %d contacts\n",counter);
	printf("coordination:%.10g\n",(double)countContacts()/DOUBLE_N );
	networkFile = fopen(networkFileName,"wb");
	fprintf(networkFile,"%.15g\t0.0\t-1\t-1\n",L);
	for (i=0;i<N;i++)
		fprintf(networkFile,"%.15g\t%.15g\t-1\t-1\n",rx[i],ry[i]);
	
	fprintf(networkFile,"%d\t-1\t-1\t-1\n",counter);
	for (i=0;i<N;i++){
		for (k=0;k<numOfNetworkNebz[i];k++){
			if ( networkNebz[i][k] > i )
				fprintf(networkFile,"%d\t%d\t%.15g\t-1\n",i,networkNebz[i][k],pairwiseDistance[i][k]);
		}
	}
	counter = 0;
	for (i=0;i<N;i++){
		m = numOfNetworkNebz[i];
		//write first triple
		theta_ijk = TWO_PI - angle[i][m-1] + angle[i][0];
		fprintf(networkFile,"%d\t%d\t%d\t%.15g\n", networkNebz[i][m-1], i, networkNebz[i][0], theta_ijk);
		counter++;
		for (k=0;k<m-1;k++){ 
			theta_ijk = angle[i][k+1]-angle[i][k];
			fprintf(networkFile,"%d\t%d\t%d\t%.15g\n", networkNebz[i][k], i, networkNebz[i][k+1], theta_ijk);
			counter++;
		}
	}
	printf("angles: %d\n", counter);
	fclose(networkFile);
	return;
}


int main(int n,char **inputStrings){
        int q,serial,bonds;
	char networkFileName[1024];
	double target_coordination = 3.5;

	/* Usage: ./make_random_2d_networks [target_coordination] [serial]
	 * Both are optional; defaults are z=3.5, serial=0. Passing them lets
	 * you sweep z (and keep multiple samples per z distinct via serial)
	 * from a shell loop, without editing/recompiling this file each time,
	 * e.g.:
	 *   for z in 2.4 3.0 3.5 3.9; do ./make_random_2d_networks $z 0; done */
	if ( n >= 2 ) target_coordination = atof(inputStrings[1]);
	serial = (n >= 3) ? atoi(inputStrings[2]) : 0;
	if ( n < 2 )
		printf("No target_coordination given -- using default z = %.3f\n",
		       target_coordination);

	allocate_memory();

        /* Seed with time XOR pid, not time alone: if this binary is
         * launched several times in a loop (e.g. sweeping z), repeated
         * calls within the same wall-clock second would otherwise share
         * an RNG seed and produce the *same* initial packing/relaxation
         * for every z, silently destroying independence across the
         * z-sweep ensemble. */
        srand( (unsigned int)(time(NULL) ^ (unsigned int)getpid()) );
	for (q=0;q<1;q++){
		initializeGrid(1.0); //1.0 is the initial packing fraction
		T = 1.0;
		run(20.0); //just mix it up a bit
		if ( fire() ){
			makeInitialNetwork();
			bonds = countContacts();
			printf("initial network has %d bonds, coordination %g\n",bonds, (double)bonds/DOUBLE_N);
			makeNetwork(target_coordination);
			sprintf(networkFileName,"floppy_network_N%d_z%.3f_s%.5d.dat",N,target_coordination,serial);
			saveNetworkWithLengths(networkFileName);
			//saveNetworkWithLengthsAndAngles(networkFileName);
			printf("saved %s\n",networkFileName);
		}
	}
	free_everything();
	return 0;
}
