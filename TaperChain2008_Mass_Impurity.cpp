//Rewriting of the Taperchain programs by Adam Sokolow 
//(velocity verlet algorithm taken from original code written by Jan Pfannes)
//In an attempt to upgrade the code...
//7/8/03
//Now actual date: 11/20/04
//Updated again 2.22.05
//5/19/07 Another revision for the breathing study, some preprogrammed driving functions added.
//10/13/08 adding precompression functionality for user and specification of initial compressions to individual grains
//20/06/15--- MP added the ability to put mass impurities in, and also added driving ability. See modified manual for more information.


//Reverting back to using the large particle as particle N-1, and small as 0	  7/9/03
 
/***************************************************
READ THIS PLEASE:


This code has been altered so that it will run from a file...
any parameter that you would like to change should be done so
through this built in file interface... otherwise this program
may produce undesired results.. or no results at all.

Also if you have little to no programming experience, keep in mind
that every little thing does matter, and if something gets accidentally
removed, the code may still compile and run, and the results you could get
could be completely believable, showing no obvious sign of error... 

Look... but be very careful. (and save a backup)

******************************************************/

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <iomanip>
//#include <malloc.h>
#include <stdio.h>
#include <time.h>




//***********************CONSTANTS/VARIABLES BETWEEN RUNS

 const double PI = 4 * atan(1.0);
 int nptles=20;                   // total number of particles
 int restart=0;                   // resume variable. here we will keep it set to 0, but it can be reset in the parameter file.



 //const bool Wall = true;
 int DEFAULTPRECISION = 10;
 bool leftWall=true;
 bool rightWall=true;
//const double rho = 3.2;				//SiC
//const double D = 0.00326603139013;

//TiAlV material constants are below...(just swap which is commented out)

//const double rho=4.42,//TiAlV            3.2 /* SiC (mg/mm^3) */, 
// D = 0.01206;//D=0.00326603139013 /* (mm^2/kN) */;

 double rho = 7.82; //Steel grains
//We need two density values now since I will be putting a mass impurity in the chain. I will use rubber walls and put a central rubber mass. One the same size as the steel
//grains and one 100 times larger.
 double rhoImpure = 1.00;
//Parameter for the host chain grain and impurity grain interactions- calculated for steel and rubber 
 double D1 = 19.0010245142;
//Parameter for the host chain grains (calculated for stainless steel)
 double D2 = 0.0070490285;
//Parameter for the impurity grain interactions (calculated for rubber-rubber): 
 double D3 = 38.48;
 int ImpGr = 19;  //This is the position of the impurity grain- we want it at the center (position 20) of 39 grain chain, but indexing starts at 0.
 double ImpRF = 100;  //This is the factor by which the radius is larger of the mass impurity.
 int NumImp = 0;  //This is the number of impurities in the chain. Impurity grains will be ImpGr + 1 --> ImpGr + 1 + NumImp (mp 06/15/15)
// 02.10.14- MP: I commented out the lines below and put my own D value in . 
//  double rho = 1.00, D = 0.1529; //rho(mg/mm^3), D(mm^/kN)
 
//const double rho = 7.780,//Stainless Steel
//D = 0.00825758;//D = 0.008327;//D = 0.0066619;
//D = 0.00196;


 double rlarge = 0.5;           // (radius of large ptle (mm))
 double q = 0;      // (tapering factor (%), q=0: monodisp.)	 

//changed below in code...

 double xn = 2.5;                   // (exponent in potential)
 double dt = 0.00001;               // (timestepwidth (musec))
 unsigned long long int nsteps = 500000;//12500000;// (# steps integration loop)

 double timeSpits = .5;  //in the output files, how frequent do we want output... 
 double ForceSpits = timeSpits;
 int realspit = int(timeSpits/dt);
 int realFspit = int(ForceSpits/dt);

 int window=0;
 int windowleft=0;
 int windowright=0;

//const double noiseDump = 1.0e-30;  // (max. unconsid. signal: dump)	//these aren't used
//const double noiseInter = 1.0e-15; // (max. unconsid. signal: int.)					  //same

//const double almostZero = 1.0e-50;	  //these aren't used

 double AMP, PER;
 double konstForce=0.0; //constant force applied to the first grain.
 double konstForceLast=0.0; //constant force applied to the last grain.
 double konstForceSym=0.0;
 double randForceVal=0.0; 
 double randForceValLast=0.0; //random force applied to the last grain.
 double randForceValFL=0.0; //random force applied symmetrically to first and last grains.
 int cursign=-1;
 int cursignSym=-1;
 int cursignR=-1;
 int asymLeft = 0;   //Controls the program flow--- if asymLeft = 0 it means no force applied to left edge of chain-- if == 1 then left asymmetric force.
 int asymRight = 0;
 int symLeftRight = 0; // If symLeftRight == 1, it means that a symmetric driving force is being applied to both ends of the chain.
 int konstLength = 0;  //length of driving intervals (in micro-s) for constant applied force. First is left edge of chain, second one is right edge of chain-- if doing asymmetric driving.
 int konstLengthLast = 0;
 int konstLengthSym = 0; //length of driving intervals (in micro-s) for constant applied force for the symmetrically driven chain.


//************* Function at end of code called 	addForcesomethignsomething... upgrades this
 double SmallInitialVelocity = 0.0;       // (initial v /small/ ptle (mm/musec))
 double LargeInitialVelocity = 0.00;//-1*0.149*0.001*1;//-.01      // initial v /large/ ptle (mm/musec))
//************

 double epsilon = 1; // ((1 - restitution factor) all ptles)
 int EP = 100; //change this variable... do not change the one above...

//**********Compression - set the initial overlap between the wall and the large particle
//from this a force balance is made... and all the particles are at rest (except minor vibrations
//due to inprecision in the decimals cause some crazy very small energy noise)
 double InitialOverlapLargeandWall=0.00;  
 double loadingForce=0;
//******************

//************dont alter these...
 double ChainLength;
 double ForceInChain;
//************************
 int count;
 double Force=0.0;
 double ForceLast=0.0;
 double ForceFirstLast=0.0;
 double storeForce=Force;
 double storeForceLast=ForceLast;
 double storeForceFL=ForceFirstLast;
//Given any single particle... it has the following:
struct particle
{
	double relativeLocation;
	double currentVelocity;
	double currentAccel;
	double radius;
	double currentKE;
	double mass;
	double absolutePosition;
};

//now we have an array,(the chain) of the above particles
particle* Chain;

double* deltas;//the overlaps
double* smalla;			  //the small-a in the calcs...
double* overbefore;  //stuff for the algorithm
double pot;
double* energy_loss;
double otime = 0.0;
double* forceBefore;

//these variables keep track of the max velocities... going in the direction of the initial impulse...
//(so ripples going back that are of larger velocity are ignored...)... this is easily alterable in the code
//the time variables correspond to what time that large velocity occurred.
double maxOverlapLargeParticle=0,maxOverlapSmallParticle=0, maxOLargeTime = 0, maxOSmallTime = 0;
double maxVelocityLargeParticle=0,maxVelocitySmallParticle=0, maxVLargeTime = 0, maxVSmallTime = 0;

std::ofstream KEouts, ETot, Vouts, XRelOuts, ForceAll, AppliedForce, AppliedForceLast, AppliedForceSym;  //file streams


void initializeAllParticles(); //set up the chain
void makeReadmeFile();		   //make a readme file
void ChainRelPositions();	   //set up more stuff
void CompressChain(double &force, double &length); //apply a compression

void velocityVerletStep();
void computeAccelerations();

double addForceDueToWave(double ttime);		//new addition... cerca 8/30/03
double addForceDueToWaveLast(double ttime);		//new addition... MP 20/06/15
double addForceDueToWaveFirstLast(double ttime);		//new addition... MP 20/06/15

void spitKEs(double ttime);	 //spit to the file the KEs of the particles at a given time
void spitVs(double ttime);	 //similar but velocities
void spitXs(double ttime);
void spitFapp(double ttime,double force);
void spitFappLast(double ttime,double forceLast);
void spitFappSym(double ttime,double forceSym);
double absolute(double val);
double sign(double val);

void loadFromFile();
void Resume(); //For resuming program from where it left off. MP 11/24/14.
void makeResumeFile();

//more of these output files can be added... but if its not needed, why output and slow the program down?

double recordTotE;

double TimeMin=0.0;
double TimeMax=nsteps*dt;   //We will set this after reading in parameters just to make sure everything is set up correctly. We specify where we want to start in
                            //the parameter file and the number of time steps, and then the program determines what the maximum time is from this. For now just leave
                            //TimeMin at 0.
std::string fileName="_"; 
bool FF=true, KK=true, VV=true, XX=true, TE=true, FA=false, FAL=false, FAS=false;

int topGrain=nptles;
int bottomGrain=1;

struct deltafunc
{
	int grain;
	double time;
	double addedVelocity;
        int howmany;   //MP-- How many times do you want to hit this thing??
	int interval;  //MP-- Interval between hits (in micro-s).
};

deltafunc* deltaVel;
deltafunc* deltaVelcopy;  //MP-- Used for writing to Readme File.
int deltafunclength=0;
int deltafunclengthcopy=0;  //MP-- Used for writing to Readme File.

struct drivingforces
{
	char type;
	double amplitude,frequency,phase,duration;
};

drivingforces* sourceforce;
drivingforces* sourceforceLast;
drivingforces* sourceforceFL;
int sourceforcelength =0;
int sourceforcelengthLast =0;
int sourceforcelengthFL =0;

int control;
double deltaW;
double deltaS;


//============================================================================================================================================================\\
//                                                                          MAIN PROGRAM                                                                      \\
//============================================================================================================================================================\\

int main()
{
	srand((unsigned)(time(0))); 

	std::string FileName;
	std::ostringstream *buffer;
	
	loadFromFile();
				
	TimeMax=nsteps*dt+TimeMin;   //Now we are setting the maximum time. If the program is not meant to be restarted, then TimeMin will be set to zero in loadFromFile();.
		
	ForceSpits = timeSpits;
	realspit = int(timeSpits/dt);
	realFspit = int(ForceSpits/dt);

	{
		std::cout<<"October 13, 2008 Version 1"<<std::endl;
		std::cout<<"number of grains: "<<nptles<<std::endl;
		std::cout<<"nsteps: "<<std::setprecision(8)<<nsteps<<std::endl;
		std::cout<<"dt: "<<dt<<std::endl;
		std::cout<<"Exponential: "<<xn<<std::endl;
		std::cout<<"q: "<<q<<std::endl;
		std::cout<<"epsilon: "<<epsilon<<std::endl;

		if(KK)
		{		 
			buffer = new std::ostringstream();
 			(*buffer) << "KE"<<fileName<<".dat";
			FileName = buffer->str();
			
			if(restart == 1){ KEouts.open(FileName.c_str(), std::ios::app); }
			else{ KEouts.open(FileName.c_str()); }
			KEouts.precision(DEFAULTPRECISION);
		 }

		 if(TE)
		 {
			 buffer = new std::ostringstream();
			 (*buffer)<<"ETot"<<fileName<<".dat";
			 FileName = buffer->str();

			 if(restart == 1){ ETot.open(FileName.c_str(), std::ios::app); }
			 else{ ETot.open(FileName.c_str()); }
			 ETot.precision(DEFAULTPRECISION);
		 }

		 if(VV)
		 {
			 buffer = new std::ostringstream();
			 (*buffer) << "VEL"<<fileName<< ".dat";
			 FileName = buffer->str();
			 
			 if(restart == 1){ Vouts.open(FileName.c_str(), std::ios::app); }
			 else{ Vouts.open(FileName.c_str()); } 
			 Vouts.precision(DEFAULTPRECISION);
		 }
		 
		 if(XX)
		 {
			 buffer = new std::ostringstream();
			 (*buffer) << "XRel"<<fileName<< ".dat";
			 FileName = buffer->str();
			 
			 if(restart == 1){ XRelOuts.open(FileName.c_str(), std::ios::app); }
			 else{ XRelOuts.open(FileName.c_str()); }
			 XRelOuts.precision(DEFAULTPRECISION);
		 }
		 
		 if(FF)
		 {
			 buffer = new std::ostringstream();
			 (*buffer) << "Force"<<fileName<< ".dat";
			 FileName = buffer->str();
			 
			 if(restart == 1){ ForceAll.open(FileName.c_str(), std::ios::app); }
			 else{ ForceAll.open(FileName.c_str()); }
			 ForceAll.precision(DEFAULTPRECISION);
		 }

		 if(FA)
		 {
			 buffer = new std::ostringstream();
			 (*buffer) << "DrivingForceLeft"<<fileName<< ".dat";
			 FileName = buffer->str();

			 if(restart == 1){ AppliedForce.open(FileName.c_str(), std::ios::app); }
 			 else{ AppliedForce.open(FileName.c_str()); }
			 AppliedForce.precision(DEFAULTPRECISION);
		 }
		 
		 if(FAL)  //right edge grain asymmetric driving
		 {
			 buffer = new std::ostringstream();
			 (*buffer) << "DrivingForceRight"<<fileName<< ".dat";
			 FileName = buffer->str();
			 
			 if(restart == 1){ AppliedForceLast.open(FileName.c_str(), std::ios::app); }
			 else{ AppliedForceLast.open(FileName.c_str()); } 
			 AppliedForceLast.precision(DEFAULTPRECISION);
		 }

		 if(FAS)  //symmetric driving
		 {
			 buffer = new std::ostringstream();
			 (*buffer) << "DrivingForceSymmetric"<<fileName<< ".dat";
			 FileName = buffer->str();

			 if(restart == 1){ AppliedForceSym.open(FileName.c_str(), std::ios::app); }
			 else{ AppliedForceSym.open(FileName.c_str()); } 
			 AppliedForceSym.precision(DEFAULTPRECISION);
		 }


		 initializeAllParticles();

		 ChainRelPositions();
		 
		 makeReadmeFile();
		 
		 computeAccelerations();

		 if(asymLeft==1)  //asymmetric driving on the left edge grain:
		 { 
			 Chain[nptles-1].currentAccel-=addForceDueToWave(otime);
			 std::cout<<Chain[nptles-1].currentAccel<<std::endl;
		 }


		 if(asymRight==1) //asymmetric driving on the right edge grain:
		 { 
			 Chain[0].currentAccel+=addForceDueToWaveLast(otime);
			 std::cout<<Chain[0].currentAccel<<std::endl;
		 }


		 if(symLeftRight==1) //symmetric driving on the both edge grains:
		 { 
			 Chain[0].currentAccel+=addForceDueToWaveFirstLast(otime)/Chain[0].mass;
			 Chain[nptles-1].currentAccel-=addForceDueToWaveFirstLast(otime)/Chain[nptles-1].mass;
			 std::cout<<Chain[0].currentAccel<<std::endl;
		 }


//Now that everything has been initialized as it normally would, now is the time to be resetting the necessary quantities for resuming where we left off.
//MP- 11/24/14. We will call on the "Resume" function, and in the resume function, we read in velocities and relative locations. Then we recompute the current grain
//accelerations. After this point, we are in a position to perform the simulation:

		if(restart==1){
			Resume();  	//This will reset all of the velocities and relative positions of the grains. After this, we use these positions and velocities to compute their 
					//Current accelerations. Then we can go through the Velocity Verlet steps.
			computeAccelerations();
		}


	
		 for(unsigned long long int	lcv = 0; lcv<nsteps; ++ lcv)
		 {
			 count = lcv;
			 otime = lcv * dt + TimeMin;   //MP-- added this is for the restart mode.

			 velocityVerletStep();
			 
			 for(int lcount=0;lcount<deltafunclength;++lcount)
			 {
				 control = deltaVel[lcount].howmany;
				 if(control == 1){
				 	if(deltaVel[lcount].time-otime<dt)
				 	{
						Chain[nptles-deltaVel[lcount].grain].currentVelocity+=deltaVel[lcount].addedVelocity ;
					 	deltaVel[lcount].addedVelocity=0;
					}
				 }
				 if(control == -1){
					deltaS = deltaVel[lcount].time;  //Time to start the "delta" perturbation-- MP
					deltaW = deltaS + deltaVel[lcount].interval; //Time to end the "delta" perturbation-- MP
					if((otime-deltaS<=dt)|((deltaS<otime)&&(otime<=deltaW))){
						Chain[nptles-deltaVel[lcount].grain].currentVelocity=deltaVel[lcount].addedVelocity ;
						//std::cout<<"applying delta to grain"<<deltaVel[lcount].grain<<std::endl;
					}
				}
			 }

			 if(lcv%realspit==0&&otime>=TimeMin&&otime<= TimeMax)
			 { 
				 if(KK||TE)
					 spitKEs(otime);
				 if(VV)
					 spitVs(otime);
				 if(FF||XX||FA||FAL||FAS)
					 spitXs(otime);
			 }
		
			 if((-1*Chain[0].currentVelocity)>maxVelocitySmallParticle)	//Added 6/7/03
			 { 
				 maxVelocitySmallParticle=(Chain[0].currentVelocity*-1);
				 maxVSmallTime = otime;
			 } //
			 if((-1*Chain[nptles-1].currentVelocity)>maxVelocityLargeParticle)		  //
			 { 
				 maxVelocityLargeParticle=(Chain[nptles-1].currentVelocity*-1);	
				 maxVLargeTime = otime;
			 }	//

			 if((overbefore[0])>maxOverlapSmallParticle)	//Added 6/7/03
			 { 
				 maxOverlapSmallParticle=(overbefore[0]);
				 maxOSmallTime = otime-dt;
			 } //
		
			 if(otime<100)
				 if((overbefore[nptles-1])>maxOverlapLargeParticle)		  //
				 { 
					 maxOverlapLargeParticle=(overbefore[nptles-1]);	
					 maxOLargeTime = otime-dt;
				 }	//
		 }

		 if(VV)
			 Vouts.close();
	
		 if(KK)	
			 KEouts.close();

		 if(TE)
			 ETot.close();  
	
		 if(XX)
			 XRelOuts.close();

		 if(FF)
			 ForceAll.close();

		 if(FA)
			 AppliedForce.close();

		 if(FAL)
			 AppliedForceLast.close();

		 if(FAS)
			 AppliedForceSym.close();

		
		 makeResumeFile();
		 
	}

	return 0;
}

void initializeAllParticles()
{
	Chain=(particle *) malloc(nptles*sizeof(particle));

	forceBefore = (double *)malloc((nptles+1)*sizeof(double));
	deltas= (double*)malloc((nptles+1)*sizeof(double));///the overlaps
	smalla= (double*)malloc((nptles+1)*sizeof(double));			  //the small-a in the calcs...
	overbefore = (double *)malloc((nptles+1)*sizeof(double));  //stuff for the algorithm
	energy_loss = (double*)malloc((nptles+1)*sizeof(double));

	Chain[nptles-1].radius = rlarge;
	
	maxVelocityLargeParticle=0;
	maxVelocitySmallParticle=0;
	maxVLargeTime = 0;
	maxVSmallTime = 0;

	maxOverlapLargeParticle=0;
	maxOverlapSmallParticle=0;
	maxOLargeTime = 0;
	maxOSmallTime = 0;

	for(int lcv = nptles-1;lcv>=0;lcv--)
	{
		Chain[lcv].relativeLocation=0;
		Chain[lcv].currentVelocity=0;
		Chain[lcv].currentAccel=0;
		Chain[lcv].currentKE=0;
		
		if(lcv!=nptles-1)
			Chain[lcv].radius = ((100.0-q)/100.0)*Chain[lcv+1].radius; 
 
		Chain[lcv].mass = Chain[lcv].radius*Chain[lcv].radius*Chain[lcv].radius*(4.0/3.0)*PI*rho;
 
		deltas[lcv]=0.0;
		smalla[lcv]=0.0;
		overbefore[lcv]=0.0;
		energy_loss[lcv]=0.0;
		forceBefore[lcv] = 0.0;
		Chain[lcv].absolutePosition = 0.0;
	}

//Now we determine the position ***is actually the index of the grain, where index = position -1 since indexing starts at 0*** of the first impurity grain:
	if(nptles % 2 == 0)  //chain has an even number of grains: then the number of impurities will also be an even number.
	{
		ImpGr = int(nptles/2 - NumImp/2);
	}
	else
	{
		ImpGr = int((nptles-1)/2) - int((NumImp-1)/2);
	}

//Now I want to reset the radius and mass of the impurity grains:
	if(NumImp>0)
	{
		for (int i3 = 0; i3 < NumImp; i3++)// 
		{	
			Chain[ImpGr+i3].radius = ImpRF*rlarge;
			Chain[ImpGr+i3].mass = Chain[ImpGr+i3].radius*Chain[ImpGr+i3].radius*Chain[ImpGr+i3].radius*(4.0/3.0)*PI*rhoImpure;
		}
	}

	Chain[nptles-1].currentVelocity = LargeInitialVelocity;
	Chain[0].currentVelocity = SmallInitialVelocity;
 
	deltas[nptles]=0.0;
	smalla[nptles]=0.0;
	overbefore[nptles]=0.0;
	energy_loss[nptles]=0.0;
 
	forceBefore[nptles]=0.0;
}

void makeReadmeFile()
{
	std::string FileName;
	std::ostringstream *buffer;
	std::ofstream out_file;

	buffer = new std::ostringstream();
	
    (*buffer) << "ReadMe"<<fileName<< ".dat";
    FileName = buffer->str();

	if(restart == 1){
		out_file.open(FileName.c_str(), std::ios::app);     //MP -- if we are in resume mode then we want to append to the files, not write over them.
	}
	else{
		out_file.open(FileName.c_str());
	} 
 	out_file.precision(DEFAULTPRECISION);


	if(restart == 1)
	{out_file<<"Resume mode. Resuming simulation at "<<TimeMin<<" micro-s"<<std::endl;}     //MP 11/24/14-- Added in for the resume mode. Output when we are picking up from.

	out_file<<"Number of Particles: "<<nptles<<std::endl;
	out_file<<"potential exponential: "<<xn<<std::endl;
	out_file<<"rho: "<<rho<<" (mg/mm^3)"<<std::endl;
        if(NumImp>0){
		out_file<<"rhoImpure: "<<rhoImpure<<" (mg/mm^3)"<<std::endl;
	}
        if(NumImp>0){
		out_file<<"D1, host-impurity interaction: "<<D1<<" (mm^2/kN)"<<std::endl;
	}
	out_file<<"D2, host-host interaction: "<<D2<<" (mm^2/kN)"<<std::endl;
	if(NumImp>1){	
		out_file<<"D3, impurity-impurity interaction: "<<D3<<" (mm^2/kN)"<<std::endl;
	}
	out_file<<"Tapering percent: "<<q<<std::endl;
	out_file<<"Restitution: "<<1-epsilon<<std::endl;
	out_file<<"Epsilon (1-w): "<<epsilon<<std::endl;
	out_file<<"Radius of Large Particle: "<<rlarge<<" (mm)"<<std::endl;
	
	out_file<<"Radii of all Particles: (mm)"<<std::endl;      //Recall that everything that is done in the program is swapped with the "real chain". i.e. left--> 
	for(int lcv=0; lcv<nptles;++lcv)
		out_file<<Chain[nptles-1-lcv].radius<<'\t';
	out_file<<std::endl;
	
	out_file<<"Masses of all Particles: (mg)"<<std::endl;
	for(int lcv2=0; lcv2<nptles;++lcv2)
		out_file<<Chain[nptles-1-lcv2].mass<<'\t';
	out_file<<std::endl;
	
	out_file<<"Total force prefactor (a) of all Particles: "<<std::endl;
	for(int lcv3=0; lcv3<nptles+1;++lcv3)
		out_file<<smalla[nptles-lcv3]<<'\t';
	out_file<<std::endl;
	
	out_file<<"dt: "<<dt<<"musec"<<std::endl;
	out_file<<"number of those steps: "<<nsteps<<std::endl;
	out_file<<"Length of run: "<<dt*nsteps<<" (musec)"<<std::endl;

	out_file<<"Initial Velocity of Large Particle: "<<-LargeInitialVelocity<<" (mm/musec)"<<std::endl;
	out_file<<"Initial Velocity of Small Particle: "<<-SmallInitialVelocity<<" (mm/musec)"<<std::endl;
		if(rightWall)
			out_file<<"Right Wall"<<std::endl;
		else
			out_file<<"No Right Wall"<<std::endl;
		if(leftWall)
			out_file<<"Left Wall"<<std::endl; 
		else
			out_file<<"No Left Wall"<<std::endl;

	if(deltafunclengthcopy>0){
		double dVel = 0.0;
		double dTime = 0.0;
		int dGrain = 0;	
		int dNum = 0;
		int dInt = 0;
		for(int lcv4=0; lcv4<deltafunclengthcopy; lcv4++){
			dVel = -deltaVelcopy[lcv4].addedVelocity;       //Negative because everything is negated when read into the program.
			dGrain = deltaVelcopy[lcv4].grain;
			dTime = double(deltaVelcopy[lcv4].time);
			dNum = deltaVelcopy[lcv4].howmany;
			dInt = deltaVelcopy[lcv4].interval;
			if(dNum == -1){
				out_file<<"Constant velocity perturbation applied to grain "<< dGrain <<" of magnitude "<< dVel <<" (mm/musec), beginning at t="<<dTime<<" musec, for a window of "<< dInt 				<<" musec."<<std::endl;
			}
			else{
				if(dNum == 1){
					out_file<<"Delta velocity perturbation applied to grain "<< dGrain <<" of magnitude "<< dVel <<" (mm/musec) at t="<< dTime <<" musecs."<<std::endl;
				}
				else{
					out_file<< dNum <<" delta velocity perturbations applied to grain "<< dGrain <<" of magnitude "<< dVel <<" (mm/musec). First perturbation at t="<< dTime 
					<<" musecs, and the remainder at "<< dInt <<" musec intervals."<<std::endl;
				}	
			}
		}
	}

	
	if(loadingForce>0)
		out_file<<"Initial Loading of Chain "<<loadingForce<<" (kN)"<<std::endl;
	

	if(symLeftRight == 1)
	{
		if(konstForceSym > 0.0)
		{
			out_file<<"Constant force applied symmetrically to edge grains "<<konstForceSym<<" (kN), window: "<<konstLengthSym<<" (micro-s)."<<std::endl;
		}
		for(int c=0;c<sourceforcelengthFL;++c)
		{
			switch(sourceforceFL[c].type)
			{
				case 's':
				{	
					out_file<<"Sine wave applied symmetrically to edge grains, amplitude: "<<sourceforceFL[c].amplitude<<" (kN), angular frequency: "
						<<sourceforceFL[c].frequency<<" (rad/micro-s), phase: "<<sourceforceFL[c].phase<<" (rad), window: "
						<<sourceforceFL[c].duration<<" musecs."<<std::endl;
					break;
				}
				case 'c':
				{
					out_file<<"Cosine wave applied symmetrically to edge grains, amplitude: "<<sourceforceFL[c].amplitude<<" (kN), angular frequency: "
						<<sourceforceFL[c].frequency<<" (rad/micro-s), phase: "<<sourceforceFL[c].phase<<" (rad), window: "
						<<sourceforceFL[c].duration<<" musecs."<<std::endl;
					break;
				}
				case 't':
				{
					out_file<<"Triangle wave applied symmetrically to edge grains, amplitude: "<<sourceforceFL[c].amplitude<<" (kN), angular frequency: "
						<<sourceforceFL[c].frequency<<" (rad/micro-s), phase: "<<sourceforceFL[c].phase<<" (rad), window: "
						<<sourceforceFL[c].duration<<" musecs."<<std::endl;
					break;
				}
				case 'w':
				{
					out_file<<"Saw wave applied symmetrically to edge grains, amplitude: "<<sourceforceFL[c].amplitude<<" (kN), angular frequency: "
						<<sourceforceFL[c].frequency<<" (rad/micro-s), phase: "<<sourceforceFL[c].phase<<" (rad), window: "
						<<sourceforceFL[c].duration<<" musecs."<<std::endl;
					break;
				}
				case 'q':
				{
					out_file<<"Square wave applied symmetrically to edge grains, amplitude: "<<sourceforceFL[c].amplitude<<" (kN), angular frequency: "
						<<sourceforceFL[c].frequency<<" (rad/micro-s), phase: "<<sourceforceFL[c].phase<<" (rad), window: "
						<<sourceforceFL[c].duration<<" musecs."<<std::endl;
					break;
				}
				case 'r':
				{
					out_file<<"Random force applied symmetrically to edge grains, maximum value: "<<sourceforceFL[c].amplitude<<" (kN), minimum value: "
						<<sourceforceFL[c].frequency<<" (kN), frequency: "<<sourceforceFL[c].phase<<" (rad/micro-s), window: "
						<<sourceforceFL[c].duration<<" musecs."<<std::endl;
					break;
				}
			}
		}
	}


	if(asymLeft == 1)
	{
		if(konstForce > 0.0)
		{
			out_file<<"Constant force applied to left edge only "<<konstForce<<" (kN), window: "<<konstLength<<" (micro-s)"<<std::endl;
		}
		for(int c=0;c<sourceforcelength;++c)
		{
			switch(sourceforce[c].type)
			{
				case 's':
				{	
					out_file<<"Sine wave applied to left edge only, amplitude: "<<sourceforce[c].amplitude<<" (kN), angular frequency: "
						<<sourceforce[c].frequency<<" (rad/micro-s)"<<"phase: "<<sourceforce[c].phase<<" (rad), window: "
						<<sourceforce[c].duration<<" musecs."<<std::endl;
					break;
				}
				case 'c':
				{
					out_file<<"Cosine wave applied to left edge only, amplitude: "<<sourceforce[c].amplitude<<" (kN), angular frequency: "
						<<sourceforce[c].frequency<<" (rad/micro-s)"<<"phase: "<<sourceforce[c].phase<<" (rad), window: "
						<<sourceforce[c].duration<<" musecs."<<std::endl;
					break;
				}
				case 't':
				{
					out_file<<"Triangle wave applied to left edge only, amplitude: "<<sourceforce[c].amplitude<<" (kN), angular frequency: "
						<<sourceforce[c].frequency<<" (rad/micro-s)"<<"phase: "<<sourceforce[c].phase<<" (rad), window: "
						<<sourceforce[c].duration<<" musecs."<<std::endl;
					break;
				}
				case 'w':
				{
					out_file<<"Saw wave applied to left edge only, amplitude: "<<sourceforce[c].amplitude<<" (kN), angular frequency: "
						<<sourceforce[c].frequency<<" (rad/micro-s)"<<"phase: "<<sourceforce[c].phase<<" (rad), window: "
						<<sourceforce[c].duration<<" musecs."<<std::endl;
					break;
				}
				case 'q':
				{
					out_file<<"Square wave applied to left edge only, amplitude: "<<sourceforce[c].amplitude<<" (kN), angular frequency: "
						<<sourceforce[c].frequency<<" (rad/micro-s), "<<"phase: "<<sourceforce[c].phase<<" (rad), window: "
						<<sourceforce[c].duration<<" musecs."<<std::endl;
					break;
				}
				case 'r':
				{
					out_file<<"Random force applied to left edge only, maximum value: "<<sourceforce[c].amplitude<<" (kN), minimum value: "
						<<sourceforce[c].frequency<<" (kN), "<<" frequency: "<<sourceforce[c].phase<<" (rad/micro-s), window: "
						<<sourceforce[c].duration<<" musecs."<<std::endl;
					break;
				}
			}
		}
	}


	if(asymRight == 1)
	{
		if(konstForceLast > 0.0)
		{
			out_file<<"Constant force applied to right edge only "<<konstForceLast<<" (kN), window: "<<konstLengthLast<<" (micro-s)"<<std::endl;
		}
		for(int c=0;c<sourceforcelengthLast;++c)
		{
			switch(sourceforceLast[c].type)
			{
				case 's':
				{	
					out_file<<"Sine wave applied to right edge only, amplitude: "<<sourceforceLast[c].amplitude<<" (kN), angular frequency: "
						<<sourceforceLast[c].frequency<<" (rad/micro-s)"<<"phase: "<<sourceforceLast[c].phase<<" (rad), window: "
						<<sourceforceLast[c].duration<<" musecs."<<std::endl;
					break;
				}
				case 'c':
				{
					out_file<<"Cosine wave applied to right edge only, amplitude: "<<sourceforceLast[c].amplitude<<" (kN), angular frequency: "
						<<sourceforceLast[c].frequency<<" (rad/micro-s)"<<"phase: "<<sourceforceLast[c].phase<<" (rad), window: "
						<<sourceforceLast[c].duration<<" musecs."<<std::endl;
					break;
				}
				case 't':
				{
					out_file<<"Triangle wave applied to right edge only, amplitude: "<<sourceforceLast[c].amplitude<<" (kN), angular frequency: "
						<<sourceforceLast[c].frequency<<" (rad/micro-s)"<<"phase: "<<sourceforceLast[c].phase<<" (rad), window: "
						<<sourceforceLast[c].duration<<" musecs."<<std::endl;
					break;
				}
				case 'w':
				{
					out_file<<"Saw wave applied to right edge only, amplitude: "<<sourceforceLast[c].amplitude<<" (kN), angular frequency: "
						<<sourceforceLast[c].frequency<<" (rad/micro-s)"<<"phase: "<<sourceforceLast[c].phase<<" (rad), window: "
						<<sourceforceLast[c].duration<<" musecs."<<std::endl;
					break;
				}
				case 'q':
				{
					out_file<<"Square wave applied to right edge only, amplitude: "<<sourceforceLast[c].amplitude<<" (kN), angular frequency: "
						<<sourceforceLast[c].frequency<<" (rad/micro-s), "<<"phase: "<<sourceforceLast[c].phase<<" (rad), window: "
						<<sourceforceLast[c].duration<<" musecs."<<std::endl;
					break;
				}
				case 'r':
				{
					out_file<<"Random force applied to right edge only, maximum value: "<<sourceforceLast[c].amplitude<<" (kN), minimum value: "
						<<sourceforceLast[c].frequency<<" (kN), "<<" frequency: "<<sourceforceLast[c].phase<<" (rad/micro-s), window: "
						<<sourceforceLast[c].duration<<" musecs."<<std::endl;
					break;
				}
			}
		}
	}


	out_file<<"Files use these conventions:"<<std::endl;
	out_file<<"Time in mu-sec"<<std::endl;
	out_file<<"Force in kN"<<std::endl;
	out_file<<"Energy in J"<<std::endl;
	out_file<<"Distances in mm"<<std::endl;
	out_file<<"Velocities in mm/musec"<<std::endl;
	out_file<<"Density in mg/mm^3"<<std::endl<<std::endl;


	out_file.close();
	delete buffer;
}

void ChainRelPositions()
{
 
	double force = 0;
	double length = 0;
	
	if(loadingForce>0)	//if(InitialOverlapLargeandWall!=0) old line, replaced 10/13/08
	{	
		//Had to change the edge grain- wall parameters. Recall that there is one more interface than number of grains. Need to also set this for 
		//the grains next to the impurity. Updated June 15/15 MP-- automated. Wall/edge grain parameter being set to homogeneous chain.
		smalla[0] = (1 / (xn * D2)) * (sqrt(Chain[0].radius));
		smalla[nptles] =  (1 / (xn * D2)) * (sqrt(Chain[nptles-1].radius));//in compresschain i set these to 0 if there are no walls

		for (int i = 1; i < nptles; i++)
			smalla[i] = (1 / (xn * D2)) * (sqrt((Chain[i].radius*Chain[i-1].radius)/(Chain[i].radius+Chain[i-1].radius)));
 
		 //Here we will reset the appropriate D values (host-impurity is D1, impurity-impurity is D3):

		if(NumImp>0)
		{
			smalla[ImpGr] = (1 / (xn * D1)) * (sqrt((Chain[ImpGr].radius*Chain[ImpGr-1].radius)/(Chain[ImpGr].radius+Chain[ImpGr-1].radius))); //both set to host-impurity e.g. steel-rubber interaction
			smalla[ImpGr+NumImp] = (1 / (xn * D1)) * (sqrt((Chain[ImpGr+NumImp].radius*Chain[ImpGr+NumImp-1].radius)/(Chain[ImpGr+NumImp].radius+Chain[ImpGr+NumImp-1].radius))); 
		}

		if(NumImp>1)
		{
			for (int i4 = 1; i4 < NumImp; i4++) //Now reset all impurity-impurity (D3) interactions.  
			{ 
				smalla[ImpGr+i4] = (1 / (xn * D3)) * (sqrt((Chain[ImpGr+i4].radius*Chain[ImpGr+i4-1].radius)/(Chain[ImpGr+i4].radius+Chain[ImpGr+i4-1].radius))); 
			}
		}

		CompressChain(force, length);

	}
	else
	{	  
		for(int lcv = 0; lcv<nptles; ++lcv)
		{ 
			length = length + 2 * Chain[lcv].radius;
			Chain[lcv].absolutePosition = length - Chain[lcv].radius;
		}
		ChainLength = length;

		if(rightWall)
			smalla[0] = (1.0 / (xn * D2)) * (sqrt(Chain[0].radius));
		else
			smalla[0] = 0;
		if(leftWall)
			smalla[nptles] =  (1.0 / (xn* D2)) * (sqrt(Chain[nptles-1].radius));
		else
			smalla[nptles] = 0;
	
		for (int i = 1; i < nptles; i++)
			smalla[i] = (1.0 / (xn * D2)) * (sqrt((Chain[i].radius*Chain[i-1].radius)/(Chain[i].radius+Chain[i-1].radius)));
 
		//We will reset everything here too since it appears this is done twice!!! Homogeneous walls (described by D2). Host-impurity = D1, impurity-impurity = D3.
		if(NumImp>0)
		{	
			smalla[ImpGr] = (1 / (xn * D1)) * (sqrt((Chain[ImpGr].radius*Chain[ImpGr-1].radius)/(Chain[ImpGr].radius+Chain[ImpGr-1].radius))); //both set to host-impurity e.g. steel-rubber interaction
			smalla[ImpGr+NumImp] = (1 / (xn * D1)) * (sqrt((Chain[ImpGr+NumImp].radius*Chain[ImpGr+NumImp-1].radius)/(Chain[ImpGr+NumImp].radius+Chain[ImpGr+NumImp-1].radius))); 
		}

		if(NumImp>1)
		{
			for (int i4 = 1; i4 < NumImp; i4++) //Now reset all impurity-impurity (D3) interactions.  
			{	
				smalla[ImpGr+i4] = (1 / (xn * D3)) * (sqrt((Chain[ImpGr+i4].radius*Chain[ImpGr+i4-1].radius)/(Chain[ImpGr+i4].radius+Chain[ImpGr+i4-1].radius))); 
			}
		}

	}	

	ForceInChain = force;
	std::cout<<"Force is: "<<force<<"kN"<<std::endl;
	std::cout<<"Length is: "<<length<<"mm"<<std::endl;

}

void CompressChain(double &force, double &length)
{
	std::cout<<"loading chain..."<<std::endl;
	InitialOverlapLargeandWall = pow((loadingForce/smalla[nptles]*1/xn),(1/(xn-1)));

	deltas[nptles]=InitialOverlapLargeandWall;//set up initial overlap for beginning of recursion/loop

	force = xn * smalla[nptles] *pow(deltas[nptles],(xn-1.0));
 
	for(int lcv = nptles-1; lcv>=0;--lcv)
	{
		deltas[lcv] = pow((force/smalla[lcv]*1/xn),(1/(xn-1)));
	}

	if(!leftWall)
	{	
		deltas[nptles]=0;
		smalla[nptles]=0;
	} //added these ifs in b/c if there is no wall but the chain is precompressed

	if(!rightWall)			//these overlaps must be 0, ACS 10/13/08
	{
		deltas[0]=0;
		smalla[0]=0;
	}

	for(int lcv2 = 0; lcv2<nptles; ++lcv2)
	{ 
		length = length + 2*Chain[lcv2].radius - deltas[lcv2];
		Chain[lcv2].absolutePosition = length - Chain[lcv2].radius;
	}
	
	length = length - deltas[nptles];

	ChainLength = length;

	for(int lcv3=0; lcv3<=nptles; ++lcv3)
		overbefore[lcv3] = deltas[lcv3];

}


void velocityVerletStep()
{
	for (int j = 0; j < nptles; j++) 
	{
	
		Chain[j].relativeLocation += Chain[j].currentVelocity * dt + 0.5 * Chain[j].currentAccel * dt*dt;
		Chain[j].currentVelocity += 0.5 * Chain[j].currentAccel * dt;
	}

	computeAccelerations();  

	//Now we have to control the program flow and apply forces at edges appropriately.

	if(symLeftRight == 1)
	{ 
		if(konstForceSym > 0.0)    //If there is a nonzero constant force being applied symmetrically to the chain edges
		{
			window = otime - konstLengthSym;
			if(window > 0)  //means we are no longer in the window of applying the force symmetrically to the chain edges
			{
				konstForceSym = 0.0;
			}
			else   //means we are applying a constant force since we are inside the window.
			{
				Chain[0].currentAccel+=addForceDueToWaveFirstLast(otime)/Chain[0].mass;
				Chain[nptles-1].currentAccel-=addForceDueToWaveFirstLast(otime)/Chain[nptles-1].mass;
			}
		}
		else      //there may be another type of driving to apply, which is not constant force over a window. We also might modify this so we can specify driving time 
		{         //for a non-constant force. Will be done as above-- trivial. MP 06/20/15
			Chain[0].currentAccel+=addForceDueToWaveFirstLast(otime)/Chain[0].mass;
			Chain[nptles-1].currentAccel-=addForceDueToWaveFirstLast(otime)/Chain[nptles-1].mass;
		}
	}
	else
	{
		if(asymLeft == 1)
		{ 
			if(konstForce > 0.0)    //If there is a nonzero constant force being applied asymmetrically to the left edge of the chain.
			{
				windowleft = otime - konstLength;
				if(windowleft > 0)  //means we are no longer in the window of applying the force 
				{
					konstForce = 0.0;
				}
				else   //means we are applying a constant force to the left edge since we are inside the window.
				{
					Chain[nptles-1].currentAccel-=addForceDueToWave(otime);
				}
			}
			else      //there may be another type of driving to apply, which is not constant force over a window. We also might modify this so we can specify driving time 
			{         //for a non-constant force. Will be done as above-- trivial. MP 06/20/15
				Chain[nptles-1].currentAccel-=addForceDueToWave(otime);
			}
		}

		if(asymRight == 1)    //applying an asymmetric driving force to the right edge of the chain. This can be done alone or in combination with asymmetric left.
		{ 
			if(konstForceLast > 0.0)    //If there is a nonzero constant force being applied asymmetrically to the right edge of the chain.
			{
				windowright = otime - konstLengthLast;
				if(windowright > 0)  //means we are no longer in the window of applying the force 
				{
					konstForceLast = 0.0;
				}
				else   //means we are applying a constant force to the right edge since we are inside the window.
				{
					Chain[0].currentAccel+=addForceDueToWaveLast(otime);
				}
			}
			else      //there may be another type of driving to apply, which is not constant force over a window. We also might modify this so we can specify driving time 
			{         //for a non-constant force. Will be done as above-- trivial. MP 06/20/15
				Chain[0].currentAccel+=addForceDueToWaveLast(otime);
			}
		}
	}

	for (int j2 = 0; j2 < nptles; j2++) 
		Chain[j2].currentVelocity += 0.5 * Chain[j2].currentAccel * dt;

}


void computeAccelerations()
{

	double over,overnm1,forceBetw,forceFactor,moved_way,energy_comp,energy_deco,forceSmall,forceLarge;

	for (int i = 0; i < nptles; i++) // zeroing all acc in every call
		Chain[i].currentAccel = 0.0;       

	pot = 0.0;	

	/******* potential/force between neighboring ptles *************/
	for (int i2 = 0; i2 < nptles-1; i2++)// 
	{
		if ((-1*Chain[i2].relativeLocation + Chain[i2+1].relativeLocation - deltas[i2+1])<=0) //swap del	had -delta[i2]>
		{//added in deltas 5/27                  	// only when overlap
			over = Chain[i2].relativeLocation + deltas[i2+1] - Chain[i2+1].relativeLocation;  //added in deltas 5/27 swap del
			overnm1 = pow(over, (xn - 1.0));
			pot += over * overnm1 * smalla[i2+1];
			forceBetw = smalla[i2+1] * xn * overnm1;

			if (overbefore[i2+1] < over)         	// when compressing
				forceFactor = 1.0;
			else forceFactor = epsilon;        	// when decompressing

			forceBetw *= forceFactor;
   
			Chain[i2].currentAccel -= forceBetw;           // sign(-): towards smaller x				swapped
			Chain[i2+1].currentAccel += forceBetw;         // sign(+): towards larger x		 swapped

			forceBefore[i2+1] = forceBetw;

			/****** calculate energy loss between ptles ********/
			if (overbefore[i2+1] < over) 
			{                                            	// compressing
			        moved_way = over - overbefore[i2+1];
			        energy_comp = forceBetw * moved_way;
				energy_loss[i2+1] += energy_comp; 	// whenever loading
			}
			if (overbefore[i2+1] > over) 
			{						// deco.; (don't mind '='-case)
				moved_way = overbefore[i2+1] - over;
				energy_deco = forceBetw * moved_way;
        			energy_loss[i2+1] -= energy_deco; 	// whenever unloading
			}
				 
	      /*****************************************************/

			overbefore[i2+1] = over;          // update for next timestep 
		}
		else
		{
			overbefore[i2+1] = 0.0;//deltas[i2+1];//switched from 0.0 to deltas[i2+1] 7/16/03
			//pot += deltas[i2+1]	* pow(deltas[i2+1], (xn - 1.0)) * smalla[i2+1];		 //added 7/16/03
		}          // reset when no overlap
	}

	/** pot./force between fixed wall (small, x=0) <-> small ptle **/	//THIS IS REALLY LARGE PARTICLE?
	if ((-1*Chain[0].relativeLocation + deltas[0])>=0 ) 
	{	 //swap deltas
		over = -1*Chain[0].relativeLocation + deltas[0];		//swap deltas 
		overnm1 = pow(over, (xn - 1.0));
		pot += over * overnm1 * smalla[0];
		forceSmall = smalla[0] * xn * overnm1;

		if (overbefore[0] < over)
			forceFactor = 1.0;
		else forceFactor = epsilon;

		forceSmall *= forceFactor;

		if(rightWall)
			Chain[0].currentAccel += forceSmall;  	//***********************swap -+

		forceBefore[0] = forceSmall;
		
		/****** calculate energy loss at wall (small) ********/

		if (overbefore[0] < over) 
		{  						// compressing
			moved_way = over - overbefore[0];
			energy_comp = forceSmall * moved_way;
			energy_loss[0] += energy_comp;
		}

		if (overbefore[0] > over) 
		{
			moved_way = overbefore[0] - over;
			energy_deco = forceSmall * moved_way;
			energy_loss[0] -= energy_deco;
		}

		/*****************************************************/

		overbefore[0] = over;
	}
	else 
	{  
		overbefore[0] = 0.0;//deltas[0];//0.0; previously just 0.0  7/16
		//pot += deltas[0]	* pow(deltas[0], (xn - 1.0)) * smalla[0]; //added 7/16
	}

	/*** pot./force between fixed wall (large) <-> large ptle ******/
	if ((1*Chain[nptles-1].relativeLocation + deltas[nptles] )>=0) 
	{							//Added deltas 5/27 ACS  //swap delta 
		over = 1*Chain[nptles-1].relativeLocation + deltas[nptles];	//Added 5/27 ACS	//swap delta  
		overnm1 = pow(over, (xn - 1.0));
		pot += over * overnm1 * smalla[nptles];
		forceLarge = smalla[nptles] * xn * overnm1;

		if (overbefore[nptles] < over)
			forceFactor = 1.0;
		else forceFactor = epsilon;

		forceLarge *= forceFactor;
	
		forceBefore[nptles] = forceLarge;

		if(leftWall)
			Chain[nptles-1].currentAccel -= forceLarge;	//**************swap +-
		//else
			//Chain[nptles-1].currentAccel -= 0;

		/****** calculate energy loss at wall (large) ********/
		if (overbefore[nptles] < over) 
		{							// compressing
			moved_way = over - overbefore[nptles];
			energy_comp = forceLarge * moved_way;
			energy_loss[nptles] += energy_comp; 
		}
		if (overbefore[nptles] > over) 
		{
			moved_way = overbefore[nptles] - over;
			energy_deco = forceLarge * moved_way;
			energy_loss[nptles] -= energy_deco;
		}
		/*****************************************************/

		overbefore[nptles] = over;
	}
	else 
	{
		overbefore[nptles] = 0.0;//deltas[nptles];//0.0; changed by ACS 7/16
		//pot += deltas[nptles]	* pow(deltas[nptles], (xn - 1.0)) * smalla[nptles];	//added 7/16
  	}
  
	/***** real dim of acc: division by mass **********/
	for (int i3 = 0; i3 < nptles; i3++)
		Chain[i3].currentAccel = Chain[i3].currentAccel / Chain[i3].mass;

}


void spitKEs(double ttime)
{	
	double totKE=0;
	double tV;
	if(KK)
		KEouts<<(ttime);
	if(TE)
		ETot<<(ttime);


	for(int lcv = 0; lcv<nptles;++lcv)
	{	
		tV = Chain[lcv].currentVelocity;
		totKE = totKE +  (tV*tV*0.5*Chain[lcv].mass);
	}


	if(KK)
		for(int lcv = bottomGrain; lcv<=topGrain;++lcv)
		{
			tV = Chain[nptles-lcv].currentVelocity;
			KEouts<<'\t'<<(tV*tV*0.5*Chain[nptles-lcv].mass);
		}

	if(KK)
		KEouts<<std::endl;

	
	recordTotE = pot+totKE;

	if(TE)
		ETot<<'\t'<<pot<<'\t'<<totKE<<'\t'<<(pot+totKE)<<std::endl;

}


void spitVs(double ttime)
{	
	
	double tV;
	if(VV)
		Vouts<<(ttime);

	for(int lcv = bottomGrain; lcv<=topGrain;++lcv)
	{	
		tV = -Chain[nptles-lcv].currentVelocity;
				
		if(VV)
			Vouts<<'\t'<<(tV);	
	}

	if(VV)
		Vouts<<std::endl;

}


void spitXs(double ttime)
{	
	double tX;

	if(FA)
		spitFapp(ttime,storeForce);
	if(FAL)
		spitFappLast(ttime,storeForceLast);
	//ForceFile<<(ttime);
	if(FAS)
		spitFappSym(ttime,storeForceFL);////***
	if(XX)
		XRelOuts<<(ttime);
	if(FF)
		ForceAll<<(ttime);

	if(FF&&leftWall&&bottomGrain==1)
		ForceAll<<'\t'<<forceBefore[nptles];
	else if(FF&&bottomGrain!=1)
		ForceAll<<'\t'<<forceBefore[nptles-bottomGrain+1];

	for(int lcv = bottomGrain; lcv<topGrain;++lcv)
	{	
		tX = Chain[nptles-lcv].relativeLocation;
		if(FF)				//if(bottomGrain!=1||leftWall==true)
			ForceAll<<'\t'<<forceBefore[nptles-lcv];//(smalla[lcv] * xn * pow(overbefore[lcv], (xn - 1.0))* epsilon);
						
		if(XX)
			XRelOuts<<'\t'<<-(tX);	
	}

//right wall corresponds to grain 0
	
	if(FF&&rightWall)
		ForceAll<<'\t'<<forceBefore[nptles-topGrain];//(smalla[nptles] * xn * pow(overbefore[nptles], (xn - 1.0))* epsilon);
	
	if(XX)
		XRelOuts<<'\t'<<-Chain[nptles-topGrain].relativeLocation;
	
	if(FF)
		ForceAll<<std::endl;
	
	if(XX)
		XRelOuts<<std::endl;
	
}


void spitFapp(double ttime,double force) //Force applied to left grain (mp 06/15/15)
{
	AppliedForce<<(ttime)<<'\t'<<force<<std::endl;
}


void spitFappLast(double ttime,double forceLast) //Force applied to right grain (mp 06/15/15)
{
	AppliedForceLast<<(ttime)<<'\t'<<forceLast<<std::endl;
}


void spitFappSym(double ttime,double forceSym) //Force applied symmetrically to left/right edge grains.
{
	AppliedForceSym<<(ttime)<<'\t'<<forceSym<<std::endl;
}



//This function adds an acceleration to the large particle as you determine you want to do it..
//it is given a force you provide, and simply divides out the mass of the particle...
//that could easily be changed...

double addForceDueToWave(double ttime)
{
	double Fduration;

	Force = konstForce;
	
	for(int c=0;c<sourceforcelength;++c)
	{
		Fduration = sourceforce[c].duration;

		if(ttime<=Fduration){
			//std::cout<<sourceforce[c].type<<'\t'<<sourceforce[c].amplitude<<'\t'<<sourceforce[c].frequency<<'\t'<<sourceforce[c].phase<<'\t'<<std::endl;
			switch(sourceforce[c].type)
			{
				case 's':
				{	
					Force+=sourceforce[c].amplitude*sin(sourceforce[c].frequency*ttime -sourceforce[c].phase);
					break;
				}
				case 'c':
				{
					Force+=sourceforce[c].amplitude*cos(sourceforce[c].frequency*ttime -sourceforce[c].phase);
					break;
				}
				case 't':
				{
					Force+=sourceforce[c].amplitude*(absolute(4.0*((ttime-sourceforce[c].phase)*sourceforce[c].frequency/(2.0*PI)
							-floor(sourceforce[c].frequency/(2.0*PI)*(ttime-sourceforce[c].phase)+0.5)))-1.0);
					break;
				}
				case 'w':
				{
					Force+=sourceforce[c].amplitude*2.0*((ttime-sourceforce[c].phase)*sourceforce[c].frequency/(2.0*PI)
							-floor(sourceforce[c].frequency/(2.0*PI)*(ttime-sourceforce[c].phase)+0.5));
					break;
				}
				case 'q':
				{
					Force+=sourceforce[c].amplitude*sign(sin(sourceforce[c].frequency*ttime-sourceforce[c].phase));
					break;
				}
				case 'r':
				{
					if(cursign==sign(sin(sourceforce[c].phase/2*ttime-PI))) //We want to generate a random number between -1 and +1; 
					{							//Previously this was set to rand()/(RAND_MAX + 1) which actually gives a number between -1 and 0.
						randForceVal = double(2.0*rand()/(double)(RAND_MAX+1.0)-1.0)*(sourceforce[c].amplitude-sourceforce[c].frequency)/2.0
								+(sourceforce[c].amplitude+sourceforce[c].frequency)/2.0;
						cursign=cursign*-1;
					}
					Force+=randForceVal;
					break;
				}
				default :
					std::cout<<"unknown force"<<std::endl;break;
			}	
		}
	}			
	/*if(randomForce)
	{
	Force += double(rand()%10)/10
	}*/	
		//Force = 0;//+sin(ttime*PI/100);
		//initial force 0...
	//if(ttime<358)						//while time is under 50... apply the following force
	//Force = (AMP/1000)*pow((sin(pow(sin(ttime/10000),2)*(ttime*(PI*(PER/100))))),2);
	//Force = 0.5*AMP;
	//else
	//Force = 0;	//otherwise...its 0...
	//Force = 0.05*cos(ttime*4)+0.3*sin(ttime/10)+5*sin(ttime/30);
	//Force = 4 + sin(ttime) + 2*cos(ttime*PI) + 3*sin(ttime/2);
	//This function could easily be applied cotinuously... by removing the if statement
	//it could be switched up any way you want...basically you have the freedom to define
	//a piece wise function, and they can be anything you want, ... so its kinda nice...

	//just be careful not to spit too much energy into the system, otherwise you'll
	//get some results that may not be accurate due to over compression of the spheres...
	
	storeForce= Force;
	//if((int(ttime/dt) % int(ForceSpits/dt))==0)
	/*if((count)%realFspit == 0)
	{
		spitFapp(ttime,(Force));
		//Force = double(rand()%10)/10;
	}*/

	return Force/Chain[nptles-1].mass;
}

double addForceDueToWaveLast(double ttime)   //Same as structure above, but for the last grain in the chain.
{
	double Fduration;

	ForceLast = konstForceLast;
	
	for(int cl=0;cl<sourceforcelengthLast;++cl)
	{
		Fduration = sourceforceLast[cl].duration;

		if(ttime<=Fduration){
	
			switch(sourceforceLast[cl].type)
			{
				case 's':
				{	
					ForceLast+=sourceforceLast[cl].amplitude*sin(sourceforceLast[cl].frequency*ttime -sourceforceLast[cl].phase);
					break;
				}
				case 'c':
				{
					ForceLast+=sourceforceLast[cl].amplitude*cos(sourceforceLast[cl].frequency*ttime -sourceforceLast[cl].phase);
					break;
				}
				case 't':
				{
					ForceLast+=sourceforceLast[cl].amplitude*(absolute(4.0*((ttime-sourceforceLast[cl].phase)*sourceforceLast[cl].frequency/(2.0*PI)
							-floor(sourceforceLast[cl].frequency/	(2.0*PI)*(ttime-sourceforceLast[cl].phase)+0.5)))-1.0);
					break;
				}
				case 'w':
				{
					ForceLast+=sourceforceLast[cl].amplitude*2.0*((ttime-sourceforceLast[cl].phase)*sourceforceLast[cl].frequency/(2.0*PI)
							-floor(sourceforceLast[cl].frequency/(2.0*PI)*(ttime-sourceforceLast[cl].phase)+0.5));
					break;
				}
				case 'q':
				{
					ForceLast+=sourceforceLast[cl].amplitude*sign(sin(sourceforceLast[cl].frequency*ttime-sourceforceLast[cl].phase));
					break;
				}
				case 'r':
				{
					if(cursignR==sign(sin(sourceforceLast[cl].phase/2*ttime-PI)))
					{
						randForceValLast = double(2.0*rand()/(double)(RAND_MAX+1.0)-1.0)*(sourceforceLast[cl].amplitude-sourceforceLast[cl].frequency)/2.0
									+(sourceforceLast[cl].amplitude+sourceforceLast[cl].frequency)/2.0;
						cursignR=cursignR*-1;
					}
					ForceLast+=randForceValLast;
					break;
				}
				default :
					std::cout<<"unknown force"<<std::endl;break;
			}	
		}
	}			

	storeForceLast= ForceLast;
	return ForceLast/Chain[0].mass;
}


double addForceDueToWaveFirstLast(double ttime)   //Same as structure above, but for the last grain in the chain.
{
	double Fduration;

	ForceFirstLast = konstForceSym;
	
	for(int cfl=0;cfl<sourceforcelengthFL;++cfl)
	{
		Fduration = sourceforceFL[cfl].duration;

		if(ttime<=Fduration){
			switch(sourceforceFL[cfl].type)
			{
				case 's':
				{	
					ForceFirstLast+=sourceforceFL[cfl].amplitude*sin(sourceforceFL[cfl].frequency*ttime -sourceforceFL[cfl].phase);
					break;
				}
				case 'c':
				{
					ForceFirstLast+=sourceforceFL[cfl].amplitude*cos(sourceforceFL[cfl].frequency*ttime -sourceforceFL[cfl].phase);
					break;
				}
				case 't':
				{
					ForceFirstLast+=sourceforceFL[cfl].amplitude*(absolute(4.0*((ttime-sourceforceFL[cfl].phase)*sourceforceFL[cfl].frequency/(2.0*PI)
							-floor(sourceforceFL[cfl].frequency/(2.0*PI)*(ttime-sourceforceFL[cfl].phase)+0.5)))-1.0);
					break;
				}
				case 'w':
				{
					ForceFirstLast+=sourceforceFL[cfl].amplitude*2.0*((ttime-sourceforceFL[cfl].phase)*sourceforceFL[cfl].frequency/(2.0*PI)
							-floor(sourceforceFL[cfl].frequency/(2.0*PI)*(ttime-sourceforceFL[cfl].phase)+0.5));
					break;
				}
				case 'q':
				{
					ForceFirstLast+=sourceforceFL[cfl].amplitude*sign(sin(sourceforceFL[cfl].frequency*ttime-sourceforceFL[cfl].phase));
					break;
				}
				case 'r':
				{
					if(cursignSym==sign(sin(sourceforceFL[cfl].phase/2*ttime-PI)))
					{
						randForceValFL = double(2.0*rand()/(double)(RAND_MAX+1.0)-1.0)*(sourceforceFL[cfl].amplitude-sourceforceFL[cfl].frequency)/2.0
									+(sourceforceFL[cfl].amplitude+sourceforceFL[cfl].frequency)/2.0;
						cursignSym=cursignSym*-1;
					}
						ForceFirstLast+=randForceValFL;
					break;
				}
				default :
					std::cout<<"unknown force"<<std::endl;break;
			}
		}	
	}			

	storeForceFL= ForceFirstLast;
	return ForceFirstLast;  ///Check this line
}


double absolute(double val)
{
	if(val<0)
		return val*-1;
	else
		return val;
}


double sign(double val)
{
	if(val<0)
		return -1;
	else
		return 1;
}


//Adding a new function in here called resume. We will open a "Resume.txt" file, and it will have a list of relative positions, followed by a blank line, followed by a 
//list of velocities. We read these in and place them appropriately in the data structures.----MP 11/24/14.

void Resume()
{
	std::ifstream in_file2("resume.txt");  //Open resume.txt for processing.
	char data[50];                      //We will read the data in as a set of characters, then convert the string into a numeric format below.
	double DataRead;
	int Rcount = 0;

	while(!in_file2.eof()&&!in_file2.fail()){
		Rcount = Rcount + 1; 
		for (int i = bottomGrain; i <= topGrain; i++){
			if ( i != nptles) {
				in_file2.getline(data, 50, '\t');  //If it is not the last particle, then the data is separated by a tab.
			}
			else {in_file2.getline(data, 50); }  //If it is the last particle, then it will not be tab-delimited. We will have an end of line character, which is default. 
			DataRead = strtod(data,NULL);        //Now we convert the data into a double. Then we will place it appropriately in the arrays. If Rcount = 1, then we need to
			if(Rcount == 1){                      //fill in the relative position array, and if Rcount = 2, then we shall fill in the current velocities.
				Chain[nptles-i].relativeLocation=-DataRead;
			}   //Recall, we are filling the array backwards, and we also have to negate everything.
			if(Rcount == 2){
				Chain[nptles-i].currentVelocity=-DataRead;
			}
		}
	}
	in_file2.close();
}


void loadFromFile()
{
	bool newTimeRange=false;
	bool Specified=false;
	bool newGrains=false;
	bool grainsSpecified=false;
	std::string param;
	char cc;
	char ccl;
	char ccfl;
	double parval;
	double Fdur;
	
	std::ifstream in_file("parameters.txt");

	in_file>>param;

	while(!in_file.eof()&&!in_file.fail())
	{
		for(unsigned int i =0;i<param.length();++i)
			param[i]=tolower(param[i]);

		if(param!="files:"&&param!="filename:"&&param!="addforce:"&&param!="addforcelast:"&&param!="addforcefirstlast:")
			in_file>>parval;

		if(param=="n:")
		{
			nptles=int(parval);
			newGrains=true;
		}
		if(param=="dt:")
		{
			dt=parval;
			newTimeRange = true;
		}
		//if(param=="firstimpure:")       //Added June 15 2015 MP-- grain number of the first impurity--- Changed so the program has the capability to determine the position on it own..
		//{
		//	ImpGR=int(parval-1);
		//}
		if(param=="numimpure:")         //Added June 15 2015 MP-- number of mass impurities 
		{
			NumImp=int(parval);
		}
		if(param=="rhoimp:")         //Added June 15 2015 MP-- density of mass impurities 
		{
			rhoImpure=parval;
		}
		if(param=="rhofac:")         //Added June 15 2015 MP-- radius factor of mass impurities (R_imp = rhoFac*rlarge)
		{
			ImpRF=parval;
		}
		if(param=="nsteps:")
		{
			nsteps= (unsigned long long int) (parval); // ( ) added by YT. 06202009 
			newTimeRange = true;
		}
		if(param=="restart:")
		{ 
			restart = (unsigned int) (parval);  //If restart = 1, then we will go into a special mode to restart from where we left off.
		}
		if(param=="q:")
			q=parval;
		if(param=="w:")
			epsilon=1-parval;//std::cout<<"wfound"<<parval<<std::endl;
		if(param=="rho:")
			rho=parval;
		if(param=="precision:")
			DEFAULTPRECISION=int(parval);
		if(param =="exponential:")
			xn = parval;
		if(param =="preload:")
			loadingForce = parval;//kN
		

		if(param=="d:")
			D2=parval;     //D value for host chain grain-grain interaction (and currently grain-wall interaction-- this can be changed).
		if(param=="dhi:")  //D value for host chain and impurity grain interaction
			D1=parval;
		if(param=="dii:")  //D value for impurity-impurity interaction
			D3=parval;
		if(param=="wall:")
		{
			if(int(parval)==0)
			{
				leftWall=rightWall=false;
			}
			if(int(parval)==11)
			{
				leftWall=rightWall=true;
			}
			if(int(parval)==10)
			{
				leftWall=true;
				rightWall = false;
			}
			if(int(parval)==1)
			{
				leftWall=false;
				rightWall = true;
			}
		//	std::cout<<"Wallfound"<<parval<<std::endl;
		}
		
		if(param=="rlarge:")
			rlarge = parval;
		if(param=="timespits:")
			timeSpits = parval;
		if(param=="smallinitv:")
			SmallInitialVelocity = -1*parval;
		if(param=="largeinitv:")
			LargeInitialVelocity = -1*parval;
		
		if(param=="deltav:")
		{
			deltafunc* temp;
			temp=(deltafunc*)malloc((deltafunclength+1)*sizeof(deltafunc));
			 
			for(int c=0;c<deltafunclength;++c)
			{ 
				temp[c].addedVelocity = deltaVel[c].addedVelocity;
				temp[c].grain = deltaVel[c].grain;
				temp[c].time = deltaVel[c].time;
				temp[c].howmany = deltaVel[c].howmany;     //MP-- added in so we can specify a window for delta velocity driving
				temp[c].interval = deltaVel[c].interval;
			}
		
			temp[deltafunclength].grain = int(parval);
			in_file>>parval;
			temp[deltafunclength].time = (parval);
			in_file>>parval;
			temp[deltafunclength].addedVelocity = -(parval);
			in_file>>parval;
                        temp[deltafunclength].howmany = int(parval);
			in_file>>parval;
                        temp[deltafunclength].interval = int(parval);

			deltafunclength++;

			deltaVel = temp;
			temp = NULL;
			deltaVelcopy = deltaVel;   //MP-- make a copy now since later we modify deltaVel-- we need deltaVelcopy and deltafunclengthcopy for the readme file.
                        deltafunclengthcopy = deltafunclength;
		}

		if(param=="timemin:")
		{	if(restart == 1)
			{ TimeMin = parval;}
			else {TimeMin = 0.0;}    //MP-- At this point, if we specify restart as == 1, then the minimum time will be set to what we choose. Otherwise, we are not
		}                                //restarting from a previous simulation, and the minimum time will be set to 0. Note we are using TimeMin as the variable to control 
		                                 //the restarting of the program. If TimeMin is zero, we will start as per usual, but if TimeMin is not zero, we need to put the
		                                 //program in the final state it was in, and go from here.
		if(param=="timemax:")
		{	TimeMax = parval;
			Specified = true;
		}

		if(param == "files:")
		{
			FF=false;
			TE=false;
			KK=false;
			XX=false;
			VV=false;
			FA=false;
			FAL=false;
			FAS=false;
		
			
			if(!in_file.eof()&&!in_file.fail())
				in_file>>param;
			
			for(unsigned int i =0;i<param.length();++i)
			{
				param[i]=tolower(param[i]);
					
				if(param[i]=='f')
					FF=true;
				if(param[i]=='t')
					TE=true;
				if(param[i]=='k')
					KK=true;
				if(param[i]=='x')
					XX=true;
				if(param[i]=='v')
					VV=true;
				if(param[i]=='a')
					FA=true;
				if(param[i]=='l')
					FAL=true;
				if(param[i]=='s')
					FAS=true;
			}
		}

		if(param == "addforce:")    //MP-- This is for perturbing the first grain in the chain (i.e. the left-most grain) with a specified force.
		{
			asymLeft = 1;           //MP-- Set the flag to 1 (meaning true) for asymmetric driving on the left end of the chain.
			cc ='d';//default
			if(!in_file.eof()&&!in_file.fail())
				in_file>>cc;
			
			//std::cout<<cc<<std::endl;

			cc=tolower(cc);
			if(cc=='k')
			{
				if(!in_file.eof()&&!in_file.fail())
				{
					parval =0;
					in_file>>parval;
					konstForce +=parval;
					in_file>>parval;
					konstLength = parval;   //MP-- Now we specify the length in terms of micro-s. We will shut off the constant force after this interval.
				}
			}
			else
			{
				drivingforces* temp;
				temp=(drivingforces*)malloc((sourceforcelength+1)*sizeof(drivingforces));
			  
				for(int c=0;c<sourceforcelength;++c)
				{ 
					temp[c].type = sourceforce[c].type;
					temp[c].amplitude = sourceforce[c].amplitude;
					temp[c].frequency = sourceforce[c].frequency ;
					temp[c].phase = sourceforce[c].phase;
					temp[c].duration = sourceforce[c].duration;
				}
			
				temp[sourceforcelength].type = cc;
				in_file>>parval;
				temp[sourceforcelength].amplitude = (parval);
				in_file>>parval;
				temp[sourceforcelength].frequency = (parval);
				in_file>>parval;
				temp[sourceforcelength].phase = (parval);
				in_file>>parval;
				temp[sourceforcelength].duration = (parval);
				sourceforcelength++;

				sourceforce = temp;
				temp = NULL;
			}
		}

		if(param == "addforcelast:")   // This is for perturbing the last grain in the chain with a specified force (i.e. the right-most grain)
		{
			asymRight = 1;
			ccl ='d';//default
			if(!in_file.eof()&&!in_file.fail())
				in_file>>ccl;
			
			//std::cout<<ccl<<std::endl;

			ccl=tolower(ccl);
			if(ccl=='k')
			{
				if(!in_file.eof()&&!in_file.fail())
				{
					parval =0;
					in_file>>parval;
					konstForceLast +=parval;
					in_file>>parval;
					konstLengthLast = parval;  //Same as above for konstLength.
				}
			}
			else
			{
				drivingforces* tempLast;
     				tempLast=(drivingforces*)malloc((sourceforcelengthLast+1)*sizeof(drivingforces));
			  
				for(int cl=0;cl<sourceforcelengthLast;++cl)
				{
					tempLast[cl].type = sourceforceLast[cl].type;
					tempLast[cl].amplitude = sourceforceLast[cl].amplitude;
					tempLast[cl].frequency = sourceforceLast[cl].frequency ;
					tempLast[cl].phase = sourceforceLast[cl].phase;
					tempLast[cl].duration = sourceforceLast[cl].duration;
				}
			
				tempLast[sourceforcelengthLast].type = ccl;
				in_file>>parval;
				tempLast[sourceforcelengthLast].amplitude = (parval);
				in_file>>parval;
				tempLast[sourceforcelengthLast].frequency = (parval);
				in_file>>parval;
				tempLast[sourceforcelengthLast].phase = (parval);
				in_file>>parval;
				tempLast[sourceforcelengthLast].duration = (parval);
				sourceforcelengthLast++;

				sourceforceLast = tempLast;
				tempLast = NULL;
			}
		}


		if(param == "addforcefirstlast:")   //This is for both ends of the chain to be symmetrically perturbed by a specified force.
		{
			symLeftRight = 1;
			ccfl ='d';//default
			if(!in_file.eof()&&!in_file.fail())
				in_file>>ccfl;
			
			//std::cout<<cc<<std::endl;

			ccfl=tolower(ccfl);
			if(ccfl=='k')
			{
				if(!in_file.eof()&&!in_file.fail())
				{
					parval =0;
					in_file>>parval;
					konstForceSym +=parval;
					in_file>>parval;
					konstLengthSym = parval;   //length of interval (in micro-s) of driving the chain symmetrically with constant force.
				}
			}
			else
			{
				drivingforces* tempFL;
				tempFL=(drivingforces*)malloc((sourceforcelengthFL+1)*sizeof(drivingforces));
			  
				for(int cfl=0;cfl<sourceforcelengthFL;++cfl)
				{
					tempFL[cfl].type = sourceforceFL[cfl].type;
					tempFL[cfl].amplitude = sourceforceFL[cfl].amplitude;
					tempFL[cfl].frequency = sourceforceFL[cfl].frequency ;
					tempFL[cfl].phase = sourceforceFL[cfl].phase;
					tempFL[cfl].duration = sourceforceFL[cfl].duration;
				}
			
				tempFL[sourceforcelengthFL].type = ccfl;
				in_file>>parval;
				tempFL[sourceforcelengthFL].amplitude = (parval);
				in_file>>parval;
				tempFL[sourceforcelengthFL].frequency = (parval);
				in_file>>parval;
				tempFL[sourceforcelengthFL].phase = (parval);
				in_file>>parval;
				tempFL[sourceforcelengthFL].duration = (parval);
				sourceforcelengthFL++;

				sourceforceFL = tempFL;
				tempFL = NULL;
			}
		}


		if(param == "grains:")
		{
			bottomGrain= (int(parval));
			in_file>>parval;
			topGrain= (int(parval));
			grainsSpecified=true;
		}

		if(param == "filename:")
		{
			in_file>>param;
			fileName = "_" + param;
		}

		in_file>>param;
	}

	in_file.close();


  ////MP----This is where I am modifying the //deltaVel structure a bit: Need to run a couple of loops:

        if(deltafunclength != 0)
        {
        	int dvmax = deltafunclength;
        	std::cout << dvmax <<std::endl; 
        	int dcount = 0;
        	int c2count = 0;
		int c3count = 0;
        	int dvmax2 = 0;
        	int dindex = 0;
                int ccount = 0;
        	deltafunclength = 0;
        	for(int c2=0;c2<dvmax;++c2)
        	{   
			c3count = deltaVel[c2].howmany;
			if(c3count == -1){
				c3count = 1;
			}
        		c2count = c2count + c3count;
        	}

        	deltafunc* temp;
        	temp=(deltafunc*)malloc((c2count)*sizeof(deltafunc));
           
        	for(int c=0;c<dvmax;++c){
        		dvmax2 = deltaVel[c].howmany;
			if(dvmax2 != -1){	
           			if(c != 0){
					ccount = deltaVel[c-1].howmany;
					if(ccount == -1){ 
						ccount = 1;   	   //Can't add deltaVel[c-1].howmany since this quantity is -1 to specify a window. 
					}
                 			dcount = dcount + ccount;
               			}
               			for(int d=0;d<dvmax2;++d){   
          		        	dindex = d + dcount;
          		        	temp[dindex].addedVelocity = deltaVel[c].addedVelocity;
	  		        	temp[dindex].grain = deltaVel[c].grain;
  					temp[dindex].time = deltaVel[c].time+ d*deltaVel[c].interval;        //Have to adjust the time--- this is where running the loop becomes useful!    
	           			temp[dindex].howmany = 1;  //Added this line in oct.11/14 MP.  Now if the d loop runs 1 time, temp will have dvmax2 entries. so in the next run of the c loop we need
                   			temp[dindex].interval = 0;
					deltafunclength++;         //to take this into account. This is why we have this additional variable dcount. I wanted to generalize it to more than just 2 edge 
				}				   //grains hence the loop. 
     			}
			else{
				if(c != 0){
 					ccount = deltaVel[c-1].howmany;
					if(ccount == -1){ 
						ccount = 1;   	   //Can't add deltaVel[c-1].howmany since this quantity is -1 to specify a window. 
					}
                 			dcount = dcount + ccount;
        	       		}
        	       		dindex = dcount;
        	  		temp[dindex].addedVelocity = deltaVel[c].addedVelocity;
		  		temp[dindex].grain = deltaVel[c].grain;
  				temp[dindex].time = deltaVel[c].time;        //Have to adjust the time--- this is where running the loop becomes useful!    
		           	temp[dindex].howmany = -1;
				temp[dindex].interval = deltaVel[c].interval;   
        	           	deltafunclength++;            
        	       	}
		}		
			
		deltaVel = temp;      //Rewrite over the existing deltaVel structure, replacing it with the appropriately restructured temp, and then as in the loops above, reset temp to nothing.
		temp = NULL;

                for(int i=0; i<deltafunclength; i++){	
			std::cout<<deltaVel[i].grain<<'\t'<<deltaVel[i].time<<'\t'<<deltaVel[i].addedVelocity<<'\t'<<deltaVel[i].howmany<<'\t'<<deltaVel[i].interval<<std::endl;
		}
	}


	if(newTimeRange&&!Specified)
		TimeMax = TimeMin + nsteps*dt;
	if(newGrains&&!grainsSpecified)
		topGrain = nptles;


	if(asymLeft==1){		//MP--added June 2015-- if we want to apply a non-constant force for the entire duration of the simulation, we can specify the duration as -1.
		for(int c1=0;c1<sourceforcelength;++c1)
		{
			Fdur = sourceforce[c1].duration;
			if(Fdur==-1){
				sourceforce[c1].duration = TimeMax;
			}
		}
	}

	if(asymRight==1){
		for(int c2=0;c2<sourceforcelengthLast;++c2)
		{
			Fdur = sourceforceLast[c2].duration;
			if(Fdur==-1){
				sourceforceLast[c2].duration = TimeMax;
			}
		}
	}

	if(symLeftRight==1){
		for(int c3=0;c3<sourceforcelengthFL;++c3)
		{
			Fdur = sourceforceFL[c3].duration;
			if(Fdur==-1){
				sourceforceFL[c3].duration = TimeMax;
			}
		}
	}





}


void makeResumeFile()
{
	double tXr,tVr;
	std::string ResumeName;
	std::ostringstream *buffer;
	std::ofstream out_file2;

	buffer = new std::ostringstream();
	
    (*buffer) << "resume.txt";
    ResumeName = buffer->str();

	out_file2.open(ResumeName.c_str()); 
 	out_file2.precision(DEFAULTPRECISION);

	//Output relative locations followed by blank line, followed by grain velocities. This file will be used if one decides to restart the program later.	
	for(int lcv = bottomGrain; lcv<topGrain;++lcv)
	{	
		tXr = Chain[nptles-lcv].relativeLocation;
		out_file2<<-(tXr)<<'\t';	
	}
	out_file2<<std::endl<<std::endl;

	for(int lcv = bottomGrain; lcv<=topGrain;++lcv)
	{	
		tVr = -Chain[nptles-lcv].currentVelocity;
		out_file2<<(tVr)<<'\t';	
	}
	out_file2<<std::endl;

	out_file2.close();
	delete buffer;
}
