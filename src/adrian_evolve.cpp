#include "Cell.h"
//#include <mpi.h>

/*
DESCRIPTION: Implementation of microorganism pupulation dynamics and evolution
inspired from "omics" data and protein biophysics.

AUTHOR: Serohijos, Adrian W.R.
Chemistry & Chemical Biology, Harvard University
serohij@fas.harvard.edu

Version: 
*/


int main(int argc, char *argv[]){
  if(argc != 3){
    std::cerr <<"evolve.linux <start file> <job ID> \n";
    exit(1);
  }

  int node = atoi(argv[2]);
  //MPI_Init(&argc,&argv);
  //MPI_Comm_rank(MPI_COMM_WORLD, &node);
                                                                                                                                                                                
  /***********************
     START MPI ROUTINES
  ************************/

    char buffer[200];
    int seed;

    double T0=0, DTEVENTS = 1;
    int DT =1;

    //Generation counter
    int GENERATION_CTR = 1; 			//default, 1 --> 1st generation
    int GENERATION_MAX = GENERATION_CTR + 1; 	//default, run evolution for only one generation

    std::string outfile, eventslogfile, genelistfile;
    std::string snapfile, startsnapfile, pDDGfile, generationfile;

    std::vector <Cell> Cell_arr;
    Cell_arr.reserve(POPSIZEMAX);

    //Initialize DDG matrix
    for(int i = 0; i != max_gene; ++i)
      for(int j = 0; j != max_resi; ++j)
        for(int k = 0; k != 20; ++k)
          pDDG[i][j][k] = 1; //i.e., exp(-0/kT) = 1

    //Get input
    sprintf(buffer, "%s.%d", argv[1], node);
    std::fstream startf (buffer);
    if ( !startf.is_open() ) {
      std::cerr << "File could not be opened" << argv[1] << std::endl ;
      exit(1);
    }

    std::string line;
    while ( !startf.eof() )
    {
      getline(startf,line);
      std::string word;
      std::istringstream iss(line, std::istringstream::in);
      iss >> word;
      
      if ( word == "T0" ){
        iss >> word;
        T0 = atof(word.c_str());
      }
      else if ( word == "DT" ){
        iss >> word;
        DT = atoi(word.c_str());
        assert(DT > 0);
      }
      else if ( word == "GEN_0" ){
        iss >> word;
        GENERATION_CTR = atoi(word.c_str());
        assert(GENERATION_CTR > 0);
      }
      else if ( word == "GEN_MAX" ){
        iss >> word;
        GENERATION_MAX = atoi(word.c_str());
        assert(GENERATION_MAX > GENERATION_CTR);
      }
      else if ( word == "DTEVENTS" ){
        iss >> word;
        DTEVENTS = atof(word.c_str());
        assert( DT > 0);
      }
      else if ( word == "OUTFILE" ){
        iss >> outfile;  
      }
      else if ( word == "EVENTSLOG" ){
        iss >> eventslogfile;  
      }
      else if ( word == "SNAPSHOT" ){
        iss >> snapfile;  
      }
      else if ( word == "C" ){
        iss >> startsnapfile; 
      }
      else if ( word == "GENE_LIST" ){
        iss >> genelistfile; 
      }
      else if ( word == "PDDG" ){
        iss >> pDDGfile;//ddg matrix file
      }
      else if ( word == "XCHANGE" ){
        iss >> word;
        xrate = atof(word.c_str());
        iss >> word;
        xSD = atof(word.c_str());
      }
      else if ( word == "GENERATION" ){
        iss >> generationfile;
      }
      else if ( word == "SEED" ){
        iss >> word;
        seed = atoi(word.c_str());
      }
    }
    startf.close();


    LoadPrimordialGenes(genelistfile);



    //Extract Primordial DDG matrix
    std::fstream temp (pDDGfile.c_str());//convert std::string to char
    if ( !temp.is_open() ) {
      std::cerr << "File could not be open: "<< pDDGfile << std::endl;
      exit(1);
    }

    int gene_num = 0;
    while( !temp.eof() ){
      std::string word;
      getline(temp,line);
      std::istringstream iss(line, std::istringstream::in);
      iss >> word;

      if ( word == "Gene_NUM"){
        iss >> word;
        gene_num = atoi(word.c_str());
      }
      else if ( word == "DDG"){
        iss >> word;
        int i = atoi(word.c_str());       //residue index
        for(int j = 0; iss>>word; j++){
          double x = atof(word.c_str());//extract DDG values

         // assert( x>ddG_min && x<ddG_max);
          pDDG[gene_num][i-1][j] = exp(-x/kT);
        }
      }   
    }
    temp.close();
    //End pDDG

    //Open starting population snapshot
    std::fstream startsnap (startsnapfile.c_str(),std::ios::in|std::ios::binary);
    if ( !startsnap.is_open() ) {
      std::cerr << "File could not be open: "<< startsnapfile << std::endl;
      exit(1);
    }

    
    //header
    double fromsnap_frametime, fromsnap_realtime;
    int Total_Cell_Count;

    startsnap.read((char*)(&fromsnap_frametime),sizeof(double));
    startsnap.read((char*)(&fromsnap_realtime),sizeof(double));
    startsnap.read((char*)(&Total_Cell_Count),sizeof(int));
    
    int i=0;
    while( i<Total_Cell_Count && !startsnap.eof()){
      Cell A(startsnap, fromsnap_realtime);
      i++;
      Cell_arr.push_back(A);
    }


    std::cout << "i value: " << i << " Cell count: " << Total_Cell_Count << std::endl;

    if ( i!= Total_Cell_Count ){
      std::cerr << "ERROR: Wrong cell count in initial snap file. "<< std::endl;
      exit(1);
    }

    startsnap.close();
    //End reading of initial population snapshot

    std::cout << "Starting evolution ..." << std::endl;



    //Check input integrity  
    //Output the starting cell population
    Total_Cell_Count = (int)(Cell_arr.size());

    //Open EVENTSLOG file
    sprintf(buffer, "%s.%03d", eventslogfile.c_str(), node);
    std::fstream EVENTSLOG(buffer, std::ios::out);
    if ( !EVENTSLOG.is_open() ) {
      std::cerr << "Events log file could not be opened";
      exit(1);
    }
 

    EVENTSLOG.write((char*)(&Total_Cell_Count),sizeof(int));



    //Open GENERATION file
    sprintf(buffer, "%s.%03d", generationfile.c_str(), node);
    std::fstream GENERATIONLOG(buffer, std::ios::out);
    if ( !GENERATIONLOG.is_open() ) {
      std::cerr << "Generation file could not be opened";
      exit(1);
    }
    //Define generation tracking variables
    //NOTE: Save time when the number of birth counts (starting current time) 
    //	equals the current population.
    int birth_count = 0;
    int next_gen_count = Total_Cell_Count;

    //Random numbers
    srand( (unsigned)time(0) + seed);

    double brate = 1;
    double drate = 1;
    double brate_f = 1;
    double drate_f = 1;


    double mrate = 0;
    double urate = 0; //duplication rate

    //Implement evolution algortihm
    double dt=0, r1=0, r2=0;
    int event=0;
    
    int save_itr_events = 0;

    double TIME = T0;
    


    while((GENERATION_CTR < GENERATION_MAX) && (Total_Cell_Count > 0)){

      int max = Cell_arr.size();
      double rates[max];

      int j=0;
      rates[j] = Cell_arr[j].get_b();			//hop through birth events only

      for(j=1; j< max; j++){				//put all rates into matrix
         rates[j] = rates[j-1] + Cell_arr[j].get_b();// + xrate;
      }

      double SUM_rates_bxu = rates[max-1];
      //double SUM_rates_d = rates[max-1][1];


      //event occured
      dt = -(1/SUM_rates_bxu)*log(RandomNumber());
      TIME += dt;

      r1 = (SUM_rates_bxu) * RandomNumber();
      std::vector<Cell>::iterator i = Cell_arr.begin();    
      std::vector<Cell>::iterator curr;    
   
      int k;
      for(k=0; k < max  ; k++){
        if(r1 < rates[k]) break;
        i++; 
      }

      //std::cout << "Cell: " << k << std::endl;

      //What event occured - birth, gene duplication, or change in expression level?
      r2 = RandomNumber();
      event = 0;	//In this version, event is always birth. (*i).event_ID(r2);

      //double events_time = save_itr_events*DTEVENTS + T0;

      Total_Cell_Count = (int)(Cell_arr.size());

      //process event
      switch (event) {
        case 0://BIRTH

               {
   	       //Start Birth

                 //OUTPUT: Cell_ID
                 //std::cout << k;
               
                 //insert copy of the cell to array
                 curr = Cell_arr.insert(i,(*i));

        	 brate = (*curr).get_b();
      		 drate = (*curr).get_d();

                 //potentially mutate daughter (just added cell)
                 int Nmut = (*curr).ranmut_Gene(EVENTSLOG); 
                 if( Nmut>0 ){
                   brate_f = (*curr).get_b();
      		   drate_f = (*curr).get_d();
                   /*
                   //EVENTSLOG.write((char*)(&TIME),sizeof(double));
                   EVENTSLOG.write((char*)(&GENERATION_CTR),sizeof(int));
                   EVENTSLOG.write((char*)(&birth_count),sizeof(int));
                   EVENTSLOG.write((char*)(&k),sizeof(int));
                   EVENTSLOG.write((char*)(&brate),sizeof(double));
                   EVENTSLOG.write((char*)(&drate),sizeof(double));
                   EVENTSLOG.write((char*)(&brate_f),sizeof(double));
                   EVENTSLOG.write((char*)(&drate_f),sizeof(double));
                   */ 
                 }   

                 //potentially mutate original cell
                 curr++;
                 Nmut = (*curr).ranmut_Gene(EVENTSLOG); 
                 if( Nmut>0 ){
                   brate_f = (*curr).get_b();
      		   drate_f = (*curr).get_d();
                   /* 
                   //EVENTSLOG.write((char*)(&TIME),sizeof(double));
                   EVENTSLOG.write((char*)(&GENERATION_CTR),sizeof(int));
                   EVENTSLOG.write((char*)(&birth_count),sizeof(int));
                   EVENTSLOG.write((char*)(&k),sizeof(int));
                   EVENTSLOG.write((char*)(&brate),sizeof(double));
                   EVENTSLOG.write((char*)(&drate),sizeof(double));
                   EVENTSLOG.write((char*)(&brate_f),sizeof(double));
                   EVENTSLOG.write((char*)(&drate_f),sizeof(double));
                   */ 
                 }   
 
                 //update birth event counter
                 birth_count++;
                
                 //End Birth
                 //std::cout << std::endl; 
                 
                 //Start Death
                 //DEATH PROBABILITY IS PROPORTIONAL TO DEATH WEIGHT
                 max = Cell_arr.size();
                 double death_weights[max];
                
                 int j=0;
                 death_weights[j] = Cell_arr[j].get_d();
                
                 for(j=1; j<max; j++){				//put all rates into matrix
                    death_weights[j] = rates[j-1] + Cell_arr[j].get_d();
                 }

                 double SUM_rates_d = death_weights[max-1];
                
                 double r3 = (SUM_rates_d) * RandomNumber();
                 std::vector<Cell>::iterator dying_cell = Cell_arr.begin();    
                
                 int l=0;
                 for(l=0; l < max  ; l++){
                   if(r3 < death_weights[l]) break;
                   dying_cell++; 
                 }

  	         Cell_arr.erase(dying_cell);
                 
                 //End death
                
                 if( Cell_arr.size()==0 ){
                   return -1;     
                 }
               }

               break;

        case 1://Gene Duplication

              //Duplicate gene
              //(*i).duplicate_gene();    
              //if (TIME > events_time) (*i).dump(EVENTSLOG,k);
              break;

        case 2://Change in Expression Level
              
              (*i).change_exprlevel();

              std::cout << "XP Change" << std::endl;
              //if (TIME > events_time) (*i).dump(EVENTSLOG,k);
              break;

        default://error
  	    std::cerr << "Invalid event ID." << std::endl;
  	    exit(1);
      }//end switch


      Total_Cell_Count = (int)(Cell_arr.size());


      //Check if population have double
      if(birth_count >= next_gen_count){
         sprintf(buffer,"%15.5f %8d", TIME, Total_Cell_Count);
         GENERATIONLOG << buffer << std::endl;

         //update generation counter
         GENERATION_CTR++;
         
         //Save population snapshot every DT generations
         if( (GENERATION_CTR % DT) == 1){

           //save population snapshot
           sprintf(buffer,"%s.gen%010d.snap.%03d",snapfile.c_str(), GENERATION_CTR, node); 

           //Open snapshot file
           std::fstream OUT2(buffer, std::ios::out);
           if ( !OUT2.is_open() ) {
             std::cerr << "Snapshot file could not be opened";
             exit(1);
           }
    
           double frame_time = TIME;
   
           OUT2.write((char*)(&frame_time),sizeof(double));
           OUT2.write((char*)(&TIME),sizeof(double));
           OUT2.write((char*)(&Total_Cell_Count),sizeof(int));
    
           int l=1;
           for(std::vector<Cell>::iterator k = Cell_arr.begin(); k != Cell_arr.end(); ++k){
             (*k).dump(OUT2,l);
             l++;
           }
    
           OUT2.close();
         }

         //next generation occurs when ...
         next_gen_count = Total_Cell_Count; 
         birth_count = 0;
      }

    }//end while
   
    EVENTSLOG.close(); 
    GENERATIONLOG.close();
    //OUT.close();  
    
    std::cout << "Done." << std::endl;  

  return 0;
}//end main
