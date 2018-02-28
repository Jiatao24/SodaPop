// -*- c++ -*-
#include "global.h"
#include "rng.h"

/*SodaPop
Copyright (C) 2017 Louis Gauthier

    SodaPop is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    SodaPop is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with SodaPop.  If not, see <http://www.gnu.org/licenses/>.
 */

class Gene 
{
    private:
        int g_num_;		//numeric ID pointing to primordial gene
        int ln_;		//length nuc seq
        int la_;		//length aa seq

        int Na_;		//number of non-synonymous substitutions
        int Ns_;		//number of sysnonymous substitutions
        
        std::string nucseq_;	//nucleotide sequence
        
        double dg_;		//stability
        double f_;      //gene "fitness"

        double stochastic_shape;
        double stochastic_scale;

    public:
        double conc;	//concentration
        int stochastic_conc;
        double e;		//essentiality: 1-if directly involved in replication, 0-otherwise
   
    public:
        static double shape_;
        static double scale_;
        static double mean_;
        static double stdev_;
        static std::gamma_distribution<double> gamma_;
        static std::normal_distribution<> normal_;

        Gene();
        Gene(int, std::string, double);
        Gene(std::fstream&);
        Gene(const Gene&); //copy constructor 
        ~Gene(); 
      
        bool operator==(Gene&);
        Gene& operator=(const Gene&);

        static void initGamma();
        static double RandomGamma();
        static void initNormal();
        static double RandomNormal();
      
        double Mutate_Stabil_Gaussian(int, int);
        std::string Mutate_Stabil(int, int);
        double Mutate_Select_Dist(int, int);
        std::string Mutate_Select(int, int);

        void Update_Sequences(std::string);
     
        void ch_dg(const double a){ dg_ = a; }
        void ch_f(const double a){ f_ = a; }
        void ch_Na(const int a){ Na_ = a; }
        void ch_Ns(const int a){ Ns_ = a; }
        void ch_ln(int l){ln_ = l;}
        void ch_la(int l){la_ = l;}
        void ch_gnum(int i){g_num_ = i;}

        const int num(){return g_num_;}
        const int length(){return ln_;}
        const int AAlength(){return la_;}
        const std::string nseq(){return nucseq_;}
        const double dg(){return dg_;}
        const double f(){return f_;}
        const int Ns(){return Ns_;}
        const int Na(){return Na_;}

        double CheckDG();
        double Pnat();
        double functional(bool stochastic=false);
        double misfolded(bool stochastic=false);
        int update_stochastic_conc();
};

double Gene::shape_ = 1.0;
double Gene::scale_ = 1.0;
double Gene::mean_ = 1.0;
double Gene::stdev_ = 1.0;

/* auto engine = ProperlySeededRandomEngine(); */
std::gamma_distribution<double> Gene::gamma_ = std::gamma_distribution<double>(Gene::shape_, Gene::scale_);
std::normal_distribution<> Gene::normal_ = std::normal_distribution<>(Gene::mean_, Gene::stdev_);

Gene::Gene(){
      g_num_ = 0;
      ln_ = 0; la_ = 0;
      Na_ = 0; Ns_ = 0;
      nucseq_ = ""; 
      dg_ = 1;
      f_ = 1;
      conc = 1;
      e = 0;
}

//Input: gene number, nuc. sequence, concentration
Gene::Gene(const int i, const std::string a, double c)
{
    if((a.length() % 3) != 0){
        std::cerr << "Invalid length for nucleotide sequence: " << a.length() << std::endl;
        exit(2);
    }
    else{
        g_num_=i;
        nucseq_=a;
        ln_=a.length();
        la_=ln_/3;
        std::string aaseq = GetProtFromNuc(nucseq_);
        //check for stop codons in midsequence
        std::string::size_type loc = aaseq.find("X", 0 );
        if(loc != std::string::npos){
            std::cerr << "ERROR: DNA sequence has STOP codons in the middle"<<std::endl;
            exit(2);
        }           
        dg_= 1;
        f_= 1;
        conc = c;
        stochastic_conc = c;
        e = 0;
        Na_ = 0;
        Ns_ = 0;
    }
}

//Input: gene file
Gene::Gene(std::fstream& gene_in)
{
    std::string line;
    // Read gene file line by line
    while (!gene_in.eof()){
        getline(gene_in, line);
        std::string word;
        std::istringstream iss(line, std::istringstream::in);
        iss >> word;
        if (word == "Gene_NUM")
        {
            iss >> word; 
            g_num_ = atoi(word.c_str());
        }
        else if (word == "N_Seq")
        { 
            iss >> nucseq_;
            ln_ = nucseq_.length();
            if ((ln_ % 3) != 0)
            {
                  std::cerr << "Invalid length for nucleotide sequence: " << ln_ << std::endl;
                  exit(2);
            }
            else
            {
                  la_ = ln_/3;
                  std::string aaseq = GetProtFromNuc(nucseq_);

                  // check stop codons in midsequence
                  std::string::size_type loc = aaseq.find("X", 0 );
                  if (loc != std::string::npos)
                  {
                        std::cerr << "ERROR: DNA sequence has STOP codons in the middle"<<std::endl;
                        exit(2);
                    }           
            }
        }
        else if (word == "E")
        {
            iss >> word;
            e = atoi(word.c_str());
        }
        else if (word == "CONC")
        {
            iss >> word;
            conc = atof(word.c_str());
            stochastic_conc = conc;
        }
        else if (word == "DG")
        { 
            iss >> word;
      	    dg_ = atof(word.c_str());
            dg_ = exp(-dg_/kT);
        }
        else if (word == "F")
        { 
            iss >> word;
            f_ = atof(word.c_str());
        }
        else if (word == "STOCHASTIC_SHAPE")
        {
            iss >> word;
            stochastic_shape = atof(word.c_str());
        }
        else if (word == "STOCHASTIC_SCALE")
        {
            iss >> word;
            stochastic_scale = atof(word.c_str());
        }
        else if (word == "//"){;}//do nothing
    }
    Na_ = 0; //default
    Ns_ = 0;
}

// copy constructor
Gene::Gene(const Gene& G)
{
    g_num_ = G.g_num_;
    ln_ = G.ln_;
    la_ = G.la_;
    nucseq_ = G.nucseq_;
    dg_ = G.dg_;
    f_ = G.f_;
    conc = G.conc;
    stochastic_shape = G.stochastic_shape;
    stochastic_scale = G.stochastic_scale;
    stochastic_conc = G.stochastic_conc;
    e = G.e;
    Na_ = G.Na_;
    Ns_ = G.Ns_;
}

Gene::~Gene()
{
}

//Genes are equal if DNA sequence and concentration are equal.
bool Gene::operator== (Gene& G) 
{
    std::string temp = G.nseq();
    if ( (temp.compare(nucseq_) == 0) && (conc == G.conc) ) return true;
    else return false;
}

// assignment overloading
Gene& Gene::operator=(const Gene& A)
{ 
    if (this != &A){
        this->g_num_ = A.g_num_;
        this->ln_ = A.ln_;
        this->la_ = A.la_;
        this->dg_ = A.dg_;
        this->f_ = A.f_;
        this->stochastic_shape = A.stochastic_shape;
        this->stochastic_scale = A.stochastic_scale;
        this->stochastic_conc = A.stochastic_conc;
        this->conc = A.conc;
        this->e = A.e;
        this->Na_ = A.Na_;
        this->Ns_ = A.Ns_;
        (this->nucseq_).assign(A.nucseq_);
    }
    return *this;
}

/*
This version of the mutation function draws the DDG value from a gaussian distribution
with a shifting mean to mimic sequence depletion.
*/
double Gene::Mutate_Stabil_Gaussian(int i, int j)
{ 
    if(i>=ln_){
        std::cerr << "ERROR: Mutation site out of bounds."<< std::endl;
        exit(2);
    }       

    double ran = randomNumber();
       
    if(ran <= fNs){//non-synonymous mutation
        double temp = Ran_Gaussian(1.0, 1.7);
        double x = exp(-temp/kT);

        this->dg_ *= x;
        this->Na_ += 1;

        return x;
    }
    else{
        this->Ns_ += 1;

        return 1;
    }
}

/*
This version of the mutation function gets the DDG value from the DDG matrix
input by the user.
INPUT: 
    i -> site to mutate
    j -> bp to mutate to
*/
std::string Gene::Mutate_Stabil(int i, int j)
{ 
    // extract codon to be mutated
    int cdn_ndx = (i%3);
    int cdn_start = i - cdn_ndx; 
    int resi = cdn_start/3;

    // fetch current codon
    std::string cdn_curr = nucseq_.substr(cdn_start, 3);
    // fetch current amino acid
    int aa_curr = GetIndexFromCodon(cdn_curr);
    std::string cdn_new = cdn_curr;

    std::string s = PrimordialAASeq.at(g_num_);     
    int aa_primo = GetIndexFromAA(s.at(resi));

    // get mutated bp
    std::string bp = AdjacentBP( cdn_curr.substr(cdn_ndx, 1), j); //new BP
   
    // mutate codon
    cdn_new.replace(cdn_ndx, 1, bp);
    // check for stop codon
    cdn_new = n3_to_n3(cdn_new, cdn_curr, cdn_ndx);
    // get new amino acid
    int aa_new = GetIndexFromCodon(cdn_new);
    
    // get DDG value from matrix
    double x = matrix[g_num_][resi][aa_new-1];

    std::string mutation = std::to_string(g_num_) + '\t' + GetProtFromNuc(cdn_curr) + '\t' + std::to_string(resi) + '\t' + GetProtFromNuc(cdn_new);

    // fetch primordial amino acid

    //Ignore mutations to and from CYSTEINE
    if( (aa_new==2) || (aa_curr==2)){
        return "CYSTEINE\tNA\tNA\tNA";
    }

    //Case unphysical DDG estimate
    if( x>DDG_min || x<DDG_max){
        return "UNPHYSICAL\tNA\tNA\tNA";
    }

    if( aa_curr == aa_new){//SILENT
          nucseq_.replace(cdn_start, 3, cdn_new);
          Ns_ += 1;
          return "SILENT\tNA\tNA\tNA";
    }
    else if(aa_primo == aa_new){//REVERT TO WT

          double x_curr = matrix[g_num_][resi][aa_curr-1];
          assert( x_curr<DDG_min || x_curr>DDG_max); 
          
          dg_ /= x_curr;
          nucseq_.replace(cdn_start, 3, cdn_new);
          Na_ += 1;
          return mutation;
    }
    else{//TYPICAL NON-SYNONYMOUS

          double x_curr = matrix[g_num_][resi][aa_curr-1];
          assert( x_curr<DDG_min || x_curr>DDG_max); 

          // assign new DG value
          // division account for wildtype background
          dg_ /= x_curr;
          dg_ *= x;
          nucseq_.replace(cdn_start, 3, cdn_new);
          Na_ += 1;
          return mutation;
    }
}

/*
This version of the mutation function draws the selection coefficient value from a gamma or normal distribution
* VZ: It seems only RandomNormal is possible here?
*/
double Gene::Mutate_Select_Dist(int i, int j)
{ 
    if(i>=ln_){
        std::cerr << "ERROR: Mutation site out of bounds."<< std::endl;
        exit(2);
    }       

    double ran = randomNumber();
       
    if(ran <= fNs){//non-synonymous mutation
        double s = RandomNormal();
        double wf = f_ + s;
        f_ = wf * f_;
        Na_ += 1;
        return s;
    }
    else{
        this->Ns_ += 1;
        return 1;
    }
}

/*
This version of the mutation function gets the DDG value from the DDG matrix
input by the user.
INPUT: 
    i -> site to mutate
    j -> bp to mutate to

VZ: I believe the function description is wrong because it was copied
    from Mutate_Stabil, because in this case, the matrix stores DMS
    data, not DDG.
*/
std::string Gene::Mutate_Select(int i, int j)
{ 
    // extract codon to be mutated
    int cdn_ndx = (i%3);
    int cdn_start = i - cdn_ndx; 
    int resi = cdn_start/3;

    // fetch current codon
    std::string cdn_curr = nucseq_.substr(cdn_start, 3);
    // fetch current amino acid
    int aa_curr = GetIndexFromCodon(cdn_curr);
    std::string cdn_new = cdn_curr;

    std::string s = PrimordialAASeq.at(g_num_);     

    // get mutated bp
    std::string bp = AdjacentBP( cdn_curr.substr(cdn_ndx, 1), j); //new BP
   
    // mutate codon
    cdn_new.replace(cdn_ndx, 1, bp);
    // check for stop codon
    cdn_new = n3_to_n3(cdn_new, cdn_curr, cdn_ndx);
    // get new amino acid
    int aa_new = GetIndexFromCodon(cdn_new);
    
    // get selection coefficient from matrix
    double new_s = matrix[g_num_][resi][aa_new-1];

    std::string mutation = std::to_string(g_num_) + '\t' + GetProtFromNuc(cdn_curr) + '\t' + std::to_string(resi) + '\t' + GetProtFromNuc(cdn_new);

    // fetch primordial amino acid

    //Ignore mutations to and from CYSTEINE
    if( (aa_new==2) || (aa_curr==2)){
        return "CYSTEINE\tNA\tNA\tNA";
    }

    if( aa_curr == aa_new){//SILENT
          nucseq_.replace(cdn_start, 3, cdn_new);
          Ns_ += 1;
          return "SILENT\tNA\tNA\tNA";
    }
    else{//TYPICAL NON-SYNONYMOUS 

          // assign new fitness value
          double new_f = f_ + new_s;
          f_ = f_ * new_f;
          nucseq_.replace(cdn_start, 3, cdn_new);
          Na_ += 1;
          return mutation;
    }
}

void Gene::initGamma()
{
    Gene::gamma_.param(std::gamma_distribution<double>::param_type(Gene::shape_, Gene::scale_));
}

void Gene::initNormal()
{
    Gene::normal_.param(std::normal_distribution<>::param_type(Gene::mean_, Gene::stdev_));
}

double Gene::RandomGamma()
{
    return Gene::gamma_(g_rng);
}

double Gene::RandomNormal()
{
    return Gene::normal_(g_rng);
}


// Updates the current DNA sequence
void Gene::Update_Sequences(const std::string DNAsequence)
{ 
    int l = DNAsequence.length();

    if(l != ln_)
    {
        std::cerr << "ERROR: Replacing DNA sequence with a non-equal length DNA. "<< std::endl;
        exit(2);
    }       

    nucseq_ = DNAsequence;
}

// from Privalov 1979 (see also: Serohijos & Shakhnovich 2013)
// Boltzmann probability of the gene product to be in the native state
double Gene::Pnat()
{
    return this->dg()/(1+this->dg());
}

// Number of functional copies in the cell
// concentration * pnat
double Gene::functional(bool stochastic)
{
    if (stochastic)
    {
        return stochastic_conc * Pnat();
    }
    else
    {
        return conc * Pnat();
    }
}

// Number of misfolded copies in the cell
double Gene::misfolded(bool stochastic)
{
    if (stochastic)
    {
        return stochastic_conc * (1 - Pnat());
    }
    else
    {
        return conc * (1 - Pnat());
    }
}

// Draw a new value for stochastic concentration
int Gene::update_stochastic_conc()
{
    stochastic_conc = (int)(0.5 +
        Gene::gamma_(g_rng, std::gamma_distribution<double>::param_type(
                         stochastic_shape, stochastic_scale)));
    return stochastic_conc;
}
