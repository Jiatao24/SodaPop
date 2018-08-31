#include <cmath>
#include "nlohmann/json.hpp"

#include "Gene.h"
#include "global.h"

// for convenience
using json = nlohmann::json;

std::gamma_distribution<> Gene::gamma_ = std::gamma_distribution<>(1.0, 1.0);
std::normal_distribution<> Gene::normal_ = std::normal_distribution<>(1.0, 1.0);


Gene::Gene()
{
      g_num_ = 0;
      ln_ = 0;
      la_ = 0;
      Na_ = 0;
      Ns_ = 0;
      nucseq_ = ""; 
      dg_ = 1;
      f_ = 1;
      conc_ = 1;
      e_ = 0;
}

// Input: gene file
Gene::Gene(std::fstream& gene_in)
{
    std::string line;
    // Read gene file line by line
    while (!gene_in.eof())
    {
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
                  std::cerr << "Invalid length for nucleotide sequence: "
			    << ln_ << std::endl;
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
                        std::cerr << "ERROR: DNA sequence has STOP codons."
				  <<std::endl;
                        exit(2);
                    }           
            }
        }
        else if (word == "E")
        {
            iss >> word;
            e_ = atoi(word.c_str());
        }
        else if (word == "CONC")
        {
            iss >> word;
            conc_ = atof(word.c_str());
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
            stochastic_shape_ = atof(word.c_str());
        }
        else if (word == "STOCHASTIC_SCALE")
        {
            iss >> word;
            stochastic_scale_ = atof(word.c_str());
        }
        else if (word == "OU_TAU")
        {
            iss >> word;
            OU_tau_ = atof(word.c_str());
        }
        else if (word == "OU_DIFFUSION")
        {
            iss >> word;
            OU_diffusion_ = atof(word.c_str());
        }
        else if (word == "BIOPHYSICS")
        {
            iss >> word;
            loadBiophysicalParameters(word);
        }
        else if (word == "//")
            ;                   // no action
    }

    init_stochastic_conc();
    Na_ = 0; // default
    Ns_ = 0;
}

// // copy constructor
// // this is dangerous right now, we're not copying everything.
// Gene::Gene(const Gene& G)
// {
//     g_num_ = G.g_num_;
//     ln_ = G.ln_;
//     la_ = G.la_;
//     nucseq_ = G.nucseq_;
//     dg_ = G.dg_;
//     f_ = G.f_;
//     conc_ = G.conc_;
//     stochastic_shape_ = G.stochastic_shape_;
//     stochastic_scale_ = G.stochastic_scale_;
//     stochastic_conc_ = G.stochastic_conc_;
//     OU_tau_ = G.OU_tau_;
//     OU_diffusion_ = G.OU_diffusion_;
//     e_ = G.e_;
//     Na_ = G.Na_;
//     Ns_ = G.Ns_;
// }


// Genes are equal if DNA sequence and concentration are equal.
bool Gene::operator== (Gene& G) 
{
    std::string temp = G.nseq();
    if ( (temp.compare(nucseq_) == 0) && (conc_ == G.conc_) )
	return true;
    else
	return false;
}

// // assignment overloading
// Gene& Gene::operator=(const Gene& A)
// { 
//     if (this != &A){
//         this->g_num_ = A.g_num_;
//         this->ln_ = A.ln_;
//         this->la_ = A.la_;
//         this->dg_ = A.dg_;
//         this->f_ = A.f_;
//         this->stochastic_shape_ = A.stochastic_shape_;
//         this->stochastic_scale_ = A.stochastic_scale_;
//         this->stochastic_conc_ = A.stochastic_conc_;
//         this->OU_tau_ = A.OU_tau_;
//         this->OU_diffusion_ = A.OU_diffusion_;
//         this->conc_ = A.conc_;
//         this->e_ = A.e_;
//         this->Na_ = A.Na_;
//         this->Ns_ = A.Ns_;
//         (this->nucseq_).assign(A.nucseq_);
//     }
//     return *this;
// }


// Mutate, updating biophysical parameters
// This is the version of the function that simply changes residue.
std::string Gene::mutate(int site, int bp)
{
    std::string oldGenotype = genotype_;

    // Codon to be mutated indexing.
    int codon_index = (site % 3); // Position within the codon
    int codon_start = site - codon_index; 
    int resi = codon_start / 3; // Residue index (0-based)

    // Fetch current codon.
    std::string codon_current = nucseq_.substr(codon_start, 3);
    // Fetch current amino acid.
    int aa_current = GetIndexFromCodon(codon_current);
    std::string codon_new = codon_current; // New codon.

    // Get mutated bp.
    std::string new_bp = AdjacentBP(codon_current.substr(codon_index, 1), bp);
   
    // Mutate codon
    codon_new.replace(codon_index, 1, new_bp);
    // Check for stop codon
    codon_new = n3_to_n3(codon_new, codon_current, codon_index);
    // get new amino acid
    int aa_new = GetIndexFromCodon(codon_new);
    
    std::string mutation;

    // Ignore mutations to and from CYSTEINE.
    if ( (aa_new==2) || (aa_current==2))
    {
        mutation = "CYSTEINE\tNA\tNA\tNA";
    }

    if (aa_current == aa_new)
    {
        // SILENT
        nucseq_.replace(codon_start, 3, codon_new);
        Ns_ += 1;
        mutation = "SILENT\tNA\tNA\tNA";
    }
    // We ignore this situation where we're reverting to the primordial seq.
    // else if (aa_primo == aa_new)
    // {
    //     //REVERT TO WT

    //     double x_current = matrix[g_num_][resi][aa_current-1];
    //     assert( x_current<DDG_min || x_current>DDG_max); 
          
    //     dg_ /= x_current;
    //     nucseq_.replace(codon_start, 3, codon_new);
    //     Na_ += 1;
    //     return mutation;
    // }
    else
    {
        //TYPICAL NON-SYNONYMOUS
        nucseq_.replace(codon_start, 3, codon_new);
        Na_ += 1;
        mutation = std::to_string(g_num_) + '\t' + GetProtFromNuc(codon_current)
            + '\t' + std::to_string(resi + 1) + '\t' + GetProtFromNuc(codon_new);
    }

    identifyGenotype();
    if (genotype_ != oldGenotype)
    {
        mutation = std::to_string(g_num_) + '\t' + oldGenotype
            + '\t' + std::to_string(resi + 1) + '\t' + genotype_;
        return mutation;
    }
    else
    {
        return "";
    }
}


std::string Gene::mutate(int resid, std::string resname)
{
    int codon_start = (resid - 1) * 3;
    std::string codon_old = nucseq_.substr(codon_start, 3);

    // Validity check: resname is a valid 1-letter AA
    auto search = prot_to_num::pnum.find(resname);
    if (search == prot_to_num::pnum.end())
    {
        std::cerr
            << "Gene::mutate(resid, resname) error: invalid 1-letter AA name: "
            << resname << std::endl;
        exit(2);
    }

    // Find new codon: search through codon to resname map and take
    //   first codon found that gives desired resname.
    std::string codon_new;
    bool found_codon = false;
    for (auto& it : codon_to_prot::cprot)
    {
        if (it.second == resname)
        {
            codon_new = it.first;
            found_codon = true;
            break;
        }
    }

    if (!found_codon)
    {
        // This should not happen since we already checked the
        // validity of resname.
        std::cerr << "Gene::mutate(resid, resname) error: "
                  << "codon for resname '" << resname << "' not found"
                  << std::endl;
        exit(2);
    }

    // Do the replacement
    nucseq_.replace(codon_start, 3, codon_new);
    identifyGenotype();

    std::string mutation = std::to_string(g_num_) + '\t' + GetProtFromNuc(codon_old)
            + '\t' + std::to_string(resid) + '\t' + GetProtFromNuc(codon_new);

    return mutation;
}


/*
This version of the mutation function draws the DDG value from a gaussian distribution
with a shifting mean to mimic sequence depletion.
*/
double Gene::Mutate_Stabil_Gaussian(int i, int j)
{ 
    if (i >= ln_)
    {
        std::cerr << "ERROR: Mutation site out of bounds." << std::endl;
        exit(2);
    }       

       
    if (randomNumber() <= fNs)
    {
        // non-synonymous mutation
        double temp = Ran_Gaussian(1.0, 1.7); // why use Ran_Gaussian?
        double x = exp(-temp/kT);

        dg_ *= x;
        Na_ += 1;

        return x;
    }
    else
    {
        Ns_ += 1;
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

void Gene::setGammaParams(double shape, double scale)
{
    Gene::gamma_.param(std::gamma_distribution<>::param_type(shape, scale));
}

void Gene::setNormalParams(double mean, double stddev)
{
    Gene::normal_.param(std::normal_distribution<>::param_type(mean, stddev));
}

double Gene::randomGamma()
{
    return Gene::gamma_(g_rng);
}

double Gene::randomNormal()
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
        // std::cerr << "Old: " << nucseq_ << std::endl;
        // std::cerr << "New: " << DNAsequence << std::endl;
        exit(2);
    }       

    nucseq_ = DNAsequence;
}


// Parse json file of biophysical parameters.
void Gene::loadBiophysicalParameters(std::string path)
{
    // load json...
    std::ifstream in(path);
    json data;
    in >> data;
    keyResidueNumbers_ = data["key_residue_numbers"].get<std::vector<int>>();
    unsigned int numKeyResidues = keyResidueNumbers_.size();
    drug_penetration_factor_ = data["drug_penetration"];
    for (auto& it : data["genotypes"])
    {
        std::string key;
        // Assert length of identifier
        if (it["identifier"].size() != numKeyResidues)
        {
            throw "Invalid identifier length!";
            // All identifiers should be the same length.
        }
        for (auto& character : it["identifier"])
        {
            key += character.get<std::string>();
        }
        biophysicalParam param(it["k_cat/K_m"], it["K_i"]);
        biophysicalParameters_.insert({key, param});
    }

    identifyGenotype();
}


// Determine the biophysical genotype of the current sequence.
void Gene::identifyGenotype()
{
    std::string newGenotype;
    // Determine identity of key residues.
    for (auto& resid : keyResidueNumbers_)
    {
        int residueIndex = resid - 1;
        int codon_start = residueIndex * 3;
        std::string codon = nucseq_.substr(codon_start, 3);
        newGenotype += GetProtFromNuc(codon);
    }
    genotype_ = newGenotype;

    auto it = biophysicalParameters_.find(genotype_);
    if (it == biophysicalParameters_.end())
    {
        // Genotypes not enumerated have 0 enzymatic output (and 0 fitness).
        rel_kcat_over_km_ = 0;
        k_i_ = 1;
    }
    else
    {
        biophysicalParam param = it->second;
        rel_kcat_over_km_ = param.first;
        k_i_ = param.second;
    }
}


// from Privalov 1979 (see also: Serohijos & Shakhnovich 2013)
// Boltzmann probability of the gene product to be in the native state
double Gene::Pnat()
{
    return dg_ / (1 + dg_);
}

// Number of functional copies in the cell
// concentration * pnat
double Gene::functional(bool stochastic)
{
    if (stochastic)
    {
        return stochastic_conc_ * Pnat();
    }
    else
    {
        return conc_ * Pnat();
    }
}

// Number of misfolded copies in the cell
// concentration * (1 - pnat)
double Gene::misfolded(bool stochastic)
{
    if (stochastic)
    {
        return stochastic_conc_ * (1 - Pnat());
    }
    else
    {
        return conc_ * (1 - Pnat());
    }
}

// Returns enzymatic output, which depends on abundance and
// biophysical catalytic properties.
double Gene::enzymaticFlux(bool stochastic)
{
    double abundance;
    if (stochastic)
        abundance = stochastic_conc_ * Pnat();
    else
        abundance = conc_ * Pnat();
    // note that conc_ and stochastic_conc_ are actually abundances
    // Multiply by 1e9 to convert to nM.
    double concentration = 1e9 * abundance / AVOGADRO / CELL_VOLUME;
    double flux = (rel_kcat_over_km_ * concentration) /
        (1 + drug_penetration_factor_ * DRUG_CONCENTRATION / k_i_);

    return flux;
}


double Gene::init_stochastic_conc()
{
    // If OU_variance is non-zero, then stochastic_concentration_
    //   is drawn from a normal distribution. 
    double OU_variance = OU_tau_ * OU_diffusion_ / 2.;
    if (OU_variance != 0)
    {
        stochastic_conc_ = Gene::normal_(
            g_rng,
            std::normal_distribution<double>::param_type(
                conc_, sqrt(OU_variance)));
        if (stochastic_conc_ < 0)
            stochastic_conc_ = 0;
    }
    else
    {
        stochastic_conc_ = conc_;
    }
    return stochastic_conc_;
}


// Draw a new value for stochastic concentration
double Gene::update_stochastic_conc_gamma()
{
    stochastic_conc_ = (int)(0.5
			     + Gene::gamma_(
				 g_rng, std::gamma_distribution<double>::param_type(
				     stochastic_shape_, stochastic_scale_)));
    return stochastic_conc_;
}


//
double Gene::update_stochastic_conc_OU()
{
    // std::cout << "stochastic conc: " << stochastic_conc_;
    stochastic_conc_ = conc_
        + (stochastic_conc_ - conc_)*exp(-1/OU_tau_)
        + (sqrt((OU_diffusion_ * OU_tau_ / 2)*(1-exp(-2/OU_tau_)))
           * Gene::normal_(g_rng,
                           std::normal_distribution<double>::param_type(0, 1)));
    // std::cout << "; new conc: " << stochastic_conc_ << std::endl;
    if (stochastic_conc_ < 0)
        stochastic_conc_ = 0;
    return stochastic_conc_;
}
