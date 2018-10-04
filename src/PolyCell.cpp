#include <iostream>
#include <iomanip>
#include <stdio.h>

#include "PolyCell.h"
#include "rng.h"

// By default the fitness function is set to neutral.
int PolyCell::ff_ = 5;
bool PolyCell::useDist_ = false;
bool PolyCell::fromS_ = false;


PolyCell::PolyCell()
{
}

// Construct from cell file (text description).
PolyCell::PolyCell(std::fstream& f) : Cell(f)
{
    selectFitness();
    // Update current rates
    UpdateRates(0);  
    // Fill gene length array
    FillGene_L();
}    

// Construct from cell file in binary and gene path.
PolyCell::PolyCell(std::fstream& f, const std::string& s) : Cell(f, s)
{
    selectFitness();
    // Update current rates
    UpdateRates(0);  
    // Fill gene length array
    FillGene_L();
}

void PolyCell::selectFitness()
{
    switch(PolyCell::ff_){
    // case 1: fit = &PolyCell::flux;
    //     break;
    // case 2: fit = &PolyCell::toxicity;
    //     break;
    // case 3: fit = &PolyCell::metabolicOutput;
    //     break;
    // case 4: fit = &PolyCell::multiplicative;
    //     break;
    // case 5: fit = &PolyCell::neutral;
    //     break;
    // case 6: fit = &PolyCell::stochasticExpression;
    //     break;
    case 7: fit = &PolyCell::enzymaticOutput;
        break;
    case 8: fit = &PolyCell::stochasticEnzymaticOutput;
        break;
    default:
        ;
    }
}

// // FLUX FITNESS FUNCTION
// double PolyCell::flux()
// {
//     double f = 0;
//     //sum (concentration*Pnat) over all genes
//     for (auto i = Gene_arr_.begin(); i != Gene_arr_.end(); ++i)
//     {
//         f += 1 / (i->functional());
//     }
//     return exp(A_FACTOR / f - 1);
// }

// // TOXICITY FITNESS FUNCTION
// double PolyCell::toxicity()
// {
//     double f = 0;
//     //sum (concentration*(1-Pnat)) over all genes
//     for (auto i = Gene_arr_.begin(); i != Gene_arr_.end(); ++i)
//     {
//         f +=i->misfolded();
//     }
//     double w = exp(-(COST*f));
//     return w;
// }

// // METABOLIC FLUX FITNESS FUNCTION
// // FLUX AND TOXICITY COMBINED
// double PolyCell::metabolicOutput()
// {
//     double flux = 0;
//     double toxicity = 0;
//     for (auto& it : Gene_arr_) {
//         flux += 1 / it.functional();
//         toxicity += it.misfolded();
//     }

//     flux = A_FACTOR / flux;
//     toxicity = COST * toxicity;

//     double fitness = flux - toxicity;
//     return (fitness < 0) ? 0 : fitness;
// }

// // MULTIPLICATIVE FITNESS FUNCTION
// double PolyCell::multiplicative()
// {
//     double fitness = 1;
//     for (auto i = Gene_arr_.begin(); i != Gene_arr_.end(); ++i)
//     {
//         fitness *= i->f();
//         if (fitness < 0)
//             return 0;
//     }
//     return fitness;
// }

// // NEUTRAL FITNESS FUNCTION
// double PolyCell::neutral()
// {
//     return 1;
// }

// // STOCHASTIC PROTEIN EXPRESSION (ALSO BASED ON FLUX)
// // it's basically metabolicOutput with stochastic concentration
// double PolyCell::stochasticExpression()
// {
//     double flux = 0;
//     double toxicity = 0;
//     for (auto& it : Gene_arr_)
//     {
//         it.update_stochastic_conc_gamma();
//         if (it.stochastic_conc() == 0)
//         {
// 	    // Fatal to not express any protein
//             return 0;
//         }
//         flux += 1 / it.functional(true);
//         toxicity += it.misfolded(true);
//     }

//     flux = A_FACTOR / flux;
//     toxicity = COST * toxicity;
    
//     double fitness = flux - toxicity;
//     return (fitness < 0) ? 0 : fitness;
// }

// ENZYMATIC OUTPUT (ALSO BASED ON FLUX)
// (Presumes that there is only a single gene)
double PolyCell::enzymaticOutput()
{
    auto& first_gene = Gene_arr_.front();
    double fitness = first_gene.enzymaticFlux()
        / (X_FACTOR + first_gene.enzymaticFlux());
    return fitness;
}

// STOCHASTIC ENZYMATIC OUTPUT (ALSO BASED ON FLUX)
//(Presumes that there is only a single gene)
double PolyCell::stochasticEnzymaticOutput()
{
    auto& first_gene = Gene_arr_.front();
    double fitness = first_gene.enzymaticFlux(true)
        / (X_FACTOR + first_gene.enzymaticFlux(true));
    return fitness;
}

double PolyCell::fitness()
{
    return fitness_;
}

void PolyCell::UpdateRates(double dt)
{
    for (auto &gene : Gene_arr_)
    {
        gene.update_stochastic_conc_OU(dt);
    }
    fitness_ = (this->*fit)();
}

// ctr is the generation
void PolyCell::ranmut_Gene(std::ofstream& log, int ctr)
{
    // get genome size
    int L = Gene_L_.back();
    // pick random site to mutate
    int site = (int) ( L * randomNumber());

    // find the corresponding gene
    auto j = Gene_arr_.begin();
    auto k = Gene_L_.begin();

    if(site >= (*k)){
        // random number generated is greater than
        // the cummulative sum of genes
        for(k = Gene_L_.begin(); k != Gene_L_.end(); ++k){
            if( site<(*k) ) break;
            j++; 
        }        
        k--;
        site = site - (*k);        
    }

    std::string mutation = "";

    int bp = (int) (3 * randomNumber());

    double wi = fitness();

    mutation = j->mutate(site, bp);
    UpdateRates(0);

    double wf = fitness();
    double s = wf - wi;

    // save beneficial mutations to log
    // we could save all mutations with abs(s) >= some value x

    // Barcode, Gene number, FROM, residue number, TO, Fitness difference, Generation count
    if (mutation.size() > 0)
    {
        log << barcode().c_str() << "\t";
        log << std::fixed;
        log << mutation << "\t";
        log << s << "\t";		// fitness difference wf - wi
        log << ctr << std::endl;		// generation count
    }
}

void PolyCell::ranmut_Gene()
{
    // get genome size
    int L = Gene_L_.back();
    // pick random site to mutate
    int site = (int) (L * randomNumber());

    // find the corresponding gene
    std::vector<Gene>::iterator j = Gene_arr_.begin();
    auto k = Gene_L_.begin();

    if (site >= (*k))
    {
        // random number generated is greater than
        // the cummulative sum of genes
        for (k = Gene_L_.begin(); k != Gene_L_.end(); ++k)
        {
            if (site < *k)
                break;
            j++; 
        }        
        k--;
        site = site - (*k);        
    }

    int bp = (int) (3 * randomNumber());
    j->mutate(site, bp);
    UpdateRates(0);
}

// Dump cell information to binary file
// Ordering is
// cell_index, cell_id, barcode size, barcode, fitness, mutation rate (c_mrate), n_genes.
// Then for each gene:
// gene_id, essentiality, concentration, stochastic_conc, deltaG, gene fitness,
// n_nonsynonymous, n_synonymous, DNA sequence length, DNA sequence
void PolyCell::dump(std::fstream& OUT, int cell_index)
{
    int x;
    double y;

    OUT.write((char*)(&cell_index),sizeof(int));
    
    //cell ID
    OUT.write((char*)(&ID_),sizeof(int));

    int s = barcode().size();
    OUT.write((char*)&s, sizeof(int));

    OUT.write(barcode().c_str(), s);

    y = fitness();
    OUT.write((char*)(&y),sizeof(double));
    
    y = c_mrate_;		 	 
    OUT.write((char*)(&y),sizeof(double));

    x = (int)(Gene_arr_.size());		 	 
    OUT.write((char*)(&x),sizeof(int));

    for (auto gene_it = Gene_arr_.begin(); gene_it != Gene_arr_.end();
         ++gene_it)
    {
        int gene_nid = gene_it->num();
        double s = gene_it->e();
        double c = gene_it->conc();
	double stoc_conc = gene_it->stochastic_conc();
        double dg = -kT*log(gene_it->dg());
        double f = gene_it->f();

        int Ns = gene_it->Ns();
        int Na = gene_it->Na();

        OUT.write((char*)(&gene_nid),sizeof(int));
        OUT.write((char*)(&s),sizeof(double));
        OUT.write((char*)(&c),sizeof(double));
        OUT.write((char*)(&stoc_conc),sizeof(double));
        OUT.write((char*)(&dg),sizeof(double));
        OUT.write((char*)(&f),sizeof(double));
        OUT.write((char*)(&Na),sizeof(int));
        OUT.write((char*)(&Ns),sizeof(int));

        //Save length of nucleo sequence
        std::string DNAsequence = gene_it->nseq();
        int nl = DNAsequence.length();
        OUT.write((char*)&nl, sizeof(int));
        OUT.write(DNAsequence.data(), nl);
    }
}

// Dump cell summary to binary file
void PolyCell::dumpShort(std::fstream& OUT)
{
    int x;
    double y;

    int s = barcode().size();
    OUT.write((char*)&s, sizeof(int));

    OUT.write(barcode().c_str(), s);

    x = Na();
    OUT.write((char*)(&x),sizeof(int));

    x = Ns();
    OUT.write((char*)(&x),sizeof(int));

    y = fitness();         
    OUT.write((char*)(&y),sizeof(double));
}

// Print cell information to stdout
void PolyCell::PrintCell(int cell_ndx)
{
    char buffer[140];
    sprintf(buffer,"C %6d %6d %12e %12e %d", cell_ndx, ID_, o_mrate_, mrate(), (int)Gene_arr_.size());  
    std::cout << buffer << std::endl;
    for (auto gene_it = Gene_arr_.begin(); gene_it != Gene_arr_.end();
         ++gene_it)
    {
        //cout << "X ";
        int gene_nid = gene_it->num();
        double e = gene_it->e();
        double c = gene_it->conc();
        double dg = -kT * log(gene_it->dg());
        int Ns = gene_it->Ns();
        int Na = gene_it->Na();
           
        sprintf(buffer,"G %d% 2.2f %10.8f %10.8f %d %d ", gene_nid, e, c, dg, Ns, Na);
        std::cout << buffer << std::endl;  
    }
    std::cout << std::endl;
}

void PolyCell::UpdateNsNa()
{
    int new_Na = 0;
    int new_Ns = 0;
    for (auto gene_it = Gene_arr_.begin(); gene_it != Gene_arr_.end();
         ++gene_it)
    {
        new_Na += gene_it->Na();
        new_Ns += gene_it->Ns();
    }
    Total_Na_ = new_Na;
    Total_Ns_ = new_Ns;
}
