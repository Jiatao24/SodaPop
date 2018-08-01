// -*- c++ -*-
#ifndef GENE_H
#define GENE_H

#include <random>
#include <unordered_map>

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

typedef std::pair<double, double> biophysicalParam;

//! Represents a single protein-encoding gene.
class Gene 
{
public:
    //! Simple constructor setting most values to 0.
    Gene();

    /*! Constructor
     *
     * @param gene_in  Should be file stream to text file with gene parameters.
     */
    Gene(std::fstream& gene_in);

    // //! Copy constructor (is this necessary?)
    // Gene(const Gene&); 
  
    bool operator==(Gene&);
    // Gene& operator=(const Gene&); // is this necessary?

    // Static functions to set probability distribution params or draw
    // random value.
    static void setGammaParams(double, double);
    static void setNormalParams(double, double);
    static double randomGamma();
    static double randomNormal();
  
    // Gene mutation functions

    //! Mutate, updating biophysical enzyme parameters
    std::string mutate(int site, int bp);

    //! Draw DDG from Normal(1, 1.7)
    double Mutate_Stabil_Gaussian(int, int);
    //! Draw DDG from DDG matrix.
    std::string Mutate_Stabil(int, int);

    //! Change DNA sequence to new input.
    void Update_Sequences(std::string DNAsequence);
 
    //! Parse json file of biophysical parameters.
    void loadBiophysicalParameters(std::string path);
    //! Determine the biophysical genotype of the current sequence.
    void identifyGenotype();

    // Basically a bunch of setters for private members
    void ch_dg(const double a) {dg_ = a;}
    void ch_conc(const double c) {conc_ = c;}
    void ch_f(const double a) { f_ = a; }
    void ch_Na(const int a) { Na_ = a; }
    void ch_Ns(const int a) { Ns_ = a; }

    // Basically a bunch of getters for private members
    int num() const {return g_num_;}
    int length() const {return ln_;}
    int AAlength() const {return la_;}
    int Ns() const {return Ns_;}
    int Na() const {return Na_;}
    std::string nseq() const {return nucseq_;}
    double dg() const {return dg_;}
    double f() const {return f_;}

    double conc() const {return conc_;}
    double stochastic_conc() const {return stochastic_conc_;}

    //! Pick a random value for stochastic concentration.
    double init_stochastic_conc();
    //! Draw new protein concentration value from gamma distribution
    double update_stochastic_conc_gamma();
    //! Update protein concentration value via Ornstein-Uhlenbeck process
    double update_stochastic_conc_OU();

    double e() const {return e_;}

    double CheckDG();           // this function is not implemented?
    double Pnat();

    //! Returns effective (functional) protein abundance.
    double functional(bool stochastic=false);
    double misfolded(bool stochastic=false);
    //! Returns enzymatic output, which depends on abundance and
    //! biophysical catalytic properties.
    double enzymaticFlux(bool stochastic=false);
    
private:
    int g_num_;		// numeric ID pointing to primordial gene
    int ln_;		// length nuc seq
    int la_;		// length aa seq

    int Na_;		// number of non-synonymous substitutions
    int Ns_;		// number of sysnonymous substitutions
    
    std::string nucseq_;	// nucleotide sequence
    
    double dg_;		// stability
    double f_;      // gene "fitness"

    // Concentration
    double conc_;
    double stochastic_conc_;
    double stochastic_shape_;
    double stochastic_scale_;
    double OU_tau_ = 0;
    double OU_diffusion_ = 0;

    // Biophysical parameters
    std::unordered_map<std::string, biophysicalParam> biophysicalParameters_;
    std::vector<int> keyResidueNumbers_;
    std::string genotype_;
    double drug_penetration_factor_;
    double rel_kcat_over_km_;
    double k_i_;

    static std::gamma_distribution<> gamma_;
    static std::normal_distribution<> normal_;

    double e_;		//essentiality: 1-if directly involved in replication, 0-otherwise
    // i don't see where e is being used however. also why is it a double?
   
};

#endif  // GENE_H
