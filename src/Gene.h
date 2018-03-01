// -*- c++ -*-
#ifndef GENE_H
#define GENE_H
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
public:
    Gene();
    Gene(int, std::string, double);
    Gene(std::fstream&);
    Gene(const Gene&); // Copy constructor (is this necessary?)
  
    bool operator==(Gene&);
    Gene& operator=(const Gene&); // is this necessary?

    static void setGammaParams(double, double);
    static void setNormalParams(double, double);
    static double randomGamma();
    static double randomNormal();
  
    double Mutate_Stabil_Gaussian(int, int);
    std::string Mutate_Stabil(int, int);
    double Mutate_Select_Dist(int, int);
    std::string Mutate_Select(int, int);

    void Update_Sequences(std::string);
 
    void ch_dg(const double a) {dg_ = a;}
    void ch_conc(const double c) {conc_ = c;}
    void ch_f(const double a) { f_ = a; }
    void ch_Na(const int a) { Na_ = a; }
    void ch_Ns(const int a) { Ns_ = a; }

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
    double update_stochastic_conc();

    double e() const {return e_;}

    double CheckDG();
    double Pnat();
    double functional(bool stochastic=false);
    double misfolded(bool stochastic=false);
    
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

    static std::gamma_distribution<> gamma_;
    static std::normal_distribution<> normal_;

    double e_;		//essentiality: 1-if directly involved in replication, 0-otherwise
    // i don't see where e is being used however. also why is it a double?
   
};

#endif  // GENE_H
