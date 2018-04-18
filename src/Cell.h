// -*- c++ -*-
#ifndef CELL_H
#define CELL_H

#include "global.h"
#include "Gene.h"

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

class Cell 
{
public:
    Cell();
    Cell(std::fstream&);			    
    Cell(std::fstream&, const std::string&);  

    void FillGene_L();
    int total_mutations(const int&);
    void init_gene_stochastic_concentrations();

    virtual void UpdateRates() = 0;
    virtual void dump(std::fstream&, int) = 0;
    virtual void PrintCell(int) = 0;       
             
    // Getters
    int ID() {return ID_;}
    double mrate() {return c_mrate_;}
    int gene_count() {return Gene_arr_.size();}
    int genome_size() {return Gene_L_.back();}
    std::string barcode() {return barcode_;}   

    // Setters
    void change_ID(int a) {ID_ = a;}
    void ch_barcode(std::string s) {barcode_ = s;}

protected:
    // Organism barcode
    std::string barcode_;
    // Organism ID
    int ID_;
    // Initial mutation rate.
    double o_mrate_;
    // Current mutation rate.
    double c_mrate_;	
        			
    // Array of genes.
    std::vector<Gene> Gene_arr_;
    // Cumulative sum of gene lengths (i.e. genome size).
    VectInt Gene_L_;
};


#endif  // CELL_H
