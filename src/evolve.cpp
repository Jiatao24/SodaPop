#include <unistd.h>
#include <sstream>

#include <tclap/CmdLine.h>

#include <gsl/gsl_randist.h>

#include "nlohmann/json.hpp"

#include "rng.h"
#include "global.h"
#include "PolyCell.h"

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


github branch tmp_resistance_fixation_probability
 */

using json = nlohmann::json;

//! Find all genotypes and their population levels.
void countGenotypes(std::vector<PolyCell> cells,
                    std::map<std::string, unsigned int>& genotypeCounts)
{
    for (auto& cell : cells)
    {
        std::string genotype = cell.get_genotype();
        auto search = genotypeCounts.find(genotype);
        if (search == genotypeCounts.end())
        {
            genotypeCounts.insert({genotype, 1});
        }
        else
        {
            search->second++;
        }
    }
    return;
}


//! Save a json-based snapshot of genotype populations.
void saveSnapshot(char* buffer, std::vector<PolyCell>& cells,
                  int generationNumber)
{
    // Do genotype counting
    std::map<std::string, unsigned int> genotypeCounts;
    countGenotypes(cells, genotypeCounts);

    // We also want fitnesses
    // We do repeat ourselves here
    // (fitnesses are also calculated in main)
    std::map<std::string, double> genotypeFitnesses;

    // Fill genotypeFitnesses and count genotype populations.
    for (auto& cell : cells)
    {
        std::string genotype = cell.get_genotype();
        double fitness = cell.fitness();
        auto search = genotypeFitnesses.find(genotype);
        if (search == genotypeFitnesses.end())
            genotypeFitnesses.insert({genotype, fitness});
        else
            search->second += fitness;
    }

    // Get average for each genotype by dividing by count
    for (auto& it : genotypeFitnesses)
    {
        it.second /= (double)genotypeCounts[it.first];
    }

    // Calculate the overall average fitness.
    double averageFitness(0);
    for (auto& it : genotypeFitnesses)
    {
        averageFitness += it.second * genotypeCounts[it.first];
    }
    averageFitness /= (double)cells.size();

    // Let's get protein abundances.
    int gene_count = cells.front().gene_count();
    std::vector<double> abundances(gene_count, 0);
    std::map<std::string, std::vector<double>> genotypeAbundances;

    // We get abundance data only if fitness is stochastic concentration.
    if (PolyCell::ff_ == 8)
    {
        for (auto& cell : cells)
        {
            std::vector<double> cell_abundances = cell.get_abundances();
            std::string genotype = cell.get_genotype();

            auto search = genotypeAbundances.find(genotype);
            if (search == genotypeAbundances.end())
            {
                genotypeAbundances.insert(
                    {genotype, std::vector<double>(gene_count, 0)});
            }

            // Instead of getting mean by dividing by cell count
            //   later, just divide each value by count right now.
            //   who knows if it's worse performance anyway?
            for (int i=0; i<gene_count; i++)
            {
                abundances[i] += cell_abundances[i] / (double)cells.size();
                genotypeAbundances[genotype][i] +=
                    cell_abundances[i] / (double)genotypeCounts[genotype];
            }
        }
    }

    // Open snapshot file
    std::fstream OUT2(buffer, std::ios::out);
    if (!OUT2.is_open())
    {
         std::cerr << "Output file could not be opened";
         exit(1);
    }

    json output;
    output["generation"] = generationNumber;
    output["population"] = cells.size();
    output["drug level (nM)"] = DRUG_CONCENTRATION;
    output["average fitness"] = averageFitness;
    output["genotypes"] = genotypeCounts;
    output["genotype fitnesses"] = genotypeFitnesses;
    if (PolyCell::ff_ == 8)
    {
        output["average abundances"] = abundances;
        output["genotype abundances"] = genotypeAbundances;
    }

    // Here we output.
    OUT2 << std::setw(2) << output << std::endl;
    OUT2.close();   
    
    return;
}


// MAIN MAIN MAIN
int main(int argc, char *argv[])
{
    // These variables will hold simulation parameters, some of which
    //   are modifiable by the user through command line arguments.
    int generationNumber = 0;
    int generationMax = generationNumber + 1;
    int mutationCount = 0;
    int equilibrationGens = 0;
    double populationSize=1;
    int DT = 1;
    double TIME = 0;            // Hmmm, This is never updated?
    char buffer[200];
    bool mutateArgSet = false;
    bool trackMutations = false;

    std::string geneListFile, genesPath;
    std::string outDir, startSnapFile, matrixFile;

    int mutateGNum = 0, mutateResid = 1;
    std::string mutateResname;

    // Wrap everything in a try block
    // errors in input are caught and explained to user
    try
    { 
        // Define the command line object
        TCLAP::CmdLine cmd("SodaPop: a multi-scale model of molecular evolution", ' ', "v1.0");

        // Define value arguments
        TCLAP::ValueArg<int> maxArg("m","maxgen","Maximum number of generations",false,10000,"int");
        TCLAP::ValueArg<int> popArg("n","size","Initial population size",false,1,"int");
        TCLAP::ValueArg<int> dtArg("t","dt","Time interval for snapshots",false,1,"int");

        //files
        TCLAP::ValueArg<std::string> prefixArg("o","prefix","Prefix to be used for snapshot files",false,"sim","filename");
        TCLAP::ValueArg<std::string> geneArg("g","gene-list","Gene list file",true,"null","filename");
        TCLAP::ValueArg<std::string> startArg("p","pop-desc","Population description file",false,"null","filename");
        TCLAP::ValueArg<std::string> libArg("l","gene-lib","Gene library directory",false,"files/genes/","filename");

        TCLAP::ValueArg<std::string> matrixArg("i","input","Input file defining the fitness landscape",false,"null","filename");
    
        // fitness function
        TCLAP::ValueArg<int> fitArg("f","fitness","Fitness function",false,1,"integer ID");
    
        // // boolean switch to use DDG as input type
        // TCLAP::ValueArg<std::string> inputArg("","sim-type","Define simulation type\n<s> (from selection coefficient, DMS or otherwise)\n<stability> (from DDG matrix or distribution)", false,"s","string");
    
        // boolean switch to track mutations
        TCLAP::SwitchArg eventsArg("e","track-events","Track mutation events", cmd, false);
    
        // RNG seed
        TCLAP::ValueArg<unsigned long> seedArg("", "seed", "Seed value for RNG.", false, 0, "unsigned int (64-bit)");

        // RNG stream
        TCLAP::ValueArg<unsigned long> streamArg(
            "", "stream", "Stream value for RNG (set different streams for different runs to assure different "
            "random numbers for each run. But no two runs with the same stream will be identical unless the same seed "
            "is used as well.)",
            false, 0, "unsigned int (64-bit)");

        // X_FACTOR
        TCLAP::ValueArg<double> xfactorArg("x", "x-factor", "Enzymatic output factor.", false, 0, "positive double");

        TCLAP::ValueArg<double> concentrationArg("", "concentration", "Amount of drug present (nM)", false, 0, "non-negative double");

        TCLAP::SwitchArg rampingArg("", "ramping", "Drug concentration adjusts to population fitness.", cmd, false);

        TCLAP::ValueArg<int> equilArg("", "equil", "Time before mutation", false, 0, "nonnegative int");

        TCLAP::ValueArg<std::string> mutateArg("", "mutate", "Single mutation to introduce", false, "", "string consisting of \"g_num:resid:resname\", e.g. \"0:21:L\"");


        // Add the arguments to the CmdLine object.
        cmd.add(seedArg);
        cmd.add(streamArg);
        cmd.add(maxArg);
        cmd.add(popArg);
        cmd.add(dtArg);
        cmd.add(prefixArg);
        cmd.add(geneArg);
        cmd.add(startArg);
        cmd.add(libArg);
        cmd.add(fitArg);
        cmd.add(matrixArg);
        cmd.add(xfactorArg);
        cmd.add(concentrationArg);
        cmd.add(equilArg);
        cmd.add(mutateArg);

        // Parse the argv array.
        cmd.parse(argc, argv);

        // Get values from args. 
        generationMax = maxArg.getValue();
        populationSize = popArg.getValue();
        DT = dtArg.getValue();

        geneListFile = geneArg.getValue();
        outDir = prefixArg.getValue();
        startSnapFile = startArg.getValue();
        genesPath = libArg.getValue();

        equilibrationGens = equilArg.getValue();

        if (seedArg.isSet())
            setRngSeed(seedArg.getValue());

        if (streamArg.isSet())
            setRngStream(streamArg.getValue());

        std::cout << "Begin ... " << std::endl;
        // // even though we won't use the matrix
        // std::cout << "Initializing matrix ..." << std::endl;
        // InitMatrix();
        // // even though we don't need primordial genes.
        // std::cout << "Loading primordial genes file ..." << std::endl;
        // LoadPrimordialGenes(geneListFile, genesPath);
        PolyCell::ff_ = fitArg.getValue();

        trackMutations = eventsArg.getValue();

        if (xfactorArg.isSet())
        {
            X_FACTOR = xfactorArg.getValue();
            if (X_FACTOR <= 0)
            {
                std::cerr << "error: --x-factor argument not positive ("
                          << X_FACTOR << ")\n";
                exit(1);
            }
        }

        if (concentrationArg.isSet())
        {
            DRUG_CONCENTRATION = concentrationArg.getValue();
            if (DRUG_CONCENTRATION < 0)
            {
                std::cerr << "error: --concentration argument < 0 ("
                          << DRUG_CONCENTRATION << ")\n";
                exit(1);
            }
        }

        if (mutateArg.isSet())
        {
            mutateArgSet = true;
            std::istringstream stream(mutateArg.getValue());
            std::string temp;
            std::getline(stream, temp, ':');
            mutateGNum = std::stoi(temp);
            std::getline(stream, temp, ':');
            mutateResid = std::stoi(temp);
            std::getline(stream, temp, ':');
            mutateResname = temp;

            std::cout << "--mutate=\"" << mutateArg.getValue() << "\":" << std::endl;
            std::cout << "  At generation " << equilibrationGens
                      << ", mutate a single cell. Mutate gene num "
                      << mutateGNum << ", at residue "
                      << mutateResid << ", to residue "
                      << mutateResname << "." << std::endl;
        }

    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId()
                  << std::endl;
    }

    std::cout << "Opening starting population snapshot ..." << std::endl;
    std::fstream startsnap(startSnapFile.c_str(), std::ios::in|std::ios::binary);
    if (!startsnap.is_open())
    {
        std::cerr << "File could not be open: "<< startSnapFile << std::endl;
        exit(2);
    }
    
    // header
    int Total_Cell_Count;
    double frame_time;
    startsnap.read((char*)(&frame_time), sizeof(double));
    startsnap.read((char*)(&TIME), sizeof(double));
    startsnap.read((char*)(&Total_Cell_Count), sizeof(int));

    sprintf(buffer, "%s/snapshots", outDir.c_str());
    std::string outPath = buffer;
    std::cout << "Creating directory " << outPath << " ... "
              << (makePath(outPath) ? "OK" : "failed") << std::endl;

    std::vector<PolyCell> Cell_arr;
    std::vector<double> fitnesses(populationSize, 0.);
    // double w_sum = 0;

    // CREATE VECTOR WITH populationSize CELLS
    std::cout << "Creating a population of " << populationSize << " cells ..."
              << std::endl;
    PolyCell A(startsnap, genesPath);
    Cell_arr = std::vector<PolyCell>(populationSize, A);
    for (auto it = Cell_arr.begin(); it != Cell_arr.end(); ++it)
    {
        it->ch_barcode(getBarcode());
        it->init_gene_stochastic_concentrations();
        it->UpdateRates();
    } 
    startsnap.close();

    std::cout << "Saving initial genotype counts ... " << std::endl;
    // save initial population snapshot
    sprintf(buffer, "%s/%s.gen%010d.json", outPath.c_str(),
            outDir.c_str(), generationNumber); 
    saveSnapshot(buffer, Cell_arr, (int)frame_time);

    std::ofstream MUTATIONLOG;
    if (trackMutations)
    {
        // Open MUTATION LOG
        sprintf(buffer, "out/%s/MUTATION_LOG", outDir.c_str());
        MUTATIONLOG.open(buffer);
        if (!MUTATIONLOG.is_open())
        {
            std::cerr << "Mutation log file could not be opened";
            exit(1);
        }
    }
    
    std::cout << "Starting evolution ..." << std::endl;

    // Setup GSL-style RNG for multinomial distribution function.
    gsl_rng *r_gsl;
    r_gsl = gsl_rng_alloc(&gsl_rng_pcg);

    // WRIGHT-FISHER PROCESS
    while (generationNumber < generationMax)
    {
        // Vector to hold the next generation
        std::vector<PolyCell> Cell_temp;
        std::vector<unsigned int> progeny_per_cell(populationSize, 0);

        // Iterate through cells and gather fitnesses into vector.
        auto fitnesses_it = fitnesses.begin();
        for (auto cell_it = Cell_arr.begin(); cell_it != Cell_arr.end();
             ++cell_it, ++fitnesses_it)
            *fitnesses_it = cell_it->fitness();

        // Number of progeny per cell determined via sample from
        // multinomial distribution.
        gsl_ran_multinomial(r_gsl, populationSize, populationSize,
                            fitnesses.data(), progeny_per_cell.data());

        // Iterate through cells and number of progeny, filling Cell_temp.
        auto cell_it = Cell_arr.begin();
        for (auto& progeny : progeny_per_cell)
        {
            std::fill_n(std::back_inserter(Cell_temp), progeny, (*cell_it));
            cell_it++;
        }
        
        // update generation counter
        generationNumber++; 
     
        // Make the single mutation if it's time to do so.
        if (mutateArgSet && generationNumber == equilibrationGens)
        {
            auto& selectedCell = Cell_temp[(int)(randomNumber() * populationSize)];
            selectedCell.mutGene(mutateGNum, mutateResid, mutateResname);
            std::cout << "At generation " << generationNumber
                      << ", mutated a single cell." << std::endl;
        }

        // Now go through each cell to update fitness (stochastic gene expression)
        for (auto& cell : Cell_temp)
        {
            // Even if we don't mutate, update the fitness.
            // In stochastic gene expression, fitness will change.
            cell.UpdateRates();
        }

        Total_Cell_Count = (int)(Cell_temp.size());
        assert (Total_Cell_Count == populationSize);
        // swap population with initial vector
        Cell_arr.swap(Cell_temp);

        // save population snapshot every DT generations
        if ((generationNumber % DT) == 0)
        {
             sprintf(buffer, "%s/%s.gen%010d.json", outPath.c_str(),
                     outDir.c_str(), generationNumber); 

             saveSnapshot(buffer, Cell_arr, generationNumber);
             std::cout << "Generation: " << generationNumber << std::endl;
        }

        // Stop simulation if
        if (generationNumber > equilibrationGens)
        {
            std::map<std::string, unsigned int> genotypeCounts;
            countGenotypes(Cell_arr, genotypeCounts);
            if (genotypeCounts.size() == 1)
            {
                std::cout << "Number of genotypes == 1. Stopping simulation at "
                          << generationNumber << " generations" << std::endl;
                sprintf(buffer, "%s/%s.gen%010d.json", outPath.c_str(),
                        outDir.c_str(), generationNumber); 

                saveSnapshot(buffer, Cell_arr, generationNumber);
                std::cout << "Final generation: " << generationNumber << std::endl;
                break;
            }
        }

    }

    MUTATIONLOG.close();
    std::cout << "Done." << std::endl;
    std::cout << "Total number of mutation events: " << mutationCount << std::endl;

    gsl_rng_free(r_gsl);
    return 0;
}
