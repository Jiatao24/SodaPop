#include <unistd.h>

#include <tclap/CmdLine.h>

#include <gsl/gsl_randist.h>

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
 */

int main(int argc, char *argv[])
{
    // these variables will hold the parameters input (or not) by the user
    int generationNumber = 0;
    int generationMax = generationNumber + 1;
    int mutationCount = 0;
    double populationSize=1;
    int DT = 1;
    double TIME = 0;            // This is never updated?
    char buffer[200];
    bool enableAnalysis = false;
    bool trackMutations = false;
    bool createPop = false;
    bool useShort = false;

    std::string geneListFile, genesPath;
    std::string outDir, startSnapFile, matrixFile;

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

        //use gamma distribution to draw selection coefficients
        TCLAP::SwitchArg gammaArg("","gamma","Draw selection coefficients from gamma distribution", cmd, false);

        //use normal distribution to draw selection coefficients
        TCLAP::SwitchArg normalArg("","normal","Draw selection coefficients from normal distribution", cmd, false);

        //first parameter of distribution
        TCLAP::ValueArg<double> alphaArg("","alpha","Alpha parameter of distribution\nGamma -> shape\nNormal -> mean",false,1,"double");

        //second parameter of distribution
        TCLAP::ValueArg<double> betaArg("","beta","Beta parameter of distribution\nGamma -> scale\nNormal -> S.D.",false,1,"double");

        // boolean switch to create population from scratch
        TCLAP::SwitchArg initArg("c","create-single","Create initial population on the fly", cmd, false);
    
        // boolean switch to enable analysis
        TCLAP::SwitchArg analysisArg("a","analysis","Enable analysis scripts", cmd, false);
    
        // boolean switch to track mutations
        TCLAP::SwitchArg eventsArg("e","track-events","Track mutation events", cmd, false);
    
        // boolean switch to use short format for snapshots
        TCLAP::SwitchArg shortArg("s","short-format","Use short format for population snapshots", cmd, false);

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

        // Add the arguments to the CmdLine object.
        cmd.add(seedArg);
        cmd.add(streamArg);
        cmd.add(xfactorArg);
        cmd.add(maxArg);
        cmd.add(popArg);
        cmd.add(dtArg);
        cmd.add(prefixArg);
        cmd.add(geneArg);
        cmd.add(startArg);
        cmd.add(libArg);
        cmd.add(fitArg);
        cmd.add(matrixArg);
        cmd.add(alphaArg);
        cmd.add(betaArg);

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

        if (seedArg.isSet())
            setRngSeed(seedArg.getValue());

        if (streamArg.isSet())
            setRngStream(streamArg.getValue());

        std::cout << "Begin ... " << std::endl;
        // if (inputType == "s")
        // {
        //     PolyCell::fromS_ = true;
        //     PolyCell::ff_ = 4;
        //     std::cout << "Initializing matrix ..." << std::endl;
        //     InitMatrix();
        //     std::cout << "Loading primordial genes file ..." << std::endl;
        //     LoadPrimordialGenes(geneListFile,genesPath);
        //     // if matrix is given
        //     if(matrixArg.isSet())
        //     {
        //         matrixFile = matrixArg.getValue();
        //         std::cout << "Extracting DMS matrix ..." << std::endl;
        //         ExtractDMSMatrix(matrixFile.c_str());
        //     }
        //     else
        //     {
        //         PolyCell::useDist_ = true;
        //         if (gammaArg.isSet())
        //         {
        //             double shape = alphaArg.getValue();
        //             double scale = betaArg.getValue();
        //             Gene::setGammaParams(shape, scale);
        //         }
        //         else if (normalArg.isSet())
        //         {
        //             double mean = alphaArg.getValue();
        //             double stddev = betaArg.getValue();
        //             Gene::setNormalParams(mean, stddev);
        //         }
        //     }
        // }
        // else if (inputType == "stability")
        // {
        //     std::cout << "Initializing matrix ..." << std::endl;
        //     InitMatrix();
        //     std::cout << "Loading primordial genes file ..." << std::endl;
        //     LoadPrimordialGenes(geneListFile,genesPath);
        //     PolyCell::ff_ = fitArg.getValue();
        //     // if DDG matrix is given
        //     if (matrixArg.isSet())
        //     {
        //         matrixFile = matrixArg.getValue();
        //         std::cout << "Extracting PDDG matrix ..." << std::endl;
        //         ExtractPDDGMatrix(matrixFile.c_str());
        //     }
        //     else
        //     {
        //         PolyCell::useDist_ = true;
        //     }
        // }

        enableAnalysis = analysisArg.getValue();
        trackMutations = eventsArg.getValue();
        useShort = shortArg.getValue();
        createPop = initArg.getValue();

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

    sprintf(buffer, "out/%s/snapshots", outDir.c_str());
    std::string outPath = buffer;
    std::cout << "Creating directory " << outPath << " ... "
              << (makePath(outPath) ? "OK" : "failed") << std::endl;

    std::vector<PolyCell> Cell_arr;
    std::vector<double> fitnesses(populationSize, 0.);
    // double w_sum = 0;

    // IF POPULATION IS INITIALLY MONOCLONAL
    // CREATE VECTOR WITH populationSize CELLS
    if (createPop)
    {
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
    }
    else
    {
        // ELSE IT MUST BE POPULATED CELL BY CELL FROM SNAP FILE
        // Note number of cells is Total_Cell_Count
        Cell_arr.reserve(populationSize);
        int count = 0;
        std::cout << "Constructing population from source "
                  << startSnapFile.c_str() << " ..." << std::endl;
        while (count < Total_Cell_Count && !startsnap.eof())
        {
            Cell_arr.emplace_back(startsnap, genesPath);
            count++;
        }
    }
    startsnap.close();

    std::cout << "Saving initial population snapshot ... " << std::endl;
    // save initial population snapshot
    sprintf(buffer, "%s/%s.gen%010d.snap", outPath.c_str(),
            outDir.c_str(), generationNumber); 

    // Open snapshot file
    std::fstream OUT2(buffer, std::ios::out | std::ios::binary);
    if (!OUT2.is_open())
    {
         std::cerr << "Snapshot file could not be opened";
         exit(1);
    }

    Total_Cell_Count = Cell_arr.size();
    OUT2.write((char*)(&frame_time), sizeof(double));
    OUT2.write((char*)(&TIME), sizeof(double));
    OUT2.write((char*)(&Total_Cell_Count), sizeof(int));

    if (useShort)
    {
        for (auto it = Cell_arr.begin(); it != Cell_arr.end(); ++it)
        {
            it->dumpShort(OUT2);
        } 
    }
    else
    {
        int l=1;
        // dump snapshot of initial population and get sum of fitnesses
        for (auto it = Cell_arr.begin(); it != Cell_arr.end(); ++it)
        {
            it->dump(OUT2, l);
            l++;
        } 
    }
    
    OUT2.close();   

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

    // Setup GSL-style RNG for multinomial.
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
        
        // Now go through each cell for mutation
        for (auto& cell : Cell_temp)
        {
            if (cell.mrate() * cell.genome_size() > randomNumber())
            {
                mutationCount++;
                if (trackMutations)
                {
                    // mutate and write mutation to file
                    cell.ranmut_Gene(MUTATIONLOG, generationNumber);
                }
                else
                {
                    cell.ranmut_Gene();
                }       
            }
            else
            {
                // Even if we don't mutate, update the fitness.
                // In stochastic gene expression, fitness will change.
                cell.UpdateRates();
            }
        }

        Total_Cell_Count = (int)(Cell_temp.size());
        assert (Total_Cell_Count == populationSize);
        // swap population with initial vector
        Cell_arr.swap(Cell_temp);

        // update Ns and Na for each cell
        for (auto it = Cell_arr.begin(); it != Cell_arr.end(); ++it)
        {
            it->UpdateNsNa();
        }
        
        // update generation counter
        generationNumber++; 
        // save population snapshot every DT generations
        if ((generationNumber % DT) == 0)
        {
             sprintf(buffer, "%s/%s.gen%010d.snap", outPath.c_str(),
                     outDir.c_str(), generationNumber); 

             // Open snapshot file
             std::fstream OUT2(buffer, std::ios::out | std::ios::binary);
             if (!OUT2.is_open())
             {
                 std::cerr << "Snapshot file could not be opened";
                 exit(1);
             }
      
             double frame_time = generationNumber;
             OUT2.write((char*)(&frame_time), sizeof(double));
             OUT2.write((char*)(&TIME), sizeof(double)); // this is always 0?
             OUT2.write((char*)(&Total_Cell_Count), sizeof(int));

             if (useShort)
             {
                 for (auto it = Cell_arr.begin(); it != Cell_arr.end(); ++it)
                 {
                     it->dumpShort(OUT2);
                 } 
             }
             else
             {
                 int l=1;
                 for (auto it = Cell_arr.begin(); it != Cell_arr.end(); ++it)
                 {
                     it->dump(OUT2, l);
                     l++;
                 }
             }
              
             OUT2.close();
         }
    }

    MUTATIONLOG.close();
    std::cout << "Done." << std::endl;
    std::cout << "Total number of mutation events: " << mutationCount << std::endl;

    // if the user toggled analysis, call shell script
    if (enableAnalysis)
    {
        std::string script = "/path/to/barcodes.sh";
        std::string command = script+" "+outDir+" "+std::to_string(generationMax)+" "+std::to_string(populationSize)+" "+std::to_string(DT)+" "+std::to_string((int) useShort);
        std::cout << "Call analysis using the following command:" << std::endl;
        std::cout << command << std::endl;
        // const char *cmd = command.c_str();
        // system(cmd);
    }

    gsl_rng_free(r_gsl);
    return 0;
}
