
#include <iostream>
#include <map>
#include <fstream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <functional>
using namespace std;
using namespace std::chrono;
using std::transform;

//constants

//const int range =200;
//const int Nsweeps = 100;
//const int sweepsPerMag = 1;   //make sure that Nsweeps / sweepsPerMag is an integer or code will break!!!
//const double J = 0.2; //interaciton strength for neighbours (1 for milestone)
//const double B = -0.05; //external mag field strength for whole lattice
//const double kbt = 1;  // meV   26 is from the notes for room temp, lower means lower partition means less chance of random fluctuation
//const int particleN0 = pow(range, 2);
//const int dataPoints = 1 + Nsweeps / sweepsPerMag;





//if particle flips to be parallel with adjacent spins, then -J.S1.S2 is  negative as S1.S2 is positive, 
// this means E is more negative which is favorable
// positive J means particles tend to become parallel
// B decides which direction they tend to
// whichever system has lowest potential (energy) is what it tends to


void saveToFile(string filename, vector <double> xAxis, vector <double> yAxis,vector <double> yErrors, int dataPoints, string xLabel, string yLabel)
{
    ofstream file(filename);
    file << dataPoints;
    file << "\n";
    file << xLabel;
    file << "\n";
    file << yLabel;

    for (int i = 0; i < dataPoints; i++)
    {
        file << "\n";
        file << xAxis[i];
        file << ",";
        file << yAxis[i];
        file << ",";
        file << yErrors[i];
    }
    file.close();
}

double findMean(vector<double> arr, int cutoff, bool beforeCutoff) //cutoff inclusive, eg cutoff 3 means up to and including the index 3
{
    //cutoff is indec for the last item to be included ( if true)
    //therefore to do whole array need to set cutoff as array.size() - 1
    double mean = 0;
    if (beforeCutoff)
    {
        for (int i = 0; i <= cutoff; i++)
        {
            mean += arr[i];
        }
        mean = mean / (cutoff + 1);
    }
    else
    {
        for (int i = cutoff; i < arr.size(); i++)
        {
            mean += arr[i];
        }
        mean = mean / (arr.size() - cutoff);
    }
   
    return mean;
}

double findStdDeviation(vector<double> arr, int cutoff, bool before)
{
    double mean = findMean(arr, cutoff, before);
    double stdDev = 0;
    if (before)
    {
        for (int i = 0; i <= cutoff; i++)
        {

            stdDev += pow((arr[i] - mean), 2);
        }
        stdDev = sqrt((stdDev / cutoff));
    }
    else
    {
        for (int i = cutoff; i < arr.size(); i++)
        {

            stdDev += pow((arr[i] - mean), 2);
        }

        stdDev = sqrt((stdDev / (arr.size() - cutoff -1)));
    }
    
    return stdDev;
        
}

double findStdError(vector<double> arr, int cutoff, bool before)
{
    return findStdDeviation(arr, cutoff, before) / cutoff;
}  //assumes whole array

void produceGraphs(string dataFileName)
{
    string filename = "ProduceGraphs.py ";
    string command = "Python ";
    command += filename;
    command += dataFileName;
    system(command.c_str());
}

class particle
{
    
public:
    particle() {};
    int spin = 0;
    int neighbours[4] = { 0,0,0,0 };

    particle(int identifier, int range)
    {
        setInitialSpin();
        setNeighbours(identifier, range);
    }
         
    void setInitialSpin()
    {
        spin = rand() % 2; //assign a random spin to the particle
        if (spin != 1)
        {
            spin = -1;
        }
    } 

    void setNeighbours(int ID, int range)
    {
        int temporaryID = ID + 1;
        if (temporaryID % range == 0) //check if on right
        {
            neighbours[0] = temporaryID - range + 1; //assign "right" neighbour to far left
        }
        else
        {
            neighbours[0] = temporaryID + 1; //assign "right" neighbour to one to right
        }
        if (temporaryID - range <= 0) //check if on bottom 
        {
            neighbours[1] = temporaryID + pow(range, 2) - range; //assign "bottom" neighbour to very top
        }
        else
        {
            neighbours[1] = temporaryID - range; //assign "bottom neighbour" to one below.
        }
        if ((temporaryID - 1) % range == 0) //check if on left
        {
            neighbours[2] = temporaryID + range - 1; //assign "left" neighbour to one on right
        }
        else
        {
            neighbours[2] = temporaryID - 1; //assign "left neighbour to one on the left
        }
        if (temporaryID + range > pow(range, 2))//check if on top
        {
            neighbours[3] = temporaryID - pow(range, 2) + range; //assign "top" neighbour to one on bottom
        }
        else
        {
            neighbours[3] = temporaryID + range; //assign "top" neighbour to one above
        }
        for (int i = 0; i < 4; i++)
        {
            neighbours[i] = neighbours[i] - 1;
        }
    }
};

class simulation
{
private:
    double J = 0;
    double B = 0;
    double kbt = 0;
    double K = 0;
    double U = 0;
    double DUDK = 0;
    double E = 0;
    int range = 0;
    int Nsweeps = 0;
    int sweepsPerMag = 0;
    int particleN0 = 0;
    int dataPoints = 0;
    double meanMag = 0;
    double stdError = 0;
    int cutoffValue = 0;
    vector<double> energyHistory;
    vector<double> sweepsHistory;
    vector<double> magnetisationHistory;
    vector<particle> lattice;

public: 
    simulation() {};
    simulation(double Jin, double Bin, double kbtin, int rangein, int Nsweepsin, int sweepsPerMagin)
    {
        J = Jin;
        B = Bin;
        kbt = kbtin;
        range = rangein;
        Nsweeps = Nsweepsin;
        sweepsPerMag = sweepsPerMagin;
        particleN0 = pow(range, 2);
        K = J / kbt;
        U = 0;
        dataPoints = 1 + Nsweeps / sweepsPerMag;
        lattice.resize(particleN0);
        sweepsHistory.resize(dataPoints);
        magnetisationHistory.resize(dataPoints);
        energyHistory.resize(dataPoints);
        for (int i = 0; i < particleN0; i++)
        {
            lattice[i] = particle(i,range);
        }  //initialise lattice
        magnetisationHistory[0] = calculateMagnetisation();
        sweepsHistory[0] = 0;
        for (int i = 1; i < Nsweeps + 1; i++) //for each sweep
        {
            lattice = performSweep();  //do the sweep
            if (i % sweepsPerMag == 0) // calculate the magnetisation of the lattice
            {
                magnetisationHistory[i / sweepsPerMag] = calculateMagnetisation();
                sweepsHistory[i / sweepsPerMag] = i;
                energyHistory[i / sweepsPerMag] = calculateEnergy();
            }
        }
        findResults();
    }
    double calculateMagnetisation()
    {
        int counter = 0;
        for (long i = 0; i < particleN0; i++)
        {
            if (lattice[i].spin == 1)
            {
                counter++;
            }
        }
        double totalSpin = 2 * counter - particleN0;
        return totalSpin / particleN0;
    }
    double calculateEnergy()
    {
        //note energy should be = 0 if perfectly allinged, and 1 if not for each site
        double energy = 0;
        for (int i = 0; i < particleN0; i++)
        {
            particle& info = lattice[i];

            int count = 0;
            for (int j = 0; j < 4; j++)
            {
                if (lattice[info.neighbours[j]].spin == 1)
                {
                    count++;
                }
                else
                {
                    count--;
                }
            }
            //count is sum of adjacent spins
            //using 4 - here, so that if all count and spin line up the energy contribution is 0
            double energyContribution = 4 - (count * info.spin);
            energy += energyContribution;
        }
        return energy;
    }
    vector<particle> performSweep()
    {
        
        for (int i = 0; i < particleN0; i++)
        {
            particle& info = lattice[i];

            int count = 0;
            for (int j = 0; j < 4; j++)
            {
                if (lattice[info.neighbours[j]].spin == 1)
                {
                    count++;
                }
                else
                {
                    count--;
                }
            }
            //count is sum of adjacent spins
            double energyContribution = 1 - count * info.spin;
          
            double partitionContribution = exp((2 * (J * count + B) * (-1 * info.spin)) / kbt);
            if (partitionContribution > 1)
            {
                info.spin = -info.spin;
            }
            else
            {
                double r = ((double)rand() / (RAND_MAX));
                if (partitionContribution > r)
                {
                    info.spin = -info.spin;
                }
            }
        }
        
        return lattice;
    }
    void printHistory()
    {
        for (int i = 0; i < dataPoints; i++)
        {
            cout << "\n Sweep Number : ";
            cout << sweepsHistory[i];
            cout << "\n Magnetisation : ";
            cout << magnetisationHistory[i];
        }
    }
    void findCutoff()
    {
        int minimumCutoff = Nsweeps - 100;
  
        vector<double> mean(dataPoints);
        vector<double> stdDev(dataPoints);
        mean[0] = magnetisationHistory[0];
        stdDev[0] = 0;
        for (int i = 30; i < dataPoints; i++) //for every data point after the 30th one
        {
            mean[i] = findMean(magnetisationHistory, i,true); //up to and incuding the index i
            stdDev[i] = findStdDeviation(magnetisationHistory, i,true);
            if (i >= minimumCutoff && stdDev[i] < stdDev[i - 1] && stdDev[i] < stdDev[i - 2] && stdDev[i] < stdDev[i - 3] &&cutoffValue == 0)
            {
                cutoffValue = i + 1;
            }
        }
    }
    void findResults()
    {
        //findCutoff();  //cutoff is inclusive, eg 17 would use sweep 17 to the end, so for 20 sweeps (21 data points) cutoff 17 would use index 17 to 20 giving 4 data points
        cutoffValue = 100;
        meanMag = findMean(magnetisationHistory, cutoffValue, false);
        double stdDev = findStdDeviation(magnetisationHistory, cutoffValue, false);
        double usableData = dataPoints - cutoffValue;
        stdError = stdDev / (usableData);
        findU();
    }

    void findU()
    {   
        vector<double> m2 (dataPoints);
        vector<double> m4 (dataPoints);
        vector<double> m2E(dataPoints);
        vector<double> m4E(dataPoints);


        for ( int i = 0; i < dataPoints; i++)
        {
            m2[i] = magnetisationHistory[i] * magnetisationHistory[i];
            m4[i] = pow(magnetisationHistory[i], 4);
            m2E[i] = m2[i] * energyHistory[i];
            m4E[i] = m4[i] * energyHistory[i];
        }
        double meanM2 = findMean(m2, dataPoints -1, true);
        double meanM4 = findMean(m4, dataPoints - 1, true);
        double meanM2E = findMean(m2E, dataPoints - 1, true);
        double meanM4E = findMean(m4E, dataPoints - 1, true);
        double meanE = findMean(energyHistory, dataPoints - 1, true);
        U = 1 - ( (static_cast<double>(1)/3) * (meanM4 / pow(meanM2,2) ) );
        DUDK = (static_cast<double>(1) - U) * (meanE - (2 * (meanM2E / meanM2)) + (meanM4E / meanM4));
    }

    double getU()
    {
        return U;
    }

    void printResults()
    {
        cout << "\n Cutoff Value = ";
        cout << cutoffValue;
        cout << "\n final mean = ";
        cout << meanMag;
        cout << "\n final std error = ";
        cout << stdError;
        cout << "\n U = ";
        cout << U;
        cout << "\n DUDK = ";
        cout << DUDK;
    }
    double getMeanMag()
    {
        return meanMag;
    }
    double getStdError()
    {
        return stdError;
    }
    void save(string filename)
    {
        vector<double> errors(dataPoints);
        for (int i = 0; i < dataPoints; i++)
        {
            errors[i] = 0;
        }
        saveToFile(filename, sweepsHistory, energyHistory,errors, dataPoints,"Number of sweeps","Total Energy");
    }
};

class batchSimulation
{
private:
    double meanMag = 0;
    double stdMagError = 0;
    double meanU = 0;
    double stdUError = 0;
public :
    batchSimulation() {};
    batchSimulation(int Nsims, double J, double B, double kbt, int range, int Nsweeps, int sweepsPerMag)
    {
        bool valid = false;  
        while (!valid)
        {
            vector<double> magnetisations(Nsims);
            vector<double> Uvalues(Nsims);
            vector<simulation> sims(Nsims);
            
            for (int i = 0; i < Nsims; i++)
            {
                sims[i] = simulation(J, B, kbt, range, Nsweeps, sweepsPerMag);
                magnetisations[i] = abs(sims[i].getMeanMag());
                Uvalues[i] = sims[i].getU();
            }

            meanMag = findMean(magnetisations, Nsims - 1, true);
            stdMagError = findStdError(magnetisations, Nsims - 1, true);
            meanU = findMean(Uvalues, Nsims - 1, true);
            stdUError = findStdError(Uvalues, Nsims - 1, true);
            //if (stdMagError < 0.05)  //note 0.05 gives close values,  Ad if to ensure data points are close
            //{
            valid = true;
            //}
        }
       
    }
    double getMeanMag()
    {
        return meanMag;
    }
    double getStdMagError()
    {
        return stdMagError;
    }
    double getStdUError()
    {
        return stdUError;
    }
    double getMeanU()
    {
        return meanU;
    }
    
};

int main()
{
    auto start = high_resolution_clock::now();
    srand((unsigned int)time(NULL));
   
    /*
    simulation sim1 = simulation(1,0,2.26,40,2000,10);
    sim1.printResults();
    sim1.save("criticalExponent.txt");
    */


    
    int Nbatches = 50;
    int simsPerBatch = 5;
    int L = 50;
    vector<batchSimulation> allBatches(Nbatches);
    vector<double> magnetisations(Nbatches);
    vector<double> stdMagErrors(Nbatches);
    vector<double> UValues(Nbatches);
    vector<double> KbtValues(Nbatches);
    vector<double> stdUErrors(Nbatches);
    for (int i = 0; i < Nbatches; i++)
    {
        KbtValues[i] = 2 + 0.008 * i;
        allBatches[i] = batchSimulation(simsPerBatch,1, 0, KbtValues[i], L, 3000, 20);
        magnetisations[i] = allBatches[i].getMeanMag();
        stdMagErrors[i] = allBatches[i].getStdMagError();
        UValues[i] = allBatches[i].getMeanU();
        stdUErrors[i] = allBatches[i].getStdUError();
    }
    saveToFile("VaryK_L=50.txt", KbtValues, UValues,stdUErrors, Nbatches, "Kbt", "U(Kbt,L)");
    
    //produceGraphs("batchsim3.txt");
    

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "\nTime taken by program is : ";
    cout << duration.count();
    cout << " microsec " << endl;

    cout << "\n";
   
    system("pause");
}
