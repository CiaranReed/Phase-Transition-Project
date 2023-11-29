
#include <iostream>
#include <map>
#include <fstream>
#include <vector>
#include <chrono>
using namespace std;
using namespace std::chrono;
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




double findMean(vector<double> arr, int cutoff, bool beforeCutoff) //cutoff inclusive, eg cutoff 3 means up to and including the element 3
{
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
    int range = 0;
    int Nsweeps = 0;
    int sweepsPerMag = 0;
    int particleN0 = 0;
    int dataPoints = 0;
    double meanMag = 0;
    double stdError = 0;
    vector<int> sweepsHistory;
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
        dataPoints = 1 + Nsweeps / sweepsPerMag;
        lattice.resize(particleN0);
        sweepsHistory.resize(dataPoints);
        magnetisationHistory.resize(dataPoints);
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
    void saveToFile(string filename)
    {
        ofstream file(filename);
        file << dataPoints;
        for (int i = 0; i < dataPoints; i++)
        {
            file << "\n";
            file << sweepsHistory[i];
            file << ",";
            file << magnetisationHistory[i];
        }
        file.close();
    }
    int findCutoff()
    {
        vector<double> mean(dataPoints);
        vector<double> stdDev(dataPoints);
        int cutoffValue = 0;
        mean[0] = magnetisationHistory[0];
        stdDev[0] = 0;
        for (int i = 1; i < dataPoints; i++) //for every data point after the first one
        {
            mean[i] = findMean(magnetisationHistory, i,true); //up to and incuding the index i
            stdDev[i] = findStdDeviation(magnetisationHistory, i,true);
            if (stdDev[i] < stdDev[i - 1] && cutoffValue == 0)
            {
                cutoffValue = i + 5;
                return cutoffValue;
            }
        }
    }
    void findResults()
    {
        int cutoffValue = findCutoff();  //cutoff is inclusive, eg 17 would use sweep 17 to the end, so for 20 sweeps (21 data points) cutoff 17 would use index 17 to 20 giving 4 data points
        meanMag = findMean(magnetisationHistory, cutoffValue, false);
        double stdDev = findStdDeviation(magnetisationHistory, cutoffValue, false);
        double usableData = dataPoints - cutoffValue;
        stdError = stdDev / (usableData);
    }
    void printResults()
    {
        cout << "\n final mean = ";
        cout << meanMag;
        cout << "\n final std error = ";
        cout << stdError;
    }
    double getMeanMag()
    {
        return meanMag;
    }
    double getStdError()
    {
        return stdError;
    }
};

int main()
{
    auto start = high_resolution_clock::now();
    srand((unsigned int)time(NULL));
    //simulation ( J, B, Kbt, Range, Sweeps, Sweeps per magnetisaton calculation)
    int Nsims = 10;
    vector<simulation> allSims(Nsims);
    vector<double> magnetisations(Nsims);
    vector<double> stdErrors(Nsims);
    for (int i = 0; i < Nsims; i++)
    {
        allSims[i] = simulation(0.5, -0.05, 1, 100, 100, 1);
        allSims[i].printResults();
        magnetisations[i] = allSims[i].getMeanMag();
        stdErrors[i] = allSims[i].getStdError();
    }
    double totalMeanMag = findMean(magnetisations, Nsims -1 , true);
    double totalStdErorr = findStdDeviation(magnetisations, Nsims -1 , true);
    cout << "\n Total Mean Mag = ";
    cout << totalMeanMag;
    cout << "\n total Std Error = ";
    cout << totalStdErorr;

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "\nTime taken by program is : ";
    cout << duration.count();
    cout << " microsec " << endl;


    allSims[1].saveToFile("magnetisation data.txt");

    string filename = "ProduceGraphs.py";
    string command = "Python ";
    command += filename;
    system(command.c_str());

    cout << "\n";


    system("pause");
  
}
//cutoff = 4
//mean = -0.34368
//std dev = 