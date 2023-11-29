
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
};

int main()
{
    auto start = high_resolution_clock::now();
    srand((unsigned int)time(NULL));
   
    simulation sim1 = simulation(0.5, 0.05, 1, 32, 100, 1);

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "Time taken by program is : ";
    cout << duration.count();
    cout << " microsec " << endl;

    sim1.printHistory();
    sim1.saveToFile("magnetisation data.txt");

    string filename = "ProduceGraphs.py";
    string command = "Python ";
    command += filename;
    system(command.c_str());

    cout << "\n";


    system("pause");
  
}
