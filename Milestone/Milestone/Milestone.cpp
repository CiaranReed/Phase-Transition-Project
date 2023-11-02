
#include <iostream>
#include <map>
#include <fstream>
#include <vector>
using namespace std;

//constants

const int range = 32;
const int Nsweeps = 100;
const int sweepsPerMag = 1;   //make sure that Nsweeps / sweepsPerMag is an integer or code will break!!!
const double J = 0.5; //interaciton strength for neighbours (1 for milestone)
const double B = -0.05; //external mag field strength for whole lattice
const double kbt = 1;  // meV   26 is from the notes for room temp, lower means lower partition means less chance of random fluctuation
const int particleN0 = pow(range, 2);
const int dataPoints = 1 + Nsweeps / sweepsPerMag;



//if particle flips to be parallel with adjacent spins, then -J.S1.S2 is  negative as S1.S2 is positive, 
// this means E is more negative which is favorable
// positive J means particles tend to become parallel
// B decides which direction they tend to
// whichever system has lowest potential (energy) is what it tends to


class particle
{
private: 
    unsigned int ID;
    bool spin; //up == true, down == false
    int spinValue;

public:
    particle() {};
    particle(int identifier)
    {
        setSpin();
        setID(identifier);
        setNeighbours();
    }
    int neighbours[4];

    bool getSpin()
    {
        return spin;
    }

    int getID()
    {
        return ID;
    }

    int sumAdjacentSpins(vector<particle> lattice)
    {
        int count = 0;
        for (int i = 0; i < 4; i++)
        {
            if (lattice[neighbours[i]].getSpin())
            {
                count++;
            }
            else
            {
                count--;
            }
        }
        return count;
 
    }
    
    void setSpin()
    {
        spin = rand() % 2; //assign a random spin to the particle
        if (spin)
        {
            spinValue = 1;
        }
        else
        {
            spinValue = -1;
        }
    }

    void flipSpin()
    {
        if (spin)
        {
            spin = false;
            spinValue = -1;
        }
        else
        {
            spin = true;
            spinValue = 1;
        }
    }

    int getSpinValue()
    {
        return spinValue;
    }

    double calculateEnergyChange(vector<particle> lattice)  //before we flip, find energy change if we do flip
    {
        return -2 * (J * sumAdjacentSpins(lattice) + B) * (-1 * spinValue);
    }

    void setID(int identifier)
    {
        ID = identifier;
    }

    double findPartition( double E)
    {
        return exp((-1 * E) / kbt);
    }

    void sweepWithInfo(vector<particle> lattice)
    {
        cout << "\n For particle ";
        cout << ID;
        cout << "\n Spin = ";
        cout << getSpinValue();
        cout << "\n Sum of adjacent spins = ";
        cout << sumAdjacentSpins(lattice);
        double delE = calculateEnergyChange(lattice);
         cout << "\n DelE = ";
         cout << delE;
        double partition = findPartition(delE);
         cout << "\n Partition = ";
         cout << partition;
        if (partition > 1)
        {
            flipSpin();
            cout << "\n";
            cout << " as partition is > 1, its Had its spin flipped";
        }
        else
        {
            double r = ((double)rand() / (RAND_MAX));
            if (partition > r)
            {
                 flipSpin();
                 cout << "\n r = ";
                 cout << r;
                 cout << "\n as partition is > r, its Had its spin flipped";
            }
            else
            {               
                cout << "\n r = ";
                cout << r;
                cout << "\n as partition < r, it Did not have its spin flipped";
            }
        }
        cout << "\n \n";
    }

    void sweep(vector<particle> lattice)
    {
        double partition = findPartition( calculateEnergyChange(lattice));
        if (partition > 1)
        {
            flipSpin();
        }
        else
        {
            double r = ((double)rand() / (RAND_MAX));
            if (partition > r)
            {
                flipSpin();
            }        
        }
    }

    void setNeighbours()
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



double countTotalSpinUp(vector<particle> lattice)
{
    int counter = 0;
    for (long i = 0; i < particleN0; i++)
    {
        if (lattice[i].getSpin())
        {
            counter++;
        }
    }
    return counter;
}

double countTotalSpinDown(vector<particle> lattice)
{
    return (particleN0 - countTotalSpinUp(lattice));
}

void printTotalSpinUp(vector<particle> lattice)
{
    cout << "\n The total number of spin ups is : ";
    cout << countTotalSpinUp(lattice);
    cout << "\n";
}

void printTotalSpinDown(vector<particle> lattice)
{
    cout << "\n The total number of spin downs is : ";
    cout << countTotalSpinDown(lattice);
    cout << "\n";
}

double calculateMagnetisation(vector<particle> lattice)
{
    double totalSpin = countTotalSpinUp(lattice) - countTotalSpinDown(lattice);
    return totalSpin / particleN0;

}

vector<particle> performSweep(vector<particle> lattice)
{
    for (int i = 0; i < particleN0; i++)
    {
        particle& info = lattice[i];
        info.sweep(lattice);
    }
    return lattice;
}

void performSweepWithInfo(vector<particle> lattice)
{
    printTotalSpinUp(lattice);
    printTotalSpinDown(lattice);
    for (int i = 0; i < particleN0; i++)
    {
        lattice[i].sweepWithInfo(lattice);
    }
    cout << "\n After sweeping : ";
    printTotalSpinUp(lattice);
    printTotalSpinDown(lattice);
    cout << "\n Magnetisaiton = ";
    cout << calculateMagnetisation(lattice);
    cout << "\n--------------------------------------------------------------- \n ";
}

void printHistory(double mHistory[], int sHistory[], int nCalculations)
{
    for (int i = 0; i < nCalculations; i++)
    {
        cout << "\n Sweep Number : ";
        cout << sHistory[i];
        cout << "\n Magnetisation : ";
        cout << mHistory[i];
    }
}


int main()
{
    srand((unsigned int)time(NULL));
    vector<particle> lattice(particleN0);
    int sweepsHistory[dataPoints];
    double magnetisationHistory[dataPoints];

    for (int i = 0; i < particleN0; i++)
    {
        lattice[i] = particle(i);  
    }  //initialise lattice

    magnetisationHistory[0] = calculateMagnetisation(lattice);
    sweepsHistory[0] = 0;
    //set the first value of the data

    for (int i = 1; i < Nsweeps + 1; i++) //for each sweep
    {
       lattice = performSweep(lattice);
       if (i % sweepsPerMag == 0)
       {
           magnetisationHistory[i/sweepsPerMag] = calculateMagnetisation(lattice);
           sweepsHistory[i/sweepsPerMag] = i;
       }
      
    }
    printHistory(magnetisationHistory, sweepsHistory,dataPoints); //at the moment we are finding the magnetisation after every sweep
    ofstream file("magnetisation data.txt");
    file << dataPoints;
    for (int i = 0; i < dataPoints; i++)
    {
        file << "\n";
        file << sweepsHistory[i];
        file << ",";
        file << magnetisationHistory[i];
    }
    file.close();
    
    string filename = "ProduceGraphs.py";
    string command = "Python ";
    command += filename;
    system(command.c_str());

    cout << "\n";


    system("pause");
  
}
