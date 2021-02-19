#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>

//Read in Output values
//compare to "real" values
//return number

double calculateLikelihood(std::vector<double> solution, int rank){
    int numProbes = 1;
    double likelihood =  0.0;  //0.25*(ll[0] + ll[1] + ll[2] + ll[3]);

    //calculate differences at the 4 different probe points
    for(int i_probe=1; i_probe <= numProbes; i_probe++){
        std::cout << "Difference in arrival time " << solution[0]-1.85232 << std::endl;
        std::cout << "Difference in max height " << solution[1]-26.23200000001134  << std::endl;
	likelihood -= 1.0/numProbes * .5 * (std::pow(solution[0]-26.23200000001134,2)/1.0+std::pow(solution[1]-1.85232,2)/0.01);
    }
    std::ofstream ost;
    ost.open("likelihood_r"+std::to_string(rank)+".log", std::ios::app);
    ost << std::exp(likelihood) << std::endl;
    return likelihood;
}

double calculateLikelihood_endtime(std::vector<double> solution, int rank){
    int numProbes = 4;
    double likelihood =  0.0;  //0.25*(ll[0] + ll[1] + ll[2] + ll[3]);

    //calculate differences at the 4 different probe points
    for(int i_probe=1; i_probe <= numProbes; i_probe++){

        double diff;
        double dh;

        {
        std::string filestr ="Output/Reference/waveheight"+std::to_string(i_probe)+"-rank-0.probe";
        const char* filename = filestr.c_str();
        std::ifstream infile(filename);
        std::string line;
        std::getline(infile, line);
        while (std::getline(infile, line))
        {
                std::stringstream ss(line);
                double d1; ss >> d1;   // read 1
                char cc; ss >> cc;   // read 1
                double d2; ss >> d2;   // read 2
                char cc2; ss >> cc2;   // read 1
                double d3; ss >> d3;   // read 2
                char cc3; ss >> cc3;   // read 1
                double d4; ss >> d4;   // read 2
                diff = d4; 
        }
        }
	std::cout << "Solution " << i_probe - 1 << " "  << solution[i_probe-1] << std::endl;
        if(std::isfinite(solution[i_probe-1]))
		diff-=solution[i_probe-1];
	else
	   diff-= 10000;

        std::cout << "Diff     " << i_probe - 1 << " " << diff << std::endl;
	likelihood -= 1.0/numProbes * .5 * diff*diff * 500;
    }
    std::ofstream ost;
    ost.open("likelihood_r"+std::to_string(rank)+".log", std::ios::app);
    ost << std::exp(likelihood) << std::endl;
    return likelihood;
}
