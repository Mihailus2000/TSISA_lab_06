#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <tuple>
#include <iomanip>
#include <fstream>


std::ofstream noise("noisedFunction.txt");
std::ofstream filter("filteredFunction.txt");

int cnt = 0;

double DetermineFunction(double x){
    return sin(x) + 0.5;
}

double Omega(int K, std::vector<double> noisedFunction){
    double  res = 0.;
    for(int k = 1; k < K; k++){
        res += std::abs(noisedFunction.at(k) - noisedFunction.at(k - 1));
    }
    return res;
}

double Delta(int K, std::vector<double> filteredFunction, std::vector<double> noisedFunction){
    double res = 0.;
    for(int k = 0; k < K; k++){
        res += std::abs(filteredFunction.at(k) - noisedFunction.at(k));
    }
    res /= K;
    return res;
}

int DetermineAmountOfExperiments(double P, double eps){
    return static_cast<int>(round(log(1-P)/log(1 - eps/(M_PI))));
}

std::vector<double> GenerateAlpha(int r){
    std::mt19937 engine(std::random_device{}());
    if(r == 3){
        std::vector<double> alpha(3);
        auto generator = std::uniform_real_distribution<double>(0, 1);
        alpha.at(1) = generator(engine);
        alpha.at(0) = alpha.at(2) = 0.5 * (1 - alpha.at(1));
        return alpha;
    }
    if(r ==5){
        std::vector<double> alpha(5);
        auto generator = std::uniform_real_distribution<double>(0, 1);
        alpha.at(2) = generator(engine);

        auto generator2 = std::uniform_real_distribution<double>(0, 1 - alpha.at(2));
        alpha.at(1) = alpha.at(3) = 0.5 * generator2(engine);

        double sum = 0;
        for (int i = 1; i < 4; i++)
        {
            sum += alpha.at(i);
        }
        alpha.at(0) = alpha.at(4) = 0.5 * (1 - sum);
        return alpha;
    }
}

std::vector<double> Filtration(int K, int M, std::vector<double> noisedFunction, std::vector<double> alpha, std::vector<double> sequence) {
    std::vector<double> filteredFunction;
    if(M == 1)
        filteredFunction.push_back(noisedFunction.at(0));
    if(M == 2){
        filteredFunction.push_back(noisedFunction.at(0));
        filteredFunction.push_back(noisedFunction.at(1));
        filteredFunction.push_back(noisedFunction.at(2));
    }
    double res = 1.;
    for (int i = M; i <= K-M; i++){
        for (int j = i-M; j < i+M; j++) {
            res *= std::pow(noisedFunction.at(j), alpha.at(j + M - i));
        }
        filteredFunction.push_back(res);
        res = 1.;
    }
//    if(cnt == 0) {
//        for (int i = 0; i < filteredFunction.size(); i++) {
//            filter << '(' << sequence.at(i) << ";" << filteredFunction.at(i) << ')' << std::endl;
//        }
//        filter << "-----------------------------\n";
//        cnt++;
//    }
    return filteredFunction;
}

std::vector<double> GenerateNoisedFunction(std::vector<double> sequence){
    std::vector<double> noisedFunction;
//    noisedFunction.resize(originalFunction.size());
    std::mt19937 engine(std::random_device{}());
    auto generator = std::uniform_real_distribution<double>(-0.25, 0.25);
    for(int k = 0; k < sequence.size(); k++){
        double alpha = generator(engine);
        noisedFunction.push_back(DetermineFunction(sequence.at(k)) + alpha);
    }
    for (int i = 0; i < noisedFunction.size(); i++)
    {
        noise << '(' << sequence.at(i) << ";" << noisedFunction.at(i) << ')' << std::endl;
    }
    noise <<  "-----------------------------\n";
    return noisedFunction;
}

std::vector<double> DiscreteSequence(double xMin, double xMax, int K){
    std::vector<double> sequence;
    sequence.resize(K);
    for(int k = 0; k < K; k++){
            sequence.at(k) =   xMin + k * (xMax - xMin) / K;
    }

    return sequence;
}

double Distance(double omega, double delta){
    return std::abs(omega) + std::abs(delta);
}

std::tuple<std::vector<double >,double,double, double> FindOptimalForLambda(int n, int lambda, int K, int r, std::vector<double> sequence, std::vector<double> noisedFunction){
    double minimum = 1000.;
    std::tuple<std::vector<double >,double,double, double> optimal; // {vectorAlpha, J, omega, delta}
    std::vector<double> filteredFunction;
    for(int i = 0; i < n; i++){
        std::vector<double> alphaVector = GenerateAlpha(r);
        filteredFunction = Filtration(K, std::floor((r - 1) / 2), noisedFunction, alphaVector, sequence);
        double omega = Omega(K,filteredFunction);
        double delta = Delta(K,filteredFunction,noisedFunction);

        double J = lambda * omega + (1-lambda) * delta;

        if(J < minimum){
            minimum = J;
            optimal = {alphaVector, J, omega, delta};
        }

        cnt++;
    }
    if(cnt == 0) {
        for (int i = 0; i < filteredFunction.size(); i++) {
            filter << '(' << sequence.at(i) << ";" << filteredFunction.at(i) << ')' << std::endl;
        }
        filter << "-----------------------------\n";
        cnt++;
    }

    return optimal;
}

double FindOptimalLambda(int L, int K, int r, std::vector<double> sequence){
    double minimumJ = 100000000.;
    double optimalLambda = -1.;

    auto noisedFunction = GenerateNoisedFunction(sequence);
    double dist;
    for(int ind = 0 ; ind <= L; ind++){
        double lambda = static_cast<double>(ind * 1.0) / L;
        std::vector<double> alpha;
        double omega, delta, j;
        std::tie(alpha,j,omega,delta) = FindOptimalForLambda(DetermineAmountOfExperiments(0.95, 0.01), lambda, K, r, sequence, noisedFunction);
        dist = Distance(omega,delta);
        if (dist < minimumJ) {
            minimumJ = dist;
            optimalLambda = lambda;
        }

        std::cout << std::fixed << "| " << std::setprecision(4) << lambda << " | " << j << " |";
        std::cout << " [";
        for(auto alpEl = alpha.begin(); alpEl != alpha.end(); alpEl++){
            std::cout << std::fixed << std::setprecision(4) << *alpEl;
            if(alpEl != alpha.end()-1)
                std::cout << ", ";
        }
        std::cout << "] | ";

        std::cout << std::fixed << std::setprecision(4) << omega << " | " << delta << " |\n";
    }
    std::cout << "\n\n";
    return optimalLambda;
}


int main() {
    std::cout << "+--------+--------+--------------------------+---------+--------+\n"
                 "| lambda |    J   |            alpha         |    w    |    d   |\n"
                 "+--------+--------+--------------------------+---------+--------+\n";
    auto optimalLambda = FindOptimalLambda(10,10,3,DiscreteSequence(0.,M_PI, 100));
    cnt = 0;
    std::cout << "Optimal lambda for 'r' = 3 : " << std::fixed << std::setprecision(4) << optimalLambda;
    std::cout << "\n\n";
    std::cout << "+--------+--------+------------------------------------------+---------+--------+\n"
                 "| lambda |    J   |                  alpha                   |    w    |    d   |\n"
                 "+--------+--------+------------------------------------------+---------+--------+\n";
    auto optimalLambda2 = FindOptimalLambda(10,10,5,DiscreteSequence(0.,M_PI, 100));
    std::cout << "Optimal lambda for 'r' = 5 : " << std::fixed << std::setprecision(4) << optimalLambda2;
}