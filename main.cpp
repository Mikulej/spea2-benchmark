#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <regex>
#include <sstream>
#include <cfloat>

# define PI 3.14159265358979323846 

struct Solution {
    size_t id;
    std::vector<double> values;
    std::vector<double> objectiveScores;//[0] = f1(values), [1] = f2(values), ...

    size_t strength; //S(k)
    size_t weakness; //R(i)
    double density; //D(i)

    // Solution operator=(const Solution& s) {
    //     id = s.id;
    //     values = s.values;
    //     objectiveScores = s.objectiveScores;
    //     return *this;
    // }

    bool operator==(const Solution& other) const {
        return id == other.id;
    }

    double getFitness() const {
        return (double)weakness + density;
    }

};

static size_t gradeAmount1 = 0;
static size_t gradeAmount2 = 0;
static size_t generateId = 0;

std::mt19937 gen(std::time({}));
std::uniform_real_distribution<> dis(0.0, 1.0);
std::uniform_real_distribution<> disZdt4(-5.0, 5.0);

double zdt1f1(const std::vector<double>& values, double parameter) {
    gradeAmount1++;
    return values[0];
}

double zdt1f2(const std::vector<double>& values, double f1score) {
    gradeAmount2++;
    double g = 1. + (9. * values[1]);
    double h = 1. - sqrt(f1score / g);
    return g * h;
}

double zdt2f1(const std::vector<double>& values, double parameter) {
    gradeAmount1++;
    return values[0];
}

double zdt2f2(const std::vector<double>& values, double f1score) {
    gradeAmount2++;
    double g = 1. + (9. * values[1]);
    double h = 1. - ((f1score / g)*(f1score / g));
    return g * h;
}

double zdt3f1(const std::vector<double>& values, double parameter) {
    gradeAmount1++;
    return values[0];
}

double zdt3f2(const std::vector<double>& values, double f1score) {
    gradeAmount2++;
    double g = 1. + (9. * values[1]);
    double h = 1. - sqrt(f1score/g) - ((f1score/g)*sin(10.*PI*f1score));
    return g * h;
}

double zdt4f1(const std::vector<double>& values, double parameter) {
    gradeAmount1++;
    return values[0];
}

double zdt4f2(const std::vector<double>& values, double f1score) {
    gradeAmount2++;
    double g = 11. + ((values[1]*values[1])-(10.*cos(4.*PI*values[1])));
    double h = 1. - sqrt(f1score / g);
    return g * h;
}

double zdt6f1(const std::vector<double>& values, double parameter) {
    gradeAmount1++;
    return 1. - exp(-4.*values[0])*pow(sin(6.*PI*values[0]),6);
}

double zdt6f2(const std::vector<double>& values, double f1score) {
    gradeAmount2++;
    double g = 1. + (9. * pow(values[1],0.25));
    double h = 1. - ((f1score / g)*(f1score / g));
    return g * h;
}




std::vector<Solution> generateRandom(int populationSize, int dimensions, std::vector<double (*)(const std::vector<double>& values, double parameter)>& objectives) {
    std::vector<Solution> population;

    if(objectives[0] == &zdt4f1 && objectives[1] == &zdt4f2){//ZDT 4
        for (int i = 0; i < populationSize; i++) {
            Solution sol;
            sol.id = generateId++;
            sol.values.push_back(dis(gen));
            for (int j = 1; j < dimensions; j++) {
                sol.values.push_back(disZdt4(gen));
            }

            sol.objectiveScores.push_back(objectives[0](sol.values, 0.));
            sol.objectiveScores.push_back(objectives[1](sol.values, sol.objectiveScores[0]));

            population.push_back(sol);
        }
    }
    else{ //ZDT 1,2,3,6
        for (int i = 0; i < populationSize; i++) {
            Solution sol;
            sol.id = generateId++;
            for (int j = 0; j < dimensions; j++) {
                sol.values.push_back(dis(gen));
            }

            sol.objectiveScores.push_back(objectives[0](sol.values, 0.));
            sol.objectiveScores.push_back(objectives[1](sol.values, sol.objectiveScores[0]));

            population.push_back(sol);
        }
    }
    
    return population;
}

//Checks if a dominates b
bool dominates(const Solution& a, const Solution& b, std::vector<double (*)(const std::vector<double>& values, double parameter)>& objectives) {
    bool atLeastOneBetter = false;
    double aScore = 0;
    double bScore = 0;
    for (size_t i = 0; i < objectives.size(); ++i) {
        aScore = a.objectiveScores[i];
        bScore = b.objectiveScores[i];

        if (aScore > bScore) {
            return false;
        }
        if (aScore < bScore) {
            atLeastOneBetter = true;
        }
    }
    return atLeastOneBetter;
}

//Kung algorithm

std::vector<Solution> kungRecursive(std::vector<Solution>& P, std::vector<double (*)(const std::vector<double>& values, double parameter)>& objectives) {
    if (P.size() <= 1) {
        return P;
    }

    //Divide population into Top and Bottom
    size_t mid = P.size() / 2;
    std::vector<Solution> T(P.begin(), P.begin() + mid);
    std::vector<Solution> B(P.begin() + mid, P.end());

    //Recursion
    std::vector<Solution> topFront = kungRecursive(T, objectives);
    std::vector<Solution> bottomFront = kungRecursive(B, objectives);

    std::vector<Solution> merged = topFront;

    //Check if topfront solutions dominates bottomFront
    for (const auto& b : bottomFront) {
        bool isDominatedByTop = false;
        for (const auto& t : topFront) {
            if (dominates(t, b, objectives)) {
                isDominatedByTop = true;
                break;
            }
        }
        if (!isDominatedByTop) {
            merged.push_back(b);
        }
    }

    return merged;
}

std::vector<Solution> kungPareto(std::vector<Solution> population, std::vector<double (*)(const std::vector<double>& values, double parameter)>& objectives) {
    //Sort the population based on the '[](const Solution& a, const Solution& b' values
    std::sort(population.begin(), population.end(), [&objectives](const Solution& a, const Solution& b) {
        double aScore = 0;
        double bScore = 0;
        for (size_t i = 0; i < objectives.size(); ++i) {
            aScore = a.objectiveScores[i];
            bScore = b.objectiveScores[i];

            //Check if aScore and bScore are the same
            if (abs(aScore - bScore) < 0.001) {
                //Move to another objective function
                continue;
            }
            else {
                return aScore < bScore;
            }
        }
        return aScore < bScore; // use last objective function
        });

    return kungRecursive(population, objectives);
}

std::vector<Solution> getDominated(const std::vector<Solution>& population, const std::vector<Solution>& nonDominated) {
    std::vector<Solution> tempPopulation = population;
    std::vector<Solution> tempNonDominated = nonDominated;
    std::vector<Solution> dominated;

    std::sort(tempPopulation.begin(), tempPopulation.end(), [](const Solution& a, const Solution& b) {
        return a.id < b.id;
        });
    std::sort(tempNonDominated.begin(), tempNonDominated.end(), [](const Solution& a, const Solution& b) {
        return a.id < b.id;
        });

    int j = 0;
    for (int i = 0; i < tempPopulation.size(); i++) {
        if ((j < tempNonDominated.size()) && (tempNonDominated[j].id == tempPopulation[i].id)) {
            //is non-domianted,skip
            j++;
        }
        else {
            dominated.push_back(population[i]);
        }
    }
    return dominated;
}

double getDistance(const Solution& a, const Solution& b) {
    return sqrt(((a.objectiveScores[0] - b.objectiveScores[0]) * (a.objectiveScores[0] - b.objectiveScores[0])) + ((a.objectiveScores[1] - b.objectiveScores[1]) * (a.objectiveScores[1] - b.objectiveScores[1])));
}

size_t archiveTruncationProcedure(const std::vector<Solution>& archive, std::vector<double (*)(const std::vector<double>& values, double parameter)>& objectives) {
    std::vector<std::vector<Solution>> sortedArchive;
    for (size_t i = 0; i < archive.size(); i++) {
        std::vector<Solution> archiveNoI = archive;
        archiveNoI.erase(archiveNoI.begin() + i);
        sortedArchive.push_back(archiveNoI);
    }

    for (size_t i = 0; i < archive.size(); i++) {
        Solution origin = archive[i];

        std::sort(sortedArchive[i].begin(), sortedArchive[i].end(), [&origin, &objectives](const Solution& a, const Solution& b) {
            return getDistance(origin, a) < getDistance(origin, b);
            });

        //now, sortedArchive[i] is sorted by distance to origin Solution
    }

    //now, each vector in sortedArchive is sorted by distance to i'th Solution, closeset one

    // for(size_t i = 0; i < archive.size(); i++){
    //     std::cout << "i=" << i << std::endl;
    //     Solution origin = archive[i];
    //     for(size_t j = 0; j < sortedArchive[i].size(); j++){
    //         std::cout << getDistance(origin,sortedArchive[i][j]) << " ";
    //     }
    //     std::cout << std::endl;
    // }

    std::vector<size_t> skipIds;
    while (true) {
        //check which solution has the smallest distance
        std::vector<size_t> minIds;
        double min = DBL_MAX;
        for (size_t i = 0; i < archive.size(); i++) {

            if (find(skipIds.begin(), skipIds.end(), i) != skipIds.end()) {//found i in skipIds
                continue;
            }

            Solution origin = archive[i];
            double distance = getDistance(origin, sortedArchive[i][0]);
            if (min > distance) {
                min = distance;
                minIds.clear();
                minIds.push_back(i);
            }
            else if (min == distance) {
                minIds.push_back(i);
            }

        }

        if (minIds.size() == 1) {
            //one smallest solution, we delete it
            return minIds[0];
        }
        else {
            //there is more than one "smallest" solutions
            if (sortedArchive[0].size() > 1) {
                //for each archive solution that is not in minIds, add it to skipIds
                for (size_t i = 0; i < archive.size(); i++) {
                    if (find(minIds.begin(), minIds.end(), i) == minIds.end()) {// could not find i in minIds, add i to skipIds
                        skipIds.push_back(i);
                    }
                    //delete the nearest ones(first element) and do the algo again
                    sortedArchive[i].erase(sortedArchive[i].begin());
                }
                continue;
            }
            else { //compared all possible distances, solution are identical
                return minIds[0];
            }
        }
    }
}

void Evaluate(std::vector<Solution>& solutions, std::vector<Solution>& populationPlusArchive, std::vector<double (*)(const std::vector<double>& values, double parameter)>& objectives) {
    //Get S(k) for each k in populationPlusArchive
    for (size_t k = 0; k < populationPlusArchive.size(); k++) {
        populationPlusArchive[k].strength = 0;
        for (size_t j = 0; j < populationPlusArchive.size(); j++) {
            if (k == j) { continue; }
            if (dominates(populationPlusArchive[k], populationPlusArchive[j], objectives)) {
                populationPlusArchive[k].strength++;
            }
        }
    }

    //Get R(i) for each i in solutions(Population/Archive)
    for (size_t i = 0; i < solutions.size(); i++) {
        solutions[i].weakness = 0;
        for (size_t j = 0; j < populationPlusArchive.size(); j++) {
            if (dominates(populationPlusArchive[j], solutions[i], objectives)) {
                solutions[i].weakness += populationPlusArchive[j].strength;
            }
        }
    }


    //Get D(i) for each i in solutions(Population/Archive) 
    const size_t k = sqrt(populationPlusArchive.size());
    for (size_t i = 0; i < solutions.size(); i++) {
        int removeId = solutions[i].id;

        //WARNING: Assuming that neighbours are from "populationPlusArchive"
        std::vector<Solution> neighbours = populationPlusArchive;

        //Remove solutions[i] from neighbours based on solutions[i].id
        auto iter = std::find_if(neighbours.begin(), neighbours.end(), [removeId](const Solution& s) {
            return s.id == removeId;
            });
        if (iter != neighbours.end()) {
            neighbours.erase(iter);
        }
        // else{//If not found it is possible that lastArchive doesnt have any id from current populationPlusArchive, it's okay
        //     std::cout <<"Something went wrong at Evaluate"<<std::endl;
        // }


        //Sort neighbours
        Solution& origin = solutions[i];
        std::sort(neighbours.begin(), neighbours.end(), [&origin](const Solution& a, const Solution& b) {
            return getDistance(origin, a) < getDistance(origin, b);
            });

        double sigma = getDistance(origin, neighbours[k]);
        solutions[i].density = 1.0 / (sigma + 2.0);

    }
}

std::vector<Solution> tournamentSelection(const std::vector<Solution>& solutions, int winnersAmount) {
    std::vector<Solution> winners = solutions;
    std::vector<Solution> brackets;

    while (winners.size() > winnersAmount) {
        brackets = winners;

        std::shuffle(brackets.begin(), brackets.end(), gen);
        winners.clear();

        while (brackets.size() >= 2) {
            int aIndex = std::uniform_int_distribution<>(0, brackets.size() - 1)(gen);
            Solution a = brackets[aIndex];
            brackets.erase(brackets.begin() + aIndex);
            int bIndex = std::uniform_int_distribution<>(0, brackets.size() - 1)(gen);
            Solution b = brackets[bIndex];
            brackets.erase(brackets.begin() + bIndex);

            //max(a.F(i),b.F(i)) 
            if (a.getFitness() < b.getFitness()) {
                winners.push_back(a);
            }
            else {
                winners.push_back(b);
            }

        }
        if (brackets.size() > 0) {
            winners.push_back(brackets[0]);
        }

    }

    return winners;
}

std::vector<Solution> recombine(const std::vector<Solution>& populationMating, int offspringAmount, std::vector<double (*)(const std::vector<double>& values, double parameter)>& objectives) {
    std::vector<Solution> offspring;
    int dimensions = populationMating[0].values.size();
    while (offspring.size() < offspringAmount) {
        std::vector<Solution> possibleMates = populationMating;
        int aIndex = std::uniform_int_distribution<>(0, possibleMates.size() - 1)(gen);
        Solution a = possibleMates[aIndex];
        possibleMates.erase(possibleMates.begin() + aIndex);
        int bIndex = std::uniform_int_distribution<>(0, possibleMates.size() - 1)(gen);
        Solution b = possibleMates[bIndex];
        possibleMates.erase(possibleMates.begin() + bIndex);

        //Two point crossover
        int point1 = std::uniform_int_distribution<>(0, dimensions - 1)(gen);
        int point2 = std::uniform_int_distribution<>(0, dimensions - 1)(gen);

        Solution child;
        child.id = generateId++;

        if (point1 < point2) {
            for (int i = 0; i < dimensions; i++) {
                if (point1 <= i && i < point2) {
                    child.values.push_back(a.values[i]);
                }
                else {
                    child.values.push_back(b.values[i]);
                }
            }
        }
        else if (point1 > point2) {
            for (int i = 0; i < dimensions; i++) {
                if (point2 <= i && i < point1) {
                    child.values.push_back(b.values[i]);
                }
                else {
                    child.values.push_back(a.values[i]);
                }
            }
        }
        else {//point1 == point2, One point crossover
            int whichSolutionFirst = std::uniform_int_distribution<>(0, 1)(gen);
            if (whichSolutionFirst) {
                //guarantee 1 value from a
                child.values.push_back(a.values[0]);
                for (int i = 1; i < dimensions; i++) {
                    if (i < point1) {
                        child.values.push_back(a.values[i]);
                    }
                    else {
                        child.values.push_back(b.values[i]);
                    }
                }
            }
            else {
                //guarantee 1 value from b
                child.values.push_back(b.values[0]);
                for (int i = 1; i < dimensions; i++) {
                    if (i < point1) {
                        child.values.push_back(b.values[i]);
                    }
                    else {
                        child.values.push_back(a.values[i]);
                    }
                }
            }

        }

        child.objectiveScores.push_back(objectives[0](child.values, 0.));
        child.objectiveScores.push_back(objectives[1](child.values, child.objectiveScores[0]));
        offspring.push_back(child);

    }
    return offspring;
}

void mutate(std::vector<Solution>& solutions, int mutationAmount) {
    std::normal_distribution<> randomOffset(0., 0.3);
    // std::uniform_real_distribution randomOffset(-0.3,0.3);
    int dimensions = solutions[0].values.size();
    for (Solution& s : solutions) {
        std::vector<double> nonMutated = s.values;
        for (int i = 0; i < mutationAmount; i++) {
            int mutateIndex = std::uniform_int_distribution<>(0, dimensions - 1)(gen);
            double original = nonMutated[mutateIndex];
            nonMutated.erase(nonMutated.begin() + mutateIndex);

            double mutated = std::clamp<double>(original + randomOffset(gen), 0., 1.);
            s.values[mutateIndex] = mutated;
        }
    }
}

std::vector<Solution> Spea2(const std::vector<Solution>& startPopulation, std::vector<double (*)(const std::vector<double>& values, double parameter)>& objectives) {
    constexpr size_t budget = 20000;
    std::vector<Solution> population = startPopulation;
    std::vector<Solution> archive;
    std::vector<Solution> lastArchive = generateRandom(population.size(), population[0].values.size(),objectives);
    //std::vector<Solution> populationNonDominated = kungPareto(population, objectives);


    while ((gradeAmount1 < budget) && (gradeAmount2 < budget)) {
        std::vector<Solution> populationNonDominated = kungPareto(population, objectives);

        //Copy non-dominated memebers of Population to Archive
        for (const Solution& s : populationNonDominated) {
            archive.push_back(s);
        }

        //Remove any solution duplicates
        std::vector<Solution> archiveUnique;

        for (size_t i = 0; i < archive.size(); i++) {
            bool unique = true;
            for (size_t j = i + 1; j < archive.size(); j++) {
                //Compre solution by the objective to remove duplicates
                if (archive[i].id == archive[j].id) {
                    unique = false;
                    break;
                }
            }

            if (unique) {
                archiveUnique.push_back(archive[i]);
            }
        }

        //Keep only non-dominated Solutions in Archive
        archive = kungPareto(archiveUnique, objectives);

        //Create Population + Archive vector
        std::vector<Solution> populationPlusArchive;
        populationPlusArchive.reserve(population.size() + archive.size()); // preallocate memory
        populationPlusArchive.insert(populationPlusArchive.end(), population.begin(), population.end());
        populationPlusArchive.insert(populationPlusArchive.end(), archive.begin(), archive.end());

        //Evaluate S(k), R(i), D(i), F(i) for Population and for Archive
        Evaluate(population, populationPlusArchive, objectives);
        Evaluate(archive, populationPlusArchive, objectives);
        Evaluate(lastArchive, populationPlusArchive, objectives);

        //Sort lastArchive based on F(i)
        std::sort(lastArchive.begin(), lastArchive.end(), [](const Solution& a, const Solution& b) {
            return a.getFitness() < b.getFitness();
            });

        //Make sure Archive has exactly the same size as Population 
        size_t fillUpArchiveIndex = 0;
        while (archive.size() < population.size()) //fillup archive using lastArchive
        {
            if (fillUpArchiveIndex == lastArchive.size()) {
                std::cout << "ERROR: lastArchive has insufficent amount of Solutions to fill up archive" << std::endl;
                return std::vector<Solution>();
            }

            Solution candidate = lastArchive[fillUpArchiveIndex];
            fillUpArchiveIndex++;

            bool alreadyExists = false;
            for (const auto& existing : archive) {
                if(existing == candidate){
                    alreadyExists = true;
                    break;
                }

                // if (std::abs(existing.objectiveScores[0] - candidate.objectiveScores[0]) < 1e-9 &&
                //     std::abs(existing.objectiveScores[1] - candidate.objectiveScores[1]) < 1e-9) {
                //     alreadyExists = true;
                //     break;
                // }
            }

            if (!alreadyExists) {
                archive.push_back(candidate);
            }

        }
        while (archive.size() > population.size())//archive truncation prodecure (remove smallest distances)
        {
            size_t index = archiveTruncationProcedure(archive, objectives);
            archive.erase(archive.begin() + index);
        }

        //Update Population + Archive vector, after changing archive (no need to evaluate, it was already done - results are copied)
        populationPlusArchive.clear();
        populationPlusArchive.reserve(population.size() + archive.size()); // preallocate memory
        populationPlusArchive.insert(populationPlusArchive.end(), population.begin(), population.end());
        populationPlusArchive.insert(populationPlusArchive.end(), archive.begin(), archive.end());

        //Tournament
        std::vector<Solution> populationMating = tournamentSelection(archive, population.size());

        //Recombine, create offspring
        std::vector<Solution> offspring = recombine(populationMating, population.size(), objectives);

        //Mutate created offspring
        mutate(offspring, 1);

        //Evaluate offspring
        Evaluate(offspring, populationPlusArchive, objectives);

        //Replace
        std::vector<Solution> nextPopulation = population;
        for (Solution& s : offspring) {
            nextPopulation.push_back(s);
        }

        std::sort(nextPopulation.begin(), nextPopulation.end(), [](const Solution& a, const Solution& b) {
            return a.getFitness() < b.getFitness();
            });

        for (int i = 0; i < population.size(); i++) {
            population[i] = nextPopulation[i];
        }

        //Save archive as lastArchive
        lastArchive = archive;

        std::cout << "Budget 1: " << gradeAmount1 << std::endl;
        std::cout << "Budget 2: " << gradeAmount2 << std::endl;
    }

    return archive;
}

void setupObjectives(int zdt,std::vector<double (*)(const std::vector<double>& values, double parameter)>& objectives){
    gradeAmount1 = 0;
    gradeAmount2 = 0;
    objectives.clear();
    switch (zdt)
    {
        case 1:{
            objectives.push_back(&zdt1f1);
            objectives.push_back(&zdt1f2);
            break;
        }
        case 2:{
            objectives.push_back(&zdt2f1);
            objectives.push_back(&zdt2f2);
            break;
        }
        case 3:{
            objectives.push_back(&zdt3f1);
            objectives.push_back(&zdt3f2);
            break;
        }
        case 4:{
            objectives.push_back(&zdt4f1);
            objectives.push_back(&zdt4f2);
            break;
        }
        case 6:{
            objectives.push_back(&zdt6f1);
            objectives.push_back(&zdt6f2);
            break;
        }
    }
}


int main() {
    int num = 5; //number of solutions
    int n = 2; //dimensions

    //initalize objectives
    std::vector<double (*)(const std::vector<double>& values, double parameter)> objectives;
    setupObjectives(4,objectives);

    //generating solutions
    std::vector<Solution> population = generateRandom(num, n, objectives);

    std::vector<Solution> results = Spea2(population, objectives);

    int a = 2;


    // std::vector<Solution> kungResult = kungPareto(population,objectives);
    // std::ofstream file;
    // file.open("non-dominated.txt");
    // std::cout <<"Non-dominated IDs:" <<std::endl;
    // for(int i = 0; i < kungResult.size(); i++){
    //     std::cout << kungResult[i].id << " ";
    //     for(double value : kungResult[i].values){
    //         file << value << " ";
    //     }
    //     file << std::endl;

    // }
    // file.close();
    // std::cout << std::endl;

    // //get dominated solutions
    // std::vector<Solution> dominated = getDominated(population,kungResult);
    // file.open("dominated.txt");
    // std::cout <<"Dominated IDs:" <<std::endl;
    // for(Solution s : dominated){
    //     std::cout << s.id << " ";
    //     for(double value : s.values){
    //         file << value << " ";
    //     }
    //     file << std::endl;
    // }
    // file.close();

    std::ofstream outFile("results.txt");
    if (outFile.is_open()) {
        for (const auto& s : results) {

            outFile << s.objectiveScores[0] << " " << s.objectiveScores[1] << "\n";
        }
        outFile.close();
    }

    return 0;
}