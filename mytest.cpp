// CMSC 341 - Spring 2022 - Project 4
#include "dnadb.h"
#include <random>
#include <vector>
enum RANDOM {UNIFORMINT, UNIFORMREAL, NORMAL};
class Random {
public:
    Random(int min, int max, RANDOM type=UNIFORMINT, int mean=50, int stdev=20) : m_min(min), m_max(max), m_type(type)
    {
        if (type == NORMAL){
            //the case of NORMAL to generate integer numbers with normal distribution
            m_generator = std::mt19937(m_device());
            //the data set will have the mean of 50 (default) and standard deviation of 20 (default)
            //the mean and standard deviation can change by passing new values to constructor 
            m_normdist = std::normal_distribution<>(mean,stdev);
        }
        else if (type == UNIFORMINT) {
            //the case of UNIFORMINT to generate integer numbers
            // Using a fixed seed value generates always the same sequence
            // of pseudorandom numbers, e.g. reproducing scientific experiments
            // here it helps us with testing since the same sequence repeats
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_unidist = std::uniform_int_distribution<>(min,max);
        }
        else{ //the case of UNIFORMREAL to generate real numbers
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_uniReal = std::uniform_real_distribution<double>((double)min,(double)max);
        }
    }
    void setSeed(int seedNum){
        // we have set a default value for seed in constructor
        // we can change the seed by calling this function after constructor call
        // this gives us more randomness
        m_generator = std::mt19937(seedNum);
    }

    int getRandNum(){
        // this function returns integer numbers
        // the object must have been initialized to generate integers
        int result = 0;
        if(m_type == NORMAL){
            //returns a random number in a set with normal distribution
            //we limit random numbers by the min and max values
            result = m_min - 1;
            while(result < m_min || result > m_max)
                result = m_normdist(m_generator);
        }
        else if (m_type == UNIFORMINT){
            //this will generate a random number between min and max values
            result = m_unidist(m_generator);
        }
        return result;
    }

    double getRealRandNum(){
        // this function returns real numbers
        // the object must have been initialized to generate real numbers
        double result = m_uniReal(m_generator);
        // a trick to return numbers only with two deciaml points
        // for example if result is 15.0378, function returns 15.03
        // to round up we can use ceil function instead of floor
        result = std::floor(result*100.0)/100.0;
        return result;
    }
    
    private:
    int m_min;
    int m_max;
    RANDOM m_type;
    std::random_device m_device;
    std::mt19937 m_generator;
    std::normal_distribution<> m_normdist;//normal distribution
    std::uniform_int_distribution<> m_unidist;//integer uniform distribution
    std::uniform_real_distribution<double> m_uniReal;//real uniform distribution

};
class Tester{
public:

    bool normalInsert();
    bool normalCollision();
    bool invalidInsert();
    bool duplicateInsert();
    bool normalRemove();
    bool removeCollision();
    bool errorRemove();
    bool removeOverwrite();
    bool findNormal();
    bool findCollision();
    bool findError();
    bool testRehashInsert();
    bool rehashInsert();
    bool rehashRemove();
};

unsigned int hashCode(const string str);
unsigned int staticHashFunc(const string str);
string sequencer(int size, int seedNum);

bool Tester::normalInsert() {
    vector<DNA> dataList;
    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);
    bool result = true;
    for (int i=0;i<40;i++){
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        bool insertability = true;// generating random data
        for(const DNA& dna : dataList) {
            if (hashCode(dataObj.getSequence()) == hashCode(dna.getSequence())) insertability = false;
        }

        if (insertability) {
            dataList.push_back(dataObj);
            dnadb.insert(dataObj);
        }
    }
    for (auto & it : dataList){
        result = result && (it == dnadb.getDNA(it.getSequence(), it.getLocId()));
    }
    if(result && dataList.size() == dnadb.m_currentSize) return true;
    else return false;
}

bool Tester::normalCollision() {
    vector<DNA> dataList;
    DnaDb dnadb(MINPRIME, staticHashFunc);
    bool result = true;
    for (int i=0;i<50;i++){
        DNA dataObj = DNA(sequencer(5, i), MINLOCID + i);
        dataList.push_back(dataObj);
        dnadb.insert(dataObj);
    }
    for (auto & it : dataList){
        result = result && (it == dnadb.getDNA(it.getSequence(), it.getLocId()));
    }
    return result;
}

bool Tester::invalidInsert() {
    DnaDb dnaDb(MINPRIME, staticHashFunc);
    try {
        dnaDb.insert(DNA());
    } catch (exception& exception) {
        return false;
    }
    if (dnaDb.m_currentSize != 0) return false;
    return true;
}

bool Tester::duplicateInsert() {
    DnaDb dnadb(MINPRIME, staticHashFunc);
    dnadb.insert(DNA("ACGT", MINLOCID + 1));
    dnadb.insert(DNA("ACGT", MINLOCID + 1));
    return dnadb.getDNA("ACGT", MINLOCID + 1).getLocId() == (MINLOCID + 1) && dnadb.m_currentSize == 1;
}

bool Tester::normalRemove() {
    DnaDb dnadb(MINPRIME, hashCode);
    dnadb.insert(DNA("ATCAG", MINLOCID + 2));
    dnadb.insert(DNA("ATCAG", MINLOCID + 1));
    dnadb.remove(DNA("ATCAG", MINLOCID + 1));
    if(dnadb.getDNA("ATCAG", MINLOCID + 1) == EMPTY){
        return true;
    }
    return false;
}

bool Tester::removeCollision() {
    DnaDb dnadb(MINPRIME, staticHashFunc);
    dnadb.insert(DNA("ACGT", MINLOCID + 1));
    dnadb.insert(DNA("ACGT", MINLOCID + 2));
    dnadb.remove(DNA("ACGT", MINLOCID + 1));
    return dnadb.m_currentTable[0] == DELETED && dnadb.m_currentTable[1].getLocId() 
    == MINLOCID + 2 && (dnadb.m_currentSize - dnadb.m_currNumDeleted) == 1;
}

bool Tester::errorRemove() {
    DnaDb dnadb(MINPRIME, hashCode);
    dnadb.insert(DNA("ATCAG", MINLOCID + 2));
    dnadb.insert(DNA("ATCAG", MINLOCID + 1));
    dnadb.remove(DNA("ATCAG", MINLOCID + 1));
    if(dnadb.getDNA("ATCAG", MINLOCID + 1) == EMPTY){
        return true;
    }
    return false;
}

bool Tester::removeOverwrite() {
    DnaDb dnadb(MINPRIME, staticHashFunc);
    dnadb.insert(DNA("ACGT", MINLOCID + 1));
    dnadb.remove(DNA("ACGT", MINLOCID + 1));
    dnadb.insert(DNA("ACTT", MINLOCID + 2));
    return dnadb.getDNA("ACTT", MINLOCID + 2).getLocId() != 0
    && dnadb.getDNA("ACGT", MINLOCID + 1).getLocId() == 0
    && dnadb.m_currNumDeleted == 0;
}

bool Tester::findNormal() {
    DnaDb dnadb(MINPRIME, hashCode);
    dnadb.insert(DNA("ACGT", MINLOCID + 1));
    return dnadb.getDNA("ACGT", MINLOCID + 1).getLocId() == (MINLOCID + 1);
}

bool Tester::findCollision() {
    DnaDb dnadb(MINPRIME, hashCode);
    dnadb.insert(DNA("ACGT", MINLOCID + 2));
    dnadb.insert(DNA("ACGT", MINLOCID + 1));
    return dnadb.getDNA("ACGT", MINLOCID + 1).getLocId() == (MINLOCID + 1);
}

bool Tester::findError() {
    DnaDb dnadb(MINPRIME, hashCode);
    dnadb.insert(DNA("ACGT", MINLOCID + 2));
    dnadb.insert(DNA("ACGT", MINLOCID + 1));
    return dnadb.getDNA("ACTT", MINLOCID + 3) == EMPTY;
}

bool Tester::testRehashInsert() {
    DnaDb dnadb(MINPRIME, staticHashFunc);
    for (int j=0;j<60;j++){
        DNA dataObj = DNA(sequencer(5, j), MINLOCID + j);
        dnadb.insert(dataObj);
    }
    return dnadb.currState != MAXPRIME;
}

bool Tester::rehashInsert() {
    DnaDb dnadb(MINPRIME, staticHashFunc);
    for (int j=0;j<50;j++){
        DNA dataObj = DNA(sequencer(5, j), MINLOCID + j);
        dnadb.insert(dataObj);
    }
    dnadb.currState == MAXPRIME && ((float)dnadb.m_currentSize)/((float)dnadb.m_currentCap) < 0.5 && dnadb.m_currentSize == 55
        && dnadb.m_currentCap > dnadb.nextPrime(50);
    if (dnadb.currState) return true;
    else return false;//Test if rehash is successful
}
bool Tester::rehashRemove() {
    DnaDb dnadb(MINPRIME, staticHashFunc);

    for (int i=0;i<50;i++){
        DNA dataObj = DNA(sequencer(5, i), MINLOCID + i);
        dnadb.insert(dataObj);
    }
    DNA* oldTable = dnadb.m_currentTable;
    for (int j=0;j<40;j++){
        DNA dataObj = DNA(sequencer(5, j), MINLOCID + j);
        dnadb.remove(dataObj);
    }   
    if(oldTable != nullptr){
        return true;
    }
    return false;
}

unsigned int hashCode(const string str) {
   unsigned int val = 0 ;
   const unsigned int thirtyThree = 33 ;  // magic number from textbook
   for ( int i = 0 ; i < str.length(); i++)
      val = val * thirtyThree + str[i] ;
   return val ;
}

unsigned int staticHashFunc(const string str) {
    return 0;
}

string sequencer(int size, int seedNum){
    //this function returns a random DNA sequence
    string sequence = "";
    Random rndObject(0,3);
    rndObject.setSeed(seedNum);
    for (int i=0;i<size;i++){
        sequence = sequence + ALPHA[rndObject.getRandNum()];
    }
    return sequence;
}
int main(){
    Tester testing;

    if (testing.normalInsert()){
        cout << "+ normal case for insert works" << endl;
    }else cout << "- normal case for insert does not work" << endl;
    if (testing.normalCollision()){
        cout << "+ normal case for collision works" << endl;
    }else cout << "- normal case for collision does not work" << endl;
    if (testing.invalidInsert()){
        cout << "+ inserting empty DNA does not break anything" << endl;
    }else cout << "- inserting empty DNA BREAKS EVERYTHING" << endl;
    if (testing.duplicateInsert()){
        cout << "+ duplicate insert doesnt break anything" << endl;
    }else cout << "- duplicate insert BREAKS EVERYTHING" << endl;
    if (testing.normalRemove()){
        cout << "+ normal case for remove works" << endl;
    }else cout << "- normal case for remove does not work" << endl;
    if (testing.removeCollision()){
        cout << "+ collision remove works" << endl;
    }else cout << "- collision remove does not work" << endl;
    if (testing.errorRemove()){
        cout << "+ error case for remove works" << endl;
    }else cout << "- error case for remove does not work" << endl;
    if (testing.removeOverwrite()){
        cout << "+ remove case for overwrite works" << endl;
    }else cout << "- remove case for overwrite does not work" << endl;
    if (testing.findNormal()){
        cout << "+ DNA found" << endl;
    }else cout << "- DNA not found" << endl;
    if (testing.findCollision()){
        cout << "+ DNA found with collision" << endl;
    }else cout << "- DNA found with" << endl;
    if (testing.findError()){
        cout << "+ look for nonexistent DNA" << endl;
    }else cout << "- broke looking for nonexistend DNA " << endl;
    if (testing.testRehashInsert()){
        cout << "+ rehash after insert" << endl;
    }else cout << "- rehash did not trigger after insert" << endl;
    if (testing.rehashInsert()){
        cout << "+ successful rehash after insert" << endl;
    }else cout << "- unsuccessful rehash after insert " << endl;
    if (testing.rehashRemove()){
        cout << "+ rehash works after each removal" << endl;
    }else cout << "- rehash doesn't work after each removal" << endl;

    return 0;
}