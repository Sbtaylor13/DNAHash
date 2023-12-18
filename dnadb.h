// CMSC 341 - Spring 2022 - Project 4
#ifndef DNADB_H
#define DNADB_H
#include <iostream>
#include <string>
#include "math.h"
using namespace std;
class Grader;   // forward declaration, will be used for grdaing
class Tester;   // forward declaration, will be used for testing
class DNA;      // forward declaration
class DnaDb;    // forward declaration
const int MINLOCID = 1000;
const int MAXLOCID = 9999;
const int MINPRIME = 101;   // Min size for hash table
const int MAXPRIME = 99991; // Max size for hash table
#define EMPTY DNA("")
#define DELETED DNA("DELETED")
#define DELETEDKEY "DELETED"
typedef unsigned int (*hash_fn)(string); // declaration of hash function
const int MAX = 4;
const char ALPHA[MAX] = {'A', 'C', 'G', 'T'};

class DNA{
    public:
    friend class Grader;
    friend class Tester;
    friend class DnaDb;
    DNA(string sequence="", int location=0); 
    string getSequence() const;             
    int getLocId() const;
    const DNA& operator=(const DNA& rhs);
    friend ostream& operator<<(ostream& sout, const DNA &dna );
    friend bool operator==(const DNA& lhs, const DNA& rhs);
    private:
    string m_sequence; 
    int m_location; 
};

class DnaDb{
    public:
    friend class Grader;
    friend class Tester;
    DnaDb(int size, hash_fn hash);
    ~DnaDb();
    float deletedRatio() const;
    bool insert(DNA dna);
    bool remove(DNA dna);
    DNA getDNA(string sequence, int location);
    void dump() const;

    private:
    hash_fn         m_hash;         // hash function

    DNA*            m_currentTable; // hash table
    unsigned int    m_currentCap;   // hash table size
    unsigned int    m_currentSize;  // current number of entries
                                    // m_currentSize includes deleted entries 
    unsigned int    m_currNumDeleted;// number of deleted entries

    DNA*            m_oldTable;     // hash table
    unsigned int    m_oldCap;       // hash table size
    unsigned int    m_oldSize;      // current number of entries
                                    // m_oldSize includes deleted entries
    unsigned int    m_oldNumDeleted;// number of deleted entries

    //private helper functions
    bool isPrime(int number);
    int nextPrime(int current);

    /******************************************
    * Private function declarations go here! *
    ******************************************/

    unsigned int findSlot(unsigned proposed);
    unsigned int search(const string &seq, const int &location, unsigned int current, unsigned run = 1);
    void rehashHelper();
    void moveSome();
    bool lookHelper(DNA dna, bool checkOld);
    unsigned currState = 0;
    DNA* theTarget;
    unsigned int theCap; 

};
#endif