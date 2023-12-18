// CMSC 341 - Spring 2022 - Project 4
#include "dnadb.h"

DnaDb::DnaDb(int size, hash_fn hash){
    m_hash = hash;
    m_currentCap = size;
    
    if (!isPrime((int)m_currentCap)){
        m_currentCap = nextPrime((int)m_currentCap);// if not prime go to next prime
    } 
    
    m_currentTable = new DNA[size];
    for (int i = 0; i < size; ++i) {
        m_currentTable[i] = EMPTY;
    }
    m_oldTable = nullptr;
    m_oldCap = 0;
    m_oldSize = 0;
    m_currentCap = size;
    m_currentSize = 0;
    m_currNumDeleted = 0;
    m_oldNumDeleted = 0;
    currState = MAXPRIME;
}

DnaDb::~DnaDb(){//destructor
    if (m_oldTable) delete[] m_oldTable;
    if (m_currentTable) delete[] m_currentTable;
    m_currentCap = 0;
    m_currNumDeleted = 0;
    m_currentSize = 0;
    m_oldNumDeleted = 0;
    m_oldCap = 0;
    m_oldSize = 0;
    currState = MAXPRIME;
}

bool DnaDb::insert(DNA dnaPass){
    bool result = lookHelper(dnaPass, true);

    if (((float)m_currentSize)/((float)m_currentCap) > 0.5 && currState == MAXPRIME){
        rehashHelper();//DNA object goes to hash table
    } 
    if (currState != MAXPRIME){
        moveSome();
    } 
    return result;
}

bool DnaDb::lookHelper(DNA dna, bool oldChecker){
    if (dna.m_location > MAXLOCID || dna.m_location < MINLOCID){
        return false;
    }
    theTarget = m_currentTable;
    theCap = m_currentCap;
    if (m_currentCap > 0){
        int searching = search(dna.m_sequence, dna.m_location, m_hash(dna.m_sequence) % theCap);
        if (searching != MAXPRIME){
            return false;
        } 
    }
    theTarget = m_oldTable;
    theCap = m_oldCap;
    if (oldChecker && theCap > 0)
        if (search(dna.m_sequence, dna.m_location, m_hash(dna.m_sequence) % theCap) != MAXPRIME) return false;

    if (m_currentTable && m_currentCap > 0) {
        int proposedIns = findSlot(m_hash(dna.getSequence()) % m_currentCap);
        if (m_currentTable[proposedIns] == DELETED){
            m_currNumDeleted--;
        } 
        else{
            m_currentSize++;
        } 
        m_currentTable[proposedIns] = dna;
        return true;
    }
    return false;
}

bool DnaDb::remove(DNA dna){
    theTarget = m_currentTable;
    theCap = m_currentCap;
    unsigned id = search(dna.m_sequence, dna.m_location, m_hash(dna.m_sequence) % theCap);
    if (id != MAXPRIME) {
        theTarget[id] = DELETED;//it is now part of the deleted object
        m_currNumDeleted++;
        if (deletedRatio() > 0.8 && currState == MAXPRIME){
            rehashHelper();//ratio of 80%
        } 
        if (currState != MAXPRIME) moveSome();
        return true;
    }

    theTarget = m_oldTable;
    theCap = m_oldCap;
    if(theCap >= MINPRIME){
        id = search(dna.m_sequence, dna.m_location, m_hash(dna.m_sequence) % theCap);
    }
    if (id != MAXPRIME) {
        theTarget[id] = DELETED;
        m_oldNumDeleted++;
        if (deletedRatio() > 0.8 && currState == MAXPRIME) rehashHelper();
        if (currState != MAXPRIME) moveSome();
        return true;
    }

    return false;
}

float DnaDb::deletedRatio() const {
    return ((float)m_currNumDeleted)/((float)m_currentSize);
}

unsigned DnaDb::findSlot(unsigned given) {
    unsigned run = 1;
    while(!(m_currentTable[given] == EMPTY) || m_currentTable[given] == DELETED) {
        if (run == m_currentCap) return MAXPRIME;
        given = (given + (run*run)) % m_currentCap;
        run += 1;
    }
    return given;
}

unsigned DnaDb::search(const string &seq, const int &location, unsigned int current, unsigned int run) {
    if (run >= theCap){
        return MAXPRIME;
    } 
    if ( theTarget[current].m_sequence == seq && theTarget[current].m_location == location){
        return current;
    } 
    
    return search(seq, location, (current + run*run) % theCap, run + 1);
}

void DnaDb::moveSome() {
    int pos = 0;
    for (int i = 0; i <= currState && pos < m_oldCap; ++i) {
        if (m_oldTable[pos] == DELETED || m_oldTable[pos] == EMPTY){
            i--;
            pos++;
        } else {
            lookHelper(m_oldTable[pos], false);
            m_oldTable[pos] = DELETED;
            m_oldNumDeleted++;
            pos++;
        }
    }
    if (m_oldNumDeleted == m_oldSize) {
        delete[] m_oldTable;
        m_oldTable = nullptr;
        m_oldCap = 0;
        m_oldNumDeleted = 0;
        m_oldSize = 0;
        currState = MAXPRIME;
    }//resets, frees memory
}

bool operator==(const DNA& x, const DNA& y){
    return ((x.m_sequence == y.m_sequence) && (x.m_location == y.m_location));
}
bool DnaDb::isPrime(int number){
    bool result = true;
    for (int i = 2; i <= number / 2; ++i) {
        if (number % i == 0) {
            result = false;
            break;
        }
    }
    return result;
}
DNA DnaDb::getDNA(string sequence, int location){
    if (m_currentTable && m_currentCap >= MINPRIME) {
        theTarget = m_currentTable;
        theCap = m_currentCap;
        if (search(sequence, location, m_hash(sequence) % theCap) != MAXPRIME) {
            return theTarget[search(sequence, location, m_hash(sequence) % theCap)];
        }//searches and returns sequence
    }

    if (m_oldTable && m_oldCap >= MINPRIME) {
        theTarget = m_oldTable;
        theCap = m_oldCap;
        if (search(sequence, location, m_hash(sequence) % theCap) != MAXPRIME){
            return theTarget[search(sequence, location, m_hash(sequence) % theCap)];
        }
    }
    return EMPTY;
}
int DnaDb::nextPrime(int current){
    if (current < MINPRIME) current = MINPRIME-1;
    for (int i=current; i<MAXPRIME; i++) {
        for (int j=2; j*j<=i; j++) {
            if (i % j == 0)
                break;
            else if (j+1 > sqrt(i) && i != current) {
                return i;
            }
        }
    }
    return MAXPRIME;
}

DNA::DNA(string sequence, int location) {
    if ((location >= MINLOCID && location <= MAXLOCID) ||
        (location == 0 && sequence == "DELETED")){
        m_sequence = sequence;
        m_location = location;
    }
    else{
        m_sequence = "";
        m_location = 0;
    }
}
void DnaDb::rehashHelper() {
    m_oldTable = m_currentTable;
    unsigned sizeCap = nextPrime((m_currentSize - m_currNumDeleted) * 4);
    m_currentTable = new DNA[sizeCap];
    for (int i = 0; i < sizeCap; ++i) {
        m_currentTable[i] = EMPTY;
    }
    m_currentCap = sizeCap;
    m_currentSize = 0;
    m_oldCap = m_currentCap;
    m_oldNumDeleted = m_currNumDeleted;
    m_oldSize = m_currentSize;
    currState = (unsigned)max(1.0f,floor(((float)m_oldSize)/float(4.0)));
    m_currNumDeleted = 0;
}
string DNA::getSequence() const {
    return m_sequence;
}

const DNA& DNA::operator=(const DNA& rhs){//assignment operator
    if (this != &rhs){
        m_sequence = rhs.m_sequence;
        m_location = rhs.m_location;
    }
    return *this;
}

int DNA::getLocId() const {
    return m_location;
}
ostream& operator<<(ostream& sout, const DNA &dna ) {//overloaded insertion, prints key and position
    if (!dna.m_sequence.empty())
        cout << dna.m_sequence << " (Location ID " << dna.m_location << ")";
    else
        cout << "";
  return sout;
}

void DnaDb::dump() const {
    cout << "Dump for current table: " << endl;
    if (m_currentTable != nullptr)
        for (int i = 0; i < m_currentCap; i++) {
            cout << "[" << i << "] : " << m_currentTable[i] << endl;
        }
    cout << "Dump for old table: " << endl;
    if (m_oldTable != nullptr)
        for (int i = 0; i < m_oldCap; i++) {
            cout << "[" << i << "] : " << m_oldTable[i] << endl;
        }
}