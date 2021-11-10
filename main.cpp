// Created by Eyu on 11/7/2021
#include <tuple>
#include <vector>
#include <iostream>
#include "Course.h"
#include <fstream>
#include <ctime>
#include <algorithm>
#include <cstdlib>
#include <chrono>

using namespace std;

/**
* Bubble Sort algorithm
 * Second improvement. Reduce number of comparisons in inner loop
 * as after every iteration, larger number will be sorted correctly
 * only the rest should be compared
*/
template<typename Comparable>
tuple<unsigned long int, unsigned long int> bubbleSort(vector<Comparable>& v) {
    unsigned long int numReads = 0;
    unsigned long int numWrites = 0;
    bool haveSwapped = true;
    int numPasses = 0;
    while (haveSwapped) {
        haveSwapped = false;
        for (int i = 0; i < (v.size() - numPasses - 1); ++i) {
            if (v[i] > v[i + 1]) {
                // The two elements are out of order. Swap them.
                Comparable temp = v[i];
                v[i] = v[i + 1];
                v[i + 1] = temp;
                haveSwapped = true;
                numWrites += 3;
                numReads += 3;
            }
            numReads += 2; //reads for comparison
        }
        ++numPasses;
        //printVector(v);
    }
    return make_tuple(numReads, numWrites);
}

/**
 * Selection Sort algorithm
 * after every iteration, the smaller value will be placed to the left
 **/
template<typename Comparable>
tuple<unsigned long int, unsigned long int> selectionSort(vector<Comparable>& v) {
    unsigned long int numReads = 0;
    unsigned long int numWrites = 0;
    for (int swapIndex = 0; swapIndex < v.size(); ++swapIndex) {
        int minIndex = swapIndex;
        for (int i = swapIndex + 1; i < v.size(); ++i) {
            if (v[i] < v[minIndex]) {
                minIndex = i;
            }
            numReads += 2;
        }
        Comparable temp = v[minIndex];
        v[minIndex] = v[swapIndex];
        v[swapIndex] = temp;
        numWrites += 3;
        numReads += 3;
    }
    return make_tuple(numReads, numWrites);
}

/**
 * Insertion sort
 * return the number of reads and writes*/
template<typename Comparable>
tuple<unsigned long int, unsigned long int> insertionSort(vector<Comparable> vec) {
    int unsortedStartIndex, insertIndex;
    unsigned long int numReads = 0, numWrites = 0;
    Comparable toBeInserted;
    for (unsortedStartIndex = 1; unsortedStartIndex < vec.size(); ++unsortedStartIndex) {
        toBeInserted = vec[unsortedStartIndex];
        ++numWrites;
        ++numReads;
        // Loop to shift over the larger elements
        insertIndex = unsortedStartIndex - 1;
        while (insertIndex >= 0 && vec[insertIndex] > toBeInserted) {
            vec[insertIndex + 1] = vec[insertIndex];
            numReads += 3;
            ++numWrites;
            --insertIndex;
        }
        // Put toBeInserted back into vec
        vec[insertIndex + 1] = toBeInserted;
        ++numWrites;
        ++numReads;
    }
    return make_tuple(numReads, numWrites);
}

/**
 * Insertion sort by another field
 * return the number of reads and writes*/
template<typename Comparable>
tuple<unsigned long int, unsigned long int> insertionSort_CRN(vector<Comparable> vec) {
    int unsortedStartIndex, insertIndex;
    unsigned long int numReads = 0, numWrites = 0;
    Comparable toBeInserted;
    for (unsortedStartIndex = 1; unsortedStartIndex < vec.size(); ++unsortedStartIndex) {
        toBeInserted = vec[unsortedStartIndex];
        ++numWrites;
        ++numReads;
        // Loop to shift over the larger elements
        insertIndex = unsortedStartIndex - 1;
        while (insertIndex >= 0 && vec[insertIndex].getCRN() > toBeInserted.getCRN()) {
            vec[insertIndex + 1] = vec[insertIndex];
            numReads += 3;
            ++numWrites;
            --insertIndex;
        }
        // Put toBeInserted back into vec
        vec[insertIndex + 1] = toBeInserted;
        ++numWrites;
        ++numReads;
    }
    return make_tuple(numReads, numWrites);
}

/**
 * Quicksort stable, recursive algorithm
 * return the number of reads and writes
 * note 3 vector smaller, equal and larger to place the objects and push back to the
 * original one.
 **/
template<typename Comparable>
tuple<unsigned long int, unsigned long int> quicksort(vector<Comparable>& v) {
    unsigned long int numWrites = 0;
    unsigned long int numReads = 0;
    // Recursive base case
    if (v.size() < 2) {
        return make_tuple(numReads, numWrites);
    }
    Comparable pivot = v[0];
    vector<Comparable> smaller, equal, larger;
    ++numReads;
    numWrites += 4;
    int i;  // OK to declare here. See where we copy values back to orig. vector
    for (i = 0; i < v.size(); ++i) {
        if (v[i] < pivot) {
            smaller.push_back(v[i]);
        } else if (v[i] > pivot) {
            larger.push_back(v[i]);
        } else {
            equal.push_back(v[i]);
        }
        numReads += 2; // 2 reads for comparison and push_back
        ++numWrites; // one write for push_back
    }

    auto quickSmall = quicksort(smaller);
    auto quickLarge = quicksort(larger);
    numReads += get<0>(quickSmall) + get<0>(quickLarge);
    numWrites += get<1>(quickSmall) + get<1>(quickLarge);

    // Copy everything back into vector
    for (i = 0; i < smaller.size(); ++i) {
        v[i] = smaller[i];
        ++numReads;
        ++numWrites;
    }
    // Omit initialization of i and continue...
    for (; i < smaller.size() + equal.size(); ++i) {
        v[i] = equal[i - smaller.size()];
        ++numReads;
        ++numWrites;
    }
    // Omit initialization of i and continue...
    for (; i < v.size(); ++i) {
        v[i] = larger[i - smaller.size() - equal.size()];
        ++numReads;
        ++numWrites;
    }
    return make_tuple(numReads, numWrites);
}

/**
 * Merge sort*/
template<typename Comparable>
tuple<unsigned long int, unsigned long int> mergeSortRec(vector<Comparable> &vec, int startIndex, int endIndex) {
    unsigned long int numReads = 0, numWrites = 0;
    // Recursive base case
    if (startIndex == endIndex) {
        // We have one item. There is nothing to split or sort.
        return make_tuple(numReads, numWrites);
    }

    // Recursive calls
    int centerIndex = (startIndex + endIndex) / 2;
    mergeSortRec(vec, startIndex, centerIndex);
    mergeSortRec(vec, centerIndex + 1, endIndex);

    // Merge
    vector<Comparable> temp;
    int leftIndex = startIndex;
    int rightIndex = centerIndex + 1;
     //write for temp
    while (leftIndex <= centerIndex && rightIndex <= endIndex) {
        numReads += 3; //read when comparing and push back
        ++numWrites; // write when push back
        if (vec[leftIndex] <= vec[rightIndex]) {
            temp.push_back(vec[leftIndex]);
            ++leftIndex;
        } else {
            temp.push_back(vec[rightIndex]);
            ++rightIndex;
        }
    }
    // At this point, one of the halves has been completely copied into temp but the other hasn't
    // We need to finish copying everything into temp, so we make loops for each half
    while (leftIndex <= centerIndex) {
        temp.push_back(vec[leftIndex]);
        ++numReads; //read when push back
        ++numWrites; // write when push back
        ++leftIndex;
    }
    while (rightIndex <= endIndex) {
        temp.push_back(vec[rightIndex]);
        ++numReads; //read when push back
        ++numWrites; // write when push back
        ++rightIndex;
    }
    // At this point, all of the items from startIndex to endIndex have been copied into temp
    // Copy everything from temp back into vec
    for (int i = 0; i < temp.size(); ++i) {
        vec[i + startIndex] = temp[i];
        ++numReads;
        ++numWrites;
    }

    return make_tuple(numReads, numWrites);
}

template<typename Comparable>
tuple<unsigned long int, unsigned long int> mergeSort(vector<Comparable> vec) {
    return mergeSortRec(vec, 0, vec.size() - 1);
}


/**
 * Heap sort algorithm
* Helper function for percolateDown
* Given some index i, the index of the left child
* of i is given by 2i + 1
*/
inline int leftChild(int i) {
    return 2 * i + 1;
}

/**
 * Percolate down; used by heapSort() and heapify()
 * We start at the root, and keep swapping until we find the
 * correct place for the value
 */
template<typename Comparable>
tuple<unsigned long long int, unsigned long long int> percolateDown(vector<Comparable> &vec, int i, int n, int child, Comparable temp) {
    // we start with root = start, but as we progress,
    // "root" will actually be the root of a subtree
    unsigned long long int numReads = 0, numWrites = 0;
    for (temp = vec[i]; leftChild(i) < n; i = child) {
        numReads += 1;
        child = leftChild(i);
        // choose the child with the larger value
        if (child != n - 1 && vec[child] < vec[child + 1]) {
            numReads += 2;
            ++child;
        }
        // if the parent is less than the child, swap them
        if (temp < vec[child]) {
            vec[i] = vec[child];
            numReads += 3;
            numWrites += 1;
        } else {
            // parent is >= both children. nothing more to do.
            break;
        }
    }
    vec[i] = temp;
    numReads += 1;
    numWrites += 1;
    return make_tuple(numReads, numWrites);
}


/**
 * Heap sort
 */
template<typename Comparable>
tuple<unsigned long long int, unsigned long long int> heapSort(vector<Comparable>& vec) {
    unsigned long long int numReads = 0, numWrites = 0;
    int i, j, child;
    Comparable temp1, temp2;
    // build the heap (with max value at root)
    for (i = vec.size() / 2 - 1; i >= 0; --i) {
        auto heapRecord_0=percolateDown(vec, i, vec.size(), child, temp2);
        numReads += get<0>(heapRecord_0);
        numWrites += get<1>(heapRecord_0);
    }
    // keep deleting the max
    for (j = vec.size() - 1; j > 0; --j) {
        // swap the maximum out
        temp1 = vec[0];
        vec[0] = vec[j];
        vec[j] = temp1;
        // Swap need 3 RW
        numReads += 3;
        numWrites += 3;

        // make it into a heap again
        auto heapRecord_1 = percolateDown(vec, 0, j, child, temp2);
        numReads += get<0>(heapRecord_1);
        numWrites += get<1>(heapRecord_1);
    }
    return make_tuple(numReads, numWrites);
}

/**2-sort
 * we use selection sort first and then insertion sort
 * return the number of reads and writes*/
template<typename Comparable>
tuple<unsigned long int, unsigned long int> twoSort(vector<Comparable>& v) {
    auto records_sel = selectionSort(v);
    auto records_ins = insertionSort_CRN(v);
    return make_tuple(get<0>(records_sel)+get<0>(records_ins), get<1>(records_sel)+get<1>(records_ins));
}

/**
 * function to write the data into file
 **/
void writeToFile(unsigned long int size, unsigned long int rBubb,unsigned long int wBubb,
                 unsigned long int rSele, unsigned long int wSele,
                 unsigned long int rQui, unsigned long int wQui,
                 unsigned long long int rHeap, unsigned long long int wHeap,
                 unsigned long int r2sort, unsigned long int w2sort,
                 unsigned long int rInser, unsigned long int wInser,
                 unsigned long int rMerge, unsigned long int wMerge){
    ofstream outFile;
    outFile.open("sorting_algorithm.csv", ios::app);
    outFile << size << ", "
            << rBubb << ", " << wBubb << ", "
            << rSele << ", " << wSele << ", "
            << rQui << ", " << wQui << ", "
            << rHeap << ", " << wHeap << ", "
            << r2sort << ", " << w2sort << ", "
            << rInser << ", " << wInser << ", "
            << rMerge << ", " << wMerge << ", "
            << endl;
    outFile.close();
}

/**
 * function to achieve vectors in different size
 * @param size: the size of the vector to sort
 * @param mainVector: the main vector that stores all the data
 **/
template<typename Comparable>
vector<Comparable> getUnsortedVector(int size, vector<Comparable> mainVector) {
    vector<Comparable> vectorUnsorted;
    for(int i=0; i < size; ++i) {
        vectorUnsorted.push_back(mainVector.at(i));
    }
    return vectorUnsorted;
}

/**
 * function to get the number of reads and writes from 5 algorithms with different size
 * @param vectorUnsorted: the unsorted vector ps.the vector would be sorted once algorithm finished
 **/
template<typename Comparable>
void sortedBy(int size, vector<Comparable> vectorUnsorted) {
    // define the 5 vectors for each algorithm
     vector<Comparable> vectorUnsortedCopy1 = vectorUnsorted;
     vector<Comparable> vectorUnsortedCopy2 = vectorUnsorted;
     vector<Comparable> vectorUnsortedCopy3 = vectorUnsorted;
     vector<Comparable> vectorUnsortedCopy4 = vectorUnsorted;
     vector<Comparable> vectorUnsortedCopy5 = vectorUnsorted;
     vector<Comparable> vectorUnsortedCopy6 = vectorUnsorted;
     vector<Comparable> vectorUnsortedCopy7 = vectorUnsorted;
    // get the num of reads and writes
    // timer start
     auto start1 = chrono::high_resolution_clock::now();
     auto records1 = bubbleSort(vectorUnsortedCopy1);
     auto stop1 = chrono::high_resolution_clock::now();
     auto duration1 = chrono::duration_cast<chrono::microseconds>(stop1 - start1);
     cout << "Bubble Sort of " << size << "records: " << duration1.count() << "microseconds" << endl;

     auto start2 = chrono::high_resolution_clock::now();
     auto records2 = selectionSort(vectorUnsortedCopy2);
     auto stop2 = chrono::high_resolution_clock::now();
     auto duration2 = chrono::duration_cast<chrono::microseconds>(stop2 - start2);
     cout << "Selection Sort of " << size << "records: " << duration2.count() << "microseconds" << endl;

     auto start3 = chrono::high_resolution_clock::now();
     auto records3 = quicksort(vectorUnsortedCopy3);
     auto stop3 = chrono::high_resolution_clock::now();
     auto duration3 = chrono::duration_cast<chrono::microseconds>(stop3 - start3);
     cout << "Quick Sort of " << size << "records: " << duration3.count() << "microseconds" << endl;

     auto start4 = chrono::high_resolution_clock::now();
     auto records4 = heapSort(vectorUnsortedCopy4);
     auto stop4 = chrono::high_resolution_clock::now();
     auto duration4 = chrono::duration_cast<chrono::microseconds>(stop4 - start4);
     cout << "Heap Sort of " << size << "records: " << duration4.count() << "microseconds" << endl;

     auto start5 = chrono::high_resolution_clock::now();
     auto records5 = twoSort(vectorUnsortedCopy5);
     auto stop5 = chrono::high_resolution_clock::now();
     auto duration5 = chrono::duration_cast<chrono::microseconds>(stop5 - start5);
     cout << "Two Sort of " << size << "records: " << duration5.count() << "microseconds" << endl;

     auto start6 = chrono::high_resolution_clock::now();
     auto records6 = insertionSort(vectorUnsortedCopy6);
     auto stop6 = chrono::high_resolution_clock::now();
     auto duration6 = chrono::duration_cast<chrono::microseconds>(stop6 - start6);
     cout << "Insertion Sort of " << size << "records: " << duration6.count() << "microseconds" << endl;

     auto start7 = chrono::high_resolution_clock::now();
     auto records7 = mergeSort(vectorUnsortedCopy7);
     auto stop7 = chrono::high_resolution_clock::now();
     auto duration7 = chrono::duration_cast<chrono::microseconds>(stop7 - start7);
     cout << "Merge Sort of " << size << "records: " << duration7.count() << "microseconds" << endl;

     // auto records5 = for 2-sort
     // writing into the file
    writeToFile(size, get<0>(records1), get<1>(records1),
                get<0>(records2), get<1>(records2),
                get<0>(records3), get<1>(records3),
                get<0>(records4), get<1>(records4),
                get<0>(records5), get<1>(records5),
                get<0>(records6), get<1>(records6),
                get<0>(records7), get<1>(records7));

 }

int main() {
    // Reda the input file and create vector of objects
    vector<Course> course;
    string fn = "../uvm_fall2021_courses.csv";

    if (loadFromFile(fn, course)) {
        cout << course.size()
             << " records read from file" << endl;
    } else {
        cout << "Something went wrong." << endl;
    }

    //open a file
    ofstream sortRecord;
    sortRecord.open("sorting_algorithm.csv", ios::app);
    sortRecord << "Size" << ", "
               << "Reads_Bubble" << ", " << "Writes_Bubble" << ", "
               << "Reads_Selec" << ", " << "Writes_Selec" << ", "
               << "Reads_Quick" << ", " << "Writes_Quick" << ", "
               << "Reads_Heap"  << ", " << "Writes_Heap" << ", "
               << "Reads_2" << ", " << "Writes_2" << ", "
               << "Reads_Insert" << ", " << "Writes_Insert" << ", "
               << "Reads_Merge" << ", " << "Writes_Merge" 
               << endl;
    sortRecord.close();

    // sorting a vector of size 100
    sortedBy(100, getUnsortedVector(100, course));
    sortedBy(200, getUnsortedVector(200, course));
    sortedBy(300, getUnsortedVector(300, course));
    sortedBy(400, getUnsortedVector(400, course));
    sortedBy(500, getUnsortedVector(500, course));
    sortedBy(600, getUnsortedVector(600, course));
    sortedBy(700, getUnsortedVector(700, course));
    sortedBy(800, getUnsortedVector(800, course));
    sortedBy(900, getUnsortedVector(900, course));
    sortedBy(1000, getUnsortedVector(1000, course));

    return 0;


}