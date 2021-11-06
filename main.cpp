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
tuple<int, int> bubbleSort(vector<Comparable>& v) {
    int numReads = 0;
    bool haveSwapped = true;
    int numPasses = 0;
    int numWrites = 0;
    while (haveSwapped) {
        haveSwapped = false;
        for (int i = 0; i < (v.size() - numPasses - 1); ++i) {
            numReads += 2;
            if (v[i] > v[i + 1]) {
                // The two elements are out of order. Swap them.
                Comparable temp = v[i];
                v[i] = v[i + 1];
                v[i + 1] = temp;
                haveSwapped = true;
                numWrites += 3;
                numReads += 3;
            }
        }
        ++numPasses;
        //printVector(v);
    }
    return make_tuple(numReads, numWrites);
}

/**
 * Selection Sort algorithm
 * after every iteration, the smaller value will be placed to the left*/
template<typename Comparable>
tuple<int, int> selectionSort(vector<Comparable>& v) {
    int numReads = 0;
    int numWrites = 0;
    for (int swapIndex = 0; swapIndex < v.size(); ++swapIndex) {
        int minIndex = swapIndex;
        for (int i = swapIndex + 1; i < v.size(); ++i) {
            numReads += 2;
            if (v[i] < v[minIndex]) {
                minIndex = i;
            }
        }
        Comparable temp = v[minIndex];
        v[minIndex] = v[swapIndex];
        v[swapIndex] = temp;
        numWrites += 3;
        numReads += 3;
        //printVector(v);
    }
    return make_tuple(numReads, numWrites);
}

/**
 * Quicksort stable, recursive algorithm
 * return the number of reads and writes
 * note 3 vector smaller, equal and larger to place the objects and push back to the
 * original one.
 */
template<typename Comparable>
tuple<int, int> quicksort(vector<Comparable>& v) {
    int numWrites = 0;
    int numReads = 0;
    // Recursive base case
    if (v.size() < 2) {
        return make_tuple(numReads, numWrites);
    }
    Comparable pivot = v[0];
    ++numReads;
    ++numWrites;
    vector<Comparable> smaller, equal, larger;
    int i;  // OK to declare here. See where we copy values back to orig. vector
    for (i = 0; i < v.size(); ++i) {
        numReads += 2;
        ++numWrites;
        if (v[i] < pivot) {
            smaller.push_back(v[i]);
        } else if (v[i] > pivot) {
            larger.push_back(v[i]);
        } else {
            equal.push_back(v[i]);
        }
    }

    auto quickSmall = quicksort(smaller);
    auto quickLarge = quicksort(larger);
    numReads += get<0>(quickSmall) + get<0>(quickLarge);
    numWrites += get<1>(quickSmall) + get<1>(quickLarge);

    // Copy everything back into vector
    for (i = 0; i < smaller.size(); ++i) {
        v[i] = smaller[i];
    }
    // Omit initialization of i and continue...
    for (; i < smaller.size() + equal.size(); ++i) {
        v[i] = equal[i - smaller.size()];
    }
    // Omit initialization of i and continue...
    for (; i < v.size(); ++i) {
        v[i] = larger[i - smaller.size() - equal.size()];
    }
    numReads += 2 * v.size();
    numWrites += v.size();
    return make_tuple(numReads, numWrites);
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
tuple<int, int> percolateDown(vector<Comparable>& v, int start, int end) {
    // we start with root = start, but as we progress,
    // "root" will actually be the root of a subtree
    int numReads = 0;
    int numWrites = 0;
    int root = start;
    while (leftChild(root) <= end) {
        int child = leftChild(root);
        if (child + 1 <= end && v[child] < v[child + 1]){
            // if there's a right child, and the value of the
            // left child is less than that of the right child...
            // we increment child, so we're working with
            // the right child.
            ++child;
            numReads += 2;
        }
        if (v[root] < v[child]) {
            // swap
            numWrites += 3;
            numReads += 3;
            Comparable temp = v[root];
            v[root] = v[child];
            v[child] = temp;
            root = child;
        } else {
            int a = numReads;
            int b = numWrites;
            return make_tuple(numReads, numWrites);
        }
    }
}

/**
 * heapify; used by heapSort()
 * If k is the number of levels in our heap, we start
 * at the k - 1 level (one up from the bottom) and we
 * percolate elements down, starting from the appropriate
 * node in that level -- that is, the rightmost element
 * with children. start = (size - 2) / 2 will yield
 * the index of that element. After each percolation is
 * complete, we decrement start, moving left, and then up
 * through the heap. In the last step, start = 0 and we
 * percolate down from the root to complete the process.
 */
template<typename Comparable>
tuple<int, int> heapify(vector<Comparable>& v, int size) {
    int numReadsRecord = 0;
    int numWritesRecord = 0;
    int start = (size - 2) / 2;
    while (start >= 0) {
        auto heapRecord = percolateDown(v, start, size - 1);
        numReadsRecord = get<0>(heapRecord);
        numWritesRecord = get<1>(heapRecord);
        --start;
    }
    return make_tuple(numReadsRecord, numWritesRecord);

}

/**
 * Heap sort
 */
template<typename Comparable>
tuple<int, int> heapSort(vector<Comparable>& v) {
    int size = v.size();
    int end = size - 1;
    auto heapRecord = heapify(v, size);
    int numReads = get<0>(heapRecord);
    int numWrites = get<1>(heapRecord);
    //std::cout << "Heapify" << std::endl;
    //printVector(v);
    //std::cout << "---" << std::endl;
    while (end > 0) {
        // v[end] is sorted portion; v[0] is current max
        numReads += 3;
        numWrites += 3;
        Comparable temp = v[end];
        v[end] = v[0];
        v[0] = temp;
        --end;
        auto heapRecord2 = percolateDown(v, 0, end);
        numReads += get<0>(heapRecord2);
        numWrites += get<1>(heapRecord2);
        //printVector(v);
    }
    return make_tuple(numReads, numWrites);
}

void writeToFile(int size, int rBubb,int wBubb,
                    int rSele, int wSele, int rQui, int wQui, int rHeap, int wHeap,
                    int r2sort, int w2sort){
    ofstream outFile;
    outFile.open("sorting_algorithm.csv", ios::app);
    outFile << size << ", "
            << rBubb << ", " << wBubb << ", "
            << rSele << ", " << wSele << ", "
            << rQui << ", " << wQui << ", "
            << rHeap << ", " << wHeap << ", "
            << r2sort << ", " << w2sort 
            << endl;
    outFile.close();
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
            << "Reads_2" << ", " << "Writes_2"
            << endl;
    sortRecord.close();

    // sorting a vector of size 100
    vector<Course> v1a, v1b, v1c, v1d, v1e;
    for (int i = 0; i < 100; ++i){
        v1a.push_back(course.at(i));
        v1b.push_back(course.at(i));
        v1c.push_back(course.at(i));
        v1d.push_back(course.at(i));
        v1e.push_back(course.at(i));
    }
    auto start1 = chrono::high_resolution_clock::now();
    auto bubbleSortRecord1 = bubbleSort(v1a);
    auto stop1 = chrono::high_resolution_clock::now();
    auto duration1 = chrono::duration_cast<chrono::microseconds>(stop1 - start1);
    cout << "Bubble Sort: "
              << duration1.count() << " microseconds" << endl;
    cout << v1a[0].getCRN() << ", " << v1a[5].getCRN() << ", " << v1a[10].getCRN()
            << ", " << v1a[20].getCRN() << endl;

    int rB100 = get<0>(bubbleSortRecord1);
    int wB100 = get<1>(bubbleSortRecord1);

    auto start2 = chrono::high_resolution_clock::now();
    auto selecSortRecord1 = selectionSort(v1b);
    auto stop2 = chrono::high_resolution_clock::now();
    auto duration2 = chrono::duration_cast<chrono::microseconds>(stop2 - start2);
    cout << "Selection Sort: "
         << duration2.count() << " microseconds" << endl;
    cout << v1b[0].getCRN() << ", " << v1b[5].getCRN() << ", " << v1b[10].getCRN()
         << ", " << v1b[20].getCRN() << endl;
    int rS100 = get<0>(selecSortRecord1);
    int wS100 = get<1>(selecSortRecord1);

    auto start3 = chrono::high_resolution_clock::now();
    auto quickSortRecord1 = quicksort(v1c);
    auto stop3 = chrono::high_resolution_clock::now();
    auto duration3 = chrono::duration_cast<chrono::microseconds>(stop3 - start3);
    cout << "Quick Sort: "
         << duration3.count() << " microseconds" << endl;
    cout << v1c[0].getCRN() << ", " << v1c[5].getCRN() << ", " << v1c[10].getCRN()
         << ", " << v1c[20].getCRN() << endl;
    int rQ100 = get<0>(quickSortRecord1);
    int wQ100 = get<1>(quickSortRecord1);

    auto start4 = chrono::high_resolution_clock::now();
    auto heapSortRecord1 = heapSort(v1d);
    auto stop4 = chrono::high_resolution_clock::now();
    auto duration4 = chrono::duration_cast<chrono::microseconds>(stop4 - start4);
    cout << "Heap Sort: "
         << duration4.count() << " microseconds" << endl;
    cout << v1d[0].getCRN() << ", " << v1d[5].getCRN() << ", " << v1d[10].getCRN()
         << ", " << v1d[20].getCRN() << endl;
    int rH100 = get<0>(heapSortRecord1);
    int wH100 = get<1>(heapSortRecord1);

    // 2 sort -- selection sort and then bubble sort on RowId
    auto start5 = chrono::high_resolution_clock::now();
    auto TwoSortRecord1 = selectionSort(v1e);
    int r2100 = get<0>(TwoSortRecord1);
    int w2100 = get<1>(TwoSortRecord1);
    int num2Reads1 = 0;
    int num2Writes1 = 0;
    for (int j = 0; j < v1e.size(); ++j) {
        for (int i = 0; i < v1e.size() - 1; ++i) {
            num2Reads1 += 2;
            if (v1e[i].getRowId() > v1e[i + 1].getRowId()) {
                // The two elements are out of order. Swap them.
                num2Writes1 += 3;
                num2Reads1 += 3;
                Course temp = v1e[i];
                v1e[i] = v1e[i + 1];
                v1e[i + 1] = temp;
            }
        }
    }
    auto stop5 = chrono::high_resolution_clock::now();
    auto duration5 = chrono::duration_cast<chrono::microseconds>(stop5 - start5);
    cout << "Two Sort: "
         << duration5.count() << " microseconds" << endl;
    r2100 += num2Reads1;
    w2100 += num2Writes1;

    writeToFile(100, rB100, wB100, rS100, wS100, rQ100, wQ100, rH100, wH100, r2100, w2100);

    // sorting a vector of size 200
    vector<Course> v2a, v2b, v2c, v2d, v2e;
    for (int i = 0; i < 200; ++i){
        v2a.push_back(course.at(i));
        v2b.push_back(course.at(i));
        v2c.push_back(course.at(i));
        v2d.push_back(course.at(i));
        v2e.push_back(course.at(i));
    }
    auto bubbleSortRecord2 = bubbleSort(v2a);
    int rB200 = get<0>(bubbleSortRecord2);
    int wB200 = get<1>(bubbleSortRecord2);
    auto selecSortRecord2 = selectionSort(v2b);
    int rS200 = get<0>(selecSortRecord2);
    int wS200 = get<1>(selecSortRecord2);
    auto quickSortRecord2 = quicksort(v2c);
    int rQ200 = get<0>(quickSortRecord2);
    int wQ200 = get<1>(quickSortRecord2);
    auto heapSortRecord2 = heapSort(v2d);
    int rH200 = get<0>(heapSortRecord2);
    int wH200 = get<1>(heapSortRecord2);

    // 2 sort -- selection sort and then bubble sort on RowId
    auto TwoSortRecord2 = selectionSort(v2e);
    int r2200 = get<0>(TwoSortRecord2);
    int w2200 = get<1>(TwoSortRecord2);
    int num2Reads2 = 0;
    int num2Writes2 = 0;
    for (int j = 0; j < v2e.size(); ++j) {
        for (int i = 0; i < v2e.size() - 1; ++i) {
            num2Reads2 += 2;
            if (v2e[i].getRowId() > v2e[i + 1].getRowId()) {
                // The two elements are out of order. Swap them.
                num2Writes2 += 3;
                num2Writes2 += 3;
                Course temp = v2e[i];
                v2e[i] = v2e[i + 1];
                v2e[i + 1] = temp;
            }
        }
    }
    r2200 += num2Reads2;
    w2200 += num2Writes2;

    writeToFile(200, rB200, wB200, rS200, wS200, rQ200, wQ200, rH200, wH200, r2200, w2200);

    // sorting a vector of size 300
    vector<Course> v3a, v3b, v3c, v3d, v3e;
    for (int i = 0; i < 300; ++i){
        v3a.push_back(course.at(i));
        v3b.push_back(course.at(i));
        v3c.push_back(course.at(i));
        v3d.push_back(course.at(i));
        v3e.push_back(course.at(i));
    }
    auto bubbleSortRecord3 = bubbleSort(v3a);
    int rB300 = get<0>(bubbleSortRecord3);
    int wB300 = get<1>(bubbleSortRecord3);
    auto selecSortRecord3 = selectionSort(v3b);
    int rS300 = get<0>(selecSortRecord3);
    int wS300 = get<1>(selecSortRecord3);
    auto quickSortRecord3 = quicksort(v3c);
    int rQ300 = get<0>(quickSortRecord3);
    int wQ300 = get<1>(quickSortRecord3);
    auto heapSortRecord3 = heapSort(v3d);
    int rH300 = get<0>(heapSortRecord3);
    int wH300 = get<1>(heapSortRecord3);

    // 2 sort -- selection sort and then bubble sort on RowId
    auto TwoSortRecord3 = selectionSort(v3e);
    int r2300 = get<0>(TwoSortRecord3);
    int w2300 = get<1>(TwoSortRecord3);
    int num2Reads3 = 0;
    int num2Writes3 = 0;
    for (int j = 0; j < v3e.size(); ++j) {
        for (int i = 0; i < v3e.size() - 1; ++i) {
            num2Reads3 += 2;
            if (v3e[i].getRowId() > v3e[i + 1].getRowId()) {
                // The two elements are out of order. Swap them.
                num2Writes3 += 3;
                num2Reads3 += 3;
                Course temp = v3e[i];
                v3e[i] = v3e[i + 1];
                v3e[i + 1] = temp;
            }
        }
    }
    r2300 += num2Reads3;
    w2300 += num2Writes3;

    writeToFile(300, rB300, wB300, rS300, wS300, rQ300, wQ300, rH300, wH300, r2300, w2300);

    // sorting a vector of size 400
    vector<Course> v4a, v4b, v4c, v4d, v4e;
    for (int i = 0; i < 400; ++i){
        v4a.push_back(course.at(i));
        v4b.push_back(course.at(i));
        v4c.push_back(course.at(i));
        v4d.push_back(course.at(i));
        v4e.push_back(course.at(i));
    }
    auto bubbleSortRecord4 = bubbleSort(v4a);
    int rB400 = get<0>(bubbleSortRecord4);
    int wB400 = get<1>(bubbleSortRecord4);
    auto selecSortRecord4 = selectionSort(v4b);
    int rS400 = get<0>(selecSortRecord4);
    int wS400 = get<1>(selecSortRecord4);
    auto quickSortRecord4 = quicksort(v4c);
    int rQ400 = get<0>(quickSortRecord4);
    int wQ400 = get<1>(quickSortRecord4);
    auto heapSortRecord4 = heapSort(v4d);
    int rH400 = get<0>(heapSortRecord4);
    int wH400 = get<1>(heapSortRecord4);

    // 2 sort -- selection sort and then bubble sort on RowId
    auto TwoSortRecord4 = selectionSort(v4e);
    int r2400 = get<0>(TwoSortRecord4);
    int w2400 = get<1>(TwoSortRecord4);
    int num2Reads4 = 0;
    int num2Writes4 = 0;
    for (int j = 0; j < v4e.size(); ++j) {
        for (int i = 0; i < v4e.size() - 1; ++i) {
            num2Reads4 += 2;
            if (v4e[i].getRowId() > v4e[i + 1].getRowId()) {
                // The two elements are out of order. Swap them.
                num2Writes4 += 3;
                num2Reads4 += 3;
                Course temp = v4e[i];
                v4e[i] = v4e[i + 1];
                v4e[i + 1] = temp;
            }
        }
    }
    r2400 += num2Reads4;
    w2400 += num2Writes4;

    writeToFile(400, rB400, wB400, rS400, wS400, rQ400, wQ400, rH400, wH400, r2400, w2400);

    // sorting a vector of size 500
    vector<Course> v5a, v5b, v5c, v5d, v5e;
    for (int i = 0; i < 500; ++i){
        v5a.push_back(course.at(i));
        v5b.push_back(course.at(i));
        v5c.push_back(course.at(i));
        v5d.push_back(course.at(i));
        v5e.push_back(course.at(i));
    }
    auto bubbleSortRecord5 = bubbleSort(v5a);
    int rB500 = get<0>(bubbleSortRecord5);
    int wB500 = get<1>(bubbleSortRecord5);
    auto selecSortRecord5 = selectionSort(v5b);
    int rS500 = get<0>(selecSortRecord5);
    int wS500 = get<1>(selecSortRecord5);
    auto quickSortRecord5 = quicksort(v5c);
    int rQ500 = get<0>(quickSortRecord5);
    int wQ500 = get<1>(quickSortRecord5);
    auto heapSortRecord5 = heapSort(v5d);
    int rH500 = get<0>(heapSortRecord5);
    int wH500 = get<1>(heapSortRecord5);

    // 2 sort -- selection sort and then bubble sort on RowId
    auto TwoSortRecord5 = selectionSort(v5e);
    int r2500 = get<0>(TwoSortRecord5);
    int w2500 = get<1>(TwoSortRecord5);
    int num2Reads5 = 0;
    int num2Writes5 = 0;
    for (int j = 0; j < v5e.size(); ++j) {
        for (int i = 0; i < v5e.size() - 1; ++i) {
            num2Reads5 += 2;
            if (v5e[i].getRowId() > v5e[i + 1].getRowId()) {
                // The two elements are out of order. Swap them.
                num2Writes5 += 3;
                num2Reads5 += 3;
                Course temp = v5e[i];
                v5e[i] = v5e[i + 1];
                v5e[i + 1] = temp;
            }
        }
    }
    r2500 += num2Reads5;
    w2500 += num2Writes5;

    writeToFile(500, rB500, wB500, rS500, wS500, rQ500, wQ500, rH500, wH500, r2500, w2500);

    // sorting a vector of size 600
    vector<Course> v6a, v6b, v6c, v6d, v6e;
    for (int i = 0; i < 600; ++i){
        v6a.push_back(course.at(i));
        v6b.push_back(course.at(i));
        v6c.push_back(course.at(i));
        v6d.push_back(course.at(i));
        v6e.push_back(course.at(i));
    }
    auto bubbleSortRecord6 = bubbleSort(v6a);
    int rB600 = get<0>(bubbleSortRecord6);
    int wB600 = get<1>(bubbleSortRecord6);
    auto selecSortRecord6 = selectionSort(v6b);
    int rS600 = get<0>(selecSortRecord6);
    int wS600 = get<1>(selecSortRecord6);
    auto quickSortRecord6 = quicksort(v6c);
    int rQ600 = get<0>(quickSortRecord6);
    int wQ600 = get<1>(quickSortRecord6);
    auto heapSortRecord6 = heapSort(v6d);
    int rH600 = get<0>(heapSortRecord6);
    int wH600 = get<1>(heapSortRecord6);

    // 2 sort -- selection sort and then bubble sort on RowId
    auto TwoSortRecord6 = selectionSort(v6e);
    int r2600 = get<0>(TwoSortRecord6);
    int w2600 = get<1>(TwoSortRecord6);
    int num2Reads6 = 0;
    int num2Writes6 = 0;
    for (int j = 0; j < v6e.size(); ++j) {
        for (int i = 0; i < v6e.size() - 1; ++i) {
            num2Reads6 += 2;
            if (v6e[i].getRowId() > v6e[i + 1].getRowId()) {
                // The two elements are out of order. Swap them.
                num2Writes6 += 3;
                num2Reads6 += 3;
                Course temp = v6e[i];
                v6e[i] = v6e[i + 1];
                v6e[i + 1] = temp;
            }
        }
    }
    r2600 += num2Reads6;
    w2600 += num2Writes6;

    writeToFile(600, rB600, wB600, rS600, wS600, rQ600, wQ600, rH600, wH600, r2600, w2600);

    // sorting a vector of size 700
    vector<Course> v7a, v7b, v7c, v7d, v7e;
    for (int i = 0; i < 700; ++i){
        v7a.push_back(course.at(i));
        v7b.push_back(course.at(i));
        v7c.push_back(course.at(i));
        v7d.push_back(course.at(i));
        v7e.push_back(course.at(i));
    }
    auto bubbleSortRecord7 = bubbleSort(v7a);
    int rB700 = get<0>(bubbleSortRecord7);
    int wB700 = get<1>(bubbleSortRecord7);
    auto selecSortRecord7 = selectionSort(v7b);
    int rS700 = get<0>(selecSortRecord7);
    int wS700 = get<1>(selecSortRecord7);
    auto quickSortRecord7 = quicksort(v7c);
    int rQ700 = get<0>(quickSortRecord7);
    int wQ700 = get<1>(quickSortRecord7);
    auto heapSortRecord7 = heapSort(v7d);
    int rH700 = get<0>(heapSortRecord7);
    int wH700 = get<1>(heapSortRecord7);

    // 2 sort -- selection sort and then bubble sort on RowId
    auto TwoSortRecord7 = selectionSort(v7e);
    int r2700 = get<0>(TwoSortRecord7);
    int w2700 = get<1>(TwoSortRecord7);
    int num2Reads7 = 0;
    int num2Writes7 = 0;
    for (int j = 0; j < v7e.size(); ++j) {
        for (int i = 0; i < v7e.size() - 1; ++i) {
            num2Reads7 += 2;
            if (v7e[i].getRowId() > v7e[i + 1].getRowId()) {
                // The two elements are out of order. Swap them.
                num2Writes7 += 3;
                num2Reads7 += 3;
                Course temp = v7e[i];
                v7e[i] = v7e[i + 1];
                v7e[i + 1] = temp;
            }
        }
    }
    r2700 += num2Reads7;
    w2700 += num2Writes7;

    writeToFile(700, rB700, wB700, rS700, wS700, rQ700, wQ700, rH700, wH700, r2700, w2700);

    // sorting a vector of size 800
    vector<Course> v8a, v8b, v8c, v8d, v8e;
    for (int i = 0; i < 800; ++i){
        v8a.push_back(course.at(i));
        v8b.push_back(course.at(i));
        v8c.push_back(course.at(i));
        v8d.push_back(course.at(i));
        v8e.push_back(course.at(i));
    }
    auto bubbleSortRecord8 = bubbleSort(v8a);
    int rB800 = get<0>(bubbleSortRecord8);
    int wB800 = get<1>(bubbleSortRecord8);
    auto selecSortRecord8 = selectionSort(v8b);
    int rS800 = get<0>(selecSortRecord8);
    int wS800 = get<1>(selecSortRecord8);
    auto quickSortRecord8 = quicksort(v8c);
    int rQ800 = get<0>(quickSortRecord8);
    int wQ800 = get<1>(quickSortRecord8);
    auto heapSortRecord8 = heapSort(v8d);
    int rH800 = get<0>(heapSortRecord8);
    int wH800 = get<1>(heapSortRecord8);

    // 2 sort -- selection sort and then bubble sort on RowId
    auto TwoSortRecord8 = selectionSort(v8e);
    int r2800 = get<0>(TwoSortRecord8);
    int w2800 = get<1>(TwoSortRecord8);
    int num2Reads8 = 0;
    int num2Writes8 = 0;
    for (int j = 0; j < v8e.size(); ++j) {
        for (int i = 0; i < v8e.size() - 1; ++i) {
            num2Reads8 += 2;
            if (v8e[i].getRowId() > v8e[i + 1].getRowId()) {
                // The two elements are out of order. Swap them.
                num2Writes8 += 3;
                num2Reads8 += 3;
                Course temp = v8e[i];
                v8e[i] = v8e[i + 1];
                v8e[i + 1] = temp;
            }
        }
    }
    r2800 += num2Reads8;
    w2800 += num2Writes8;

    writeToFile(800, rB800, wB800, rS800, wS800, rQ800, wQ800, rH800, wH800, r2800, w2800);

    // sorting a vector of size 900
    vector<Course> v9a, v9b, v9c,v9d, v9e;
    for (int i = 0; i < 900; ++i){
        v9a.push_back(course.at(i));
        v9b.push_back(course.at(i));
        v9c.push_back(course.at(i));
        v9d.push_back(course.at(i));
        v9e.push_back(course.at(i));
    }
    auto bubbleSortRecord9 = bubbleSort(v9a);
    int rB900 = get<0>(bubbleSortRecord9);
    int wB900 = get<1>(bubbleSortRecord9);
    auto selecSortRecord9 = selectionSort(v9b);
    int rS900 = get<0>(selecSortRecord9);
    int wS900 = get<1>(selecSortRecord9);
    auto quickSortRecord9 = quicksort(v9c);
    int rQ900 = get<0>(quickSortRecord9);
    int wQ900 = get<1>(quickSortRecord9);
    auto heapSortRecord9 = heapSort(v9d);
    int rH900 = get<0>(heapSortRecord9);
    int wH900 = get<1>(heapSortRecord9);

    // 2 sort -- selection sort and then bubble sort on RowId
    auto TwoSortRecord9 = selectionSort(v9e);
    int r2900 = get<0>(TwoSortRecord9);
    int w2900 = get<1>(TwoSortRecord9);
    int num2Reads9 = 0;
    int num2Writes9 = 0;
    for (int j = 0; j < v9e.size(); ++j) {
        for (int i = 0; i < v9e.size() - 1; ++i) {
            num2Reads9 += 2;
            if (v9e[i].getRowId() > v9e[i + 1].getRowId()) {
                // The two elements are out of order. Swap them.
                num2Writes9 += 3;
                num2Reads9 += 3;
                Course temp = v9e[i];
                v9e[i] = v9e[i + 1];
                v9e[i + 1] = temp;
            }
        }
    }
    r2900 += num2Reads9;
    w2900 += num2Writes9;

    writeToFile(900, rB900, wB900, rS900, wS900, rQ900, wQ900, rH900, wH900, r2900, w2900);

    // sorting a vector of size 1000
    vector<Course> v10a, v10b, v10c, v10d, v10e;
    for (int i = 0; i < 1000; ++i){
        v10a.push_back(course.at(i));
        v10b.push_back(course.at(i));
        v10c.push_back(course.at(i));
        v10d.push_back(course.at(i));
        v10e.push_back(course.at(i));
    }
    auto bubbleSortRecord10 = bubbleSort(v10a);
    int rB1000 = get<0>(bubbleSortRecord10);
    int wB1000 = get<1>(bubbleSortRecord10);
    auto selecSortRecord10 = selectionSort(v10b);
    int rS1000 = get<0>(selecSortRecord10);
    int wS1000 = get<1>(selecSortRecord10);
    auto quickSortRecord10 = quicksort(v10c);
    int rQ1000 = get<0>(quickSortRecord10);
    int wQ1000 = get<1>(quickSortRecord10);
    auto heapSortRecord10 = heapSort(v10d);
    int rH1000 = get<0>(heapSortRecord10);
    int wH1000 = get<1>(heapSortRecord10);

    // 2 sort -- selection sort and then bubble sort on RowId
    auto TwoSortRecord10 = selectionSort(v10e);
    int r21000 = get<0>(TwoSortRecord10);
    int w21000 = get<1>(TwoSortRecord10);
    int num2Reads10 = 0;
    int num2Writes10 = 0;
    for (int j = 0; j < v10e.size(); ++j) {
        for (int i = 0; i < v10e.size() - 1; ++i) {
            num2Reads10 += 2;
            if (v10e[i].getRowId() > v10e[i + 1].getRowId()) {
                // The two elements are out of order. Swap them.
                num2Writes10 += 3;
                num2Reads10 += 3;
                Course temp = v10e[i];
                v10e[i] = v10e[i + 1];
                v10e[i + 1] = temp;
            }
        }
    }
    r21000 += num2Reads10;
    w21000 += num2Writes10;

    writeToFile(1000, rB1000, wB1000, rS1000, wS1000, rQ1000, wQ1000, rH1000, wH1000, r21000, w21000);

    return 0;


    }