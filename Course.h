//
// Created by Eyu on 11/7/2021.
//

#ifndef PROJECT4_COURSE_H
#define PROJECT4_COURSE_H

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <string>
#include <vector>
#include <array>

/**
 * Class referring to the fall courses of uvm in 2021,
 * resource from https://serval.uvm.edu/~rgweb/batch/curr_enroll_fall.html
 *
 * Fields (naming follows abbreviations used in CSV file):
 *     rowId: the row ID for ordering
 *     subject: the subject that the course belongs to
 *     title: the title of the course
 *     CRN: the registration number for this course
 *     college: the college code that the class belongs to
 *     maxEnrollment: the maximum enrollment of the course
 *     currEnrollment: the current enrollment of this course by now
 *     instructor: the professor name of the course, in the format of FirstName, LastName
 *     email: the uvm email address of the professor
 */

class Course {

private:
    std::string subject, title, college, instructor, email;
    int CRN, maxEnrollment, currEnrollment, rowId;

public:
    // Constructor with no input param
    Course() {
    }

    // constructor with input param
    Course(int rowId, int CRN, int maxEnrollment, int currEnrollment,
           std::string subject,
           std::string title,
           std::string college,
           std::string instructor,
           std::string email
    ) {
        this->rowId = rowId;
        this->CRN = CRN;
        this->maxEnrollment = maxEnrollment;
        this->currEnrollment = currEnrollment;
        this->subject = subject;
        this->title = title;
        this->college = college;
        this->instructor = instructor;
        this->email = email;
    }

    // Getters
    int getCRN() const {
        return CRN;
    }

    int getMaxEnrollment() const {
        return maxEnrollment;
    }

    int getCurrEnrollment() const {
        return currEnrollment;
    }

    int getRowId() const {
        return rowId;
    }

    std::string getSubject() const {
        return subject;
    }

    std::string getTitle() const {
        return title;
    }

    std::string getCollege() const {
        return college;
    }

    std::string getInstructor() const {
        return instructor;
    }

    std::string getEmail() const {
        return email;
    }

    // Setters
    void setCRN(int CRN) {
        this->CRN = CRN;
    }

    void setMaxEnrollment(int maxEnrollment) {
        this->maxEnrollment = maxEnrollment;
    }

    void setCurrEnrollment(int currEnrollment) {
        this->currEnrollment = currEnrollment;
    }

    void setRowId(int rowId) {
        this->rowId = rowId;
    }

    void setSubject(std::string subject) {
        this->subject = subject;
    }

    void setTitle(std::string title) {
        this->title = title;
    }

    void setCollege(std::string college) {
        this->college = college;
    }

    void setInstructor(std::string instructor) {
        this->instructor = instructor;
    }

    void setEmail(std::string email) {
        this->email = email;
    }

    /**
     * Overload stream insertion operator for pretty printing
     */
    friend std::ostream& operator << (std::ostream& outs, const Course& cObject) {
        outs << std::fixed
             << std::left
             << std::setw(7) << cObject.getRowId()
             << std::setw(8) << cObject.getSubject()
             << std::setw(40) << cObject.getTitle()
             << std::setw(7) << cObject.getCRN()
             << std::setw(9) << cObject.getCollege()
             << std::setw(5) << cObject.getMaxEnrollment()
             << std::setw(5) << cObject.getCurrEnrollment()
             << std::setw(25) << cObject.getInstructor()
             << std::setw(25) << cObject.getEmail();
        return outs;
    }

    /**
     * Overload < operator for object comparison
     */
    friend bool operator < (const Course& lhs, const Course& rhs) {
        return lhs.getMaxEnrollment() < rhs.getMaxEnrollment();
    }

    /**
     * Overload > operator for object comparison
     */
    friend bool operator > (const Course& lhs, const Course& rhs) {
        return lhs.getMaxEnrollment() > rhs.getMaxEnrollment();
    }

    /**
     * Overload <= operator for object comparison
     */
    friend bool operator <= (const Course& lhs, const Course& rhs) {
        return lhs.getMaxEnrollment() <= rhs.getMaxEnrollment();
    }

    /**
     * Overload >= operator for object comparison
     */
    friend bool operator >= (const Course& lhs, const Course& rhs) {
        return lhs.getMaxEnrollment() >= rhs.getMaxEnrollment();
    }

    /**
     * Overload == operator for object comparison
     */
    friend bool operator == (const Course& lhs, const Course& rhs) {
        return lhs.getMaxEnrollment() == rhs.getMaxEnrollment();
    }

    /**
     * Overload != operator for object comparison
     */
    friend bool operator != (const Course& lhs, const Course& rhs) {
        return lhs.getMaxEnrollment() != rhs.getMaxEnrollment();
    }

};

// Other functions

/**
 * Read data into a vector from the downloaded CSV file.
 * Returns true if everything goes well; false otherwise.
 */
bool loadFromFile(std::string file, std::vector<Course>& objs) {
    std::ifstream fileInput;
    fileInput.open(file);

    if (!fileInput) {
        // Bad file / could not open
        return false;
    }

    // Opened file OK
    // Get rid of the header line: read it and throw it away
    std::string header = "";
    getline(fileInput, header);

    // Read in data while the file is still good
    // and we haven't yet reached the end of file (EOF).
    int counter = 1;
    while (fileInput && fileInput.peek() != EOF) {
        // We initialize variables to sentinel values --
        // values we *know* will not appear in our data.
        // If any object fields have one of these values, we
        // *know* it was incorrectly set.
        int rowId = -1, CRN = 90000;
        int maxEnrollment = 0, currEnrollment = -1;
        std::string subject = "INVALID";
        std::string title = "INVALID";
        std::string college = "INVALID";
        std::string instructor = "INVALID";
        std::string email = "INVALID";
        char comma = ',', doubleQuote = '"';

        // Using the >> (extraction) operator to read an int/double
        // Using getline() to read string
        fileInput >> rowId;
        if (rowId != counter) {
            std::cout << rowId << " != " << counter << std::endl;
            return false;
        }
        fileInput >> comma;
        getline(fileInput, subject, comma);
        getline(fileInput, title, comma);
        fileInput >> CRN;
        fileInput >> comma;
        getline(fileInput, college, comma);
        fileInput >> maxEnrollment;
        fileInput >> comma;
        fileInput >> currEnrollment;
        fileInput >> comma;
        //to read instructor name with double quotes using peek()
        if (fileInput.peek() == doubleQuote) {
            fileInput >> doubleQuote;
            getline(fileInput, instructor, doubleQuote);
            fileInput >> comma;
        }else {getline(fileInput, instructor, comma);}
        //in case of no instructor assigned to this course, the email would be empty
        if (fileInput.peek() != '\r' && fileInput.peek() != '\n') {
            getline(fileInput, email);
        }

        //create a Course object using the dta read above
        Course courseObj(rowId, CRN, maxEnrollment, currEnrollment,
                         subject, title, college, instructor, email);
        objs.push_back(courseObj);
        counter += 1;
    }
    fileInput.close();
    return true;
}

/**
 * Count the number of courses without instructor and email
 * No courses with instructor existing but no email
 */
void needInstructor(std::vector<Course>& objs) {
    int numNeedInstructor = 0;
    int numWithInstructor = 0;

    for (Course& cour : objs) {
        std::string instructorInfo = cour.getInstructor();
        std::string emailInfo = cour.getEmail();
        if (instructorInfo.compare("Staff, . ") == 0) {
            numNeedInstructor += 1;
        }
        else if (instructorInfo.empty()) {numNeedInstructor += 1; }
        else {numWithInstructor += 1; }
    }

    std::cout << std::endl;
    std::cout << std::right << "Courses need Instructor: " << std::setw(5) << numNeedInstructor << std::endl;
    std::cout << std::right << "Courses All Set: " << std::setw(5) << numWithInstructor << std::endl;
    std::cout << std::string(95, '=') << std::endl;
}

/**
 * To validate the data
 * Possible issue: incorrect row ID, empty fields(emails could be empty as no instructor assigned),
 * some instructor with no existing email
 * The total count of each issue would be output
 */
void validFields(std::vector<Course>& objs) {
    int i = 0;
    int badRowId = 0, badCRN = 0, badMaxEnrollment = 0, badCurrEnrollment = 0, badSubject = 0;
    int badTitle = 0, badCollege = 0, badInstructor = 0, badEmail = 0;
    for (Course& cour : objs) {
        // Row ID should always = i + 1
        if (cour.getRowId() != i + 1) {
            badRowId += 1;
        }
        if (cour.getCRN() < 90000) {
            badCRN += 1;
        }
        if (cour.getMaxEnrollment() == 0) {
            badMaxEnrollment += 1;
        }
        if (cour.getCurrEnrollment() < 0) {
            badCurrEnrollment += 1;
        }

        //check the string fields are good
        std::string subjectValid = cour.getSubject();
        if (subjectValid.empty()) {
            badSubject += 1;
        }
        std::string titleValid = cour.getTitle();
        if (titleValid.empty()) {
            badTitle += 1;
        }
        std::string collegeValid = cour.getCollege();
        if (collegeValid.empty()) {
            badCollege += 1;
        }
        //check if there are existing teachers forget their email
        std::string instructorValid = cour.getInstructor();
        std::string emailValid = cour.getEmail();
        if (instructorValid.empty()) {
            badInstructor += 1;
        }
        else if (!instructorValid.empty() && emailValid.empty()) {
            badEmail += 1;
        }
        i += 1;
    }
    std::cout << std::string(95, '=') << std::endl;
    std::cout << std::endl;
    if (badRowId > 0) {std::cout << "There are " << badRowId << " bad row ID. " << std::endl;}
    if (badCRN > 0) {std::cout << "There are " << badCRN << " bad CRN. " << std::endl;}
    if (badMaxEnrollment > 0) {std::cout << "There are " << badMaxEnrollment << " bad max enrollment data. " << std::endl;}
    if (badCurrEnrollment > 0) {std::cout << "There are " << badCurrEnrollment << " bad current enrollment data. " << std::endl; }
    if (badSubject > 0) {std::cout << "There are " << badSubject << " bad subject data. " << std::endl;}
    if (badCollege > 0) {std::cout << "There are " << badCollege << " bad college data. " << std::endl;}
    if (badTitle > 0) {std::cout << "There are " << badTitle << " bad title. " << std::endl;}
    if (badInstructor > 0) {std::cout << "There are " << badInstructor << " bad instructor data. " << std::endl;}
    if (badEmail > 0) {std::cout << "There are " << badEmail << " bad email data. " << std::endl;}
    std::cout << std::string(95, '=') << std::endl;
}


/**
 * Assertions are cool. You can use them to test your objects and functions.
 * You don't have to use them, but they make testing your code easy!
 * If an assertion fails at runtime, your program will raise an exception.
 */
void CourseTest() {
    Course testObject = Course();
    // Test setters and getters
    testObject.setRowId(1);
    assert (testObject.getRowId() == 1);
    testObject.setCRN(10000);
    assert (testObject.getCRN() == 10000);
    testObject.setCollege("CALS");
    assert (testObject.getCollege() == "CALS");
    testObject.setMaxEnrollment(100);
    assert (testObject.getMaxEnrollment() == 100);
    testObject.setCurrEnrollment(50);
    assert (testObject.getCurrEnrollment() == 50);
    testObject.setSubject("ABIO");
    assert (testObject.getSubject() == "ABIO");
    testObject.setTitle("Doctoral Dissertation Research");
    assert (testObject.getTitle() == "Doctoral Dissertation Research");
    testObject.setInstructor("Harvey, Jean Ruth");
    assert (testObject.getInstructor() == "Harvey, Jean Ruth");
    testObject.setEmail("Jean.Harvey@uvm.edu");
    assert (testObject.getEmail() == "Jean.Harvey@uvm.edu");

    // Test initialization (and getters)
    Course testObject2 = Course(3, 96196, 20, 0, "ABIO", "Doctoral Dissertation Research",
                                "CALS", "Harvey, Jean Ruth", "Jean.Harvey@uvm.edu");
    assert (testObject2.getRowId() == 3);
    assert (testObject2.getCRN() == 96196);
    assert (testObject2.getMaxEnrollment() == 20);
    assert (testObject2.getCurrEnrollment() == 0);
    assert (testObject2.getSubject() == "ABIO");
    assert (testObject2.getTitle() == "Doctoral Dissertation Research");
    assert (testObject2.getCollege() == "CALS");
    assert (testObject2.getInstructor() == "Harvey, Jean Ruth");
    assert (testObject2.getEmail() == "Jean.Harvey@uvm.edu");
}

#endif //PROJECT4_COURSE_H
