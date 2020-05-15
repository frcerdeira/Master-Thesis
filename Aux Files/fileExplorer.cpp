#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <experimental/filesystem>

using namespace std;
namespace fs = std::experimental::filesystem;

std::string path = fs::current_path(); // get current path

inline std::vector<std::string> getFileNames(){

    std::vector<std::string> names; // create vector with file names

    // Populate vector "names" with file names
    for (const auto & entry : fs::recursive_directory_iterator(path)){
        names.push_back(entry.path());
    }

    return names; // return vector
}

inline std::vector<std::string> filterNames(std::vector<std::string>& strings, std::string pattern){

    std::vector<std::string> namesFiltered;

    for(auto& s1:strings){
        if(s1.find(pattern) != std::string::npos){
            namesFiltered.push_back(s1);
        }
    }

    return namesFiltered;

}

inline std::string eraseFileExtension(std::string fileName){

    std::string newName;
    int pos = 0;

    pos = fileName.find(".");
    newName = fileName.erase(pos);

    return newName;
    
}

void mergeFiles(std::vector<std::string>& fileNames, std::string fileName){

    std::string line;
    std::string nameNewFile;

    nameNewFile = eraseFileExtension(fileName); // file name provided by user to what we will do shit + file extension
    nameNewFile += "_MERGED.dat"; // specifies the folder in the new file name and that an average was performed

    cout << "Merging following files:" << endl;
    ofstream merge;
        merge.open(nameNewFile);

    for(auto& openFile : fileNames){

        cout << openFile << endl;

        std::ifstream file;
        file.open(openFile); // opens file with the name given by the vector 'fileNames'

        if(!file.is_open()){ // throws error if it can't open the file
            perror("Error open");
            exit(EXIT_FAILURE);
        }
        while(getline(file, line)){ // goes from line to line of a file
            merge << line << endl;
        }
    }

    merge.close();

    std::cout << "The file '" << nameNewFile << "' was created on '" << path << "'" << endl;

}

inline std::string readFiles(std::vector<std::string>& fileNames, std::string fileName, std::string folderName){

    std::string line;
    std::string nameNewFile;
    std::string col1String; // vector for column 1
    std::string col2String; // vector for column 2
    std::string col3String; // vector for column 2
    std::vector<double> col1Double; // vector for column 1
    std::vector<double> col2Double; // vector for column 2
    // std::vector<double> col3Double; // vector for column 3

    int sizeOfVector = 0; // so I know how long my for loop needs to be to read all #

    int numberOfFiles = fileNames.size(); // Gives the number of files (Important for averages)

    int i = 0;
    int firstFile = 0; // checks if I'm reading the first file

    nameNewFile = eraseFileExtension(fileName); // file name provided by user to what we will do shit + file extension
    nameNewFile += "_" + folderName + "_AVG.dat"; // specifies the folder in the new file name and that an average was performed

    for(auto& openFile : fileNames){

        i = 0;

        std::ifstream file;
        file.open(openFile); // opens file with the name given by the vector 'fileNames'

        if(!file.is_open()){ // throws error if it can't open the file
            perror("Error open");
            exit(EXIT_FAILURE);
        }
        while(getline(file, line)){ // goes from line to line of a file
            // test << line << endl;
            stringstream ss(line);
            ss >> col1String >> col2String >> col3String; 

            // !! SEG FAULTS HERE !! //
            if(firstFile == 0){
                col1Double.push_back(std::stod(col1String));
                col2Double.push_back(std::stod(col2String));
                // col3Double.push_back(std::stod(col3String));
            }
            else if(firstFile == 1){
                col1Double[i] += std::stod(col1String);
                col2Double[i] += std::stod(col2String);
                // col3Double[i] += std::stod(col3String);
            }

            i++;

        }

        firstFile = 1;

        file.close();
    }

    ofstream test;
        test.open(nameNewFile);

    sizeOfVector = col1Double.size();

    for(i = 0; i < sizeOfVector; i++){
        test << col1Double[i]/numberOfFiles << "\t" << col2Double[i]/numberOfFiles << endl;//"\t" << col3Double[i]/numberOfFiles << endl; 
    }

    test.close();

    std::cout << "The file '" << nameNewFile << "' was created on '" << path << "'" << endl; 

    return nameNewFile; 
}

int main()
{
    std::vector<std::string> names; // vector with the names of all the content
    std::vector<std::string> namesFiltered; // vector with the filtered names
    std::vector<std::string> createdFileNames; // vector with the filtered names
    std::string folderName;
    double value = 0;
    char valueAsString[80];

    std::string inputParameterName; // user answer
    std::string inputFileName; // user answer
    double inputStart = 0;
    double inputEnd = 0;
    double inputIncrement = 0;
    char ANSuserError; // user answer (proceed or not) (user gave wrong file/folder)
    char ANSrunAgain; // user answer (proceed or not) (user wants to go again)
    char ANSmergeFiles; // user answer (proceed or not) (user wants to merge files)

    
    names = getFileNames(); // get all files in folder
    do{
        // Ask user what folder does he want (each folder uses different parameters so we don't want to mix them)
        std::cout << "\nChoose Parameter: ";
        std::cin >> inputParameterName;
        std::cout << "\nStarts at: ";
        std::cin >> inputStart;
        std::cout << "Ends at: ";
        std::cin >> inputEnd;
        std::cout << "Increment: ";
        std::cin >> inputIncrement;

        // Check if user is happy with the files that will be handled
        cout << "\nConfirm (Y/N)? ";
        cin >> ANSuserError;

    } while((ANSuserError == 'N')||(ANSuserError == 'n'));
    
    do{
        // Ask user what is the common name of the files he wants to open
        std::cout << "Choose File Name: ";
        std::cin >> inputFileName;
        inputFileName += ".dat";

        value = inputStart;

        while(value <= inputEnd + 0.001){ // 0.001 to make sure there aren't any rounding errors

            sprintf(valueAsString, "%.2f", value);

            folderName = inputParameterName + "_=" + valueAsString;

            namesFiltered = filterNames(names, folderName); // filter files - specific folder

            namesFiltered = filterNames(namesFiltered, inputFileName); // filter files - specific file name

            // Show user the selected files
            std::cout << "\n\nThe following files are being handled:" << endl;  
            for (auto& filesSelected : namesFiltered){
                std::cout << filesSelected << endl;
            }

            createdFileNames.push_back(readFiles(namesFiltered, inputFileName, folderName));

            value += inputIncrement;
        }
            
        // Check if user is happy with the files that will be handled
        cout << "\n\nDo you want to run the code again for another file (Y/N)? ";
        cin >> ANSrunAgain;

    } while((ANSrunAgain == 'Y')||(ANSrunAgain == 'y'));

    cout << "\n\nDo you want to merge the files that were created (Y/N)? ";
    cin >> ANSmergeFiles;

    if((ANSmergeFiles == 'Y')||(ANSmergeFiles == 'y')){
        mergeFiles(createdFileNames, inputFileName);
    }

   return 0; 

}

