#include <iostream>
#include <cstdlib>
#include <fstream>

using namespace std;


int main(){
    ifstream infile;

    infile.open("GSE64881_segmentation_at_30000bp.passqc.multibam.txt");

    if(!infile){
        cout << "Error Opening File\n";
        exit(0);
    }

    //1.) 

    string str = "";
    int linecount = -1;
    while(getline(infile,str)){
        linecount++;
    }
    infile.close();
    //cout << "Line Count: " << linecount << endl; //answer 1

    //2.)

    infile.open("GSE64881_segmentation_at_30000bp.passqc.multibam.txt");

    int num_nps = -2;

    string nps = "";
    getline(infile, nps);
    // cout << "nps" << nps << '\n';
    for(int i = 0; i < nps.size(); i++){
        if(nps.at(i) == '\t'){
            num_nps++;
        }
    }

    //cout << "Numer of nps" << num_nps << '\n'; //answer 2

    infile.close();


    //3.)

    infile.open("GSE64881_segmentation_at_30000bp.passqc.multibam.txt");


    string line = "";
    string trash;
    int num_ones = 0;

    getline(infile, trash);
    for(int i = 0; i < linecount; i++){
        infile >> trash;
        infile >> trash;
        infile >> trash;
        while(cin >> trash){
            if(trash == "1"){
                num_ones++;
            }
        }    
    }

    cout << "num_ones: " << num_ones << endl;





    infile.close();

    return 0;
}

