////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function that creates an ROH matrix and a Autozygosity Matrix                                    ////
////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <sstream>
#include <vector>
#include <cmath>
#include <ctime>
// Class to store index start and stop site ##
class ROH_Index
{
private:
    int Chromosome;             /* Chromsome */
    int StartPosition;          /* Start Position nucleotides */
    int EndPosition;            /* End Position Nucleotides */
    int StartIndex;             /* Start Index */
    int EndIndex;               /* End Index */
    int NumberSNP;              /* Number of SNP in ROH */
public:
    ROH_Index();
    ROH_Index(int chr = 0, int stpos = 0, int enpos = 0, int stind = 0, int enind = 0, int numsnp = 0);
    ~ROH_Index();
    int getChr(){return Chromosome;}
    int getStPos(){return StartPosition;}
    int getEnPos(){return EndPosition;}
    int getStInd(){return StartIndex;}
    int getEnInd(){return EndIndex;}
    int getNumSNP(){return NumberSNP;}
};

// Constructors ROH_Index
ROH_Index::ROH_Index()
{
    Chromosome = 0;  StartPosition = 0; EndPosition = 0; StartIndex = 0; EndIndex = 0; NumberSNP = 0;
}
ROH_Index::ROH_Index(int chr, int stpos, int enpos, int stind, int enind, int numsnp)
{
    Chromosome = chr;  StartPosition = stpos; EndPosition = enpos;
    StartIndex = stind; EndIndex = enind; NumberSNP = numsnp;
}
// Destructors
ROH_Index::~ROH_Index(){}

////Summary of ROH output////////
class Ani_ROH
{
private:
    std::string ROHChrome;      /* Chromsomes for each ROH */
    std::string ROHStart;               /* ROH Start in Base number*/
    std::string ROHEnd;                 /* ROH End in Base number */
    std::string Distance;               /* Length of ROHs in megabases*/
public:
    Ani_ROH();
    Ani_ROH(std::string rohchrome = "", std::string ws = "", std::string we = "", std::string dist_MB = "");
    ~Ani_ROH();
    std::string getrohchromeO(){return ROHChrome;}
    std::string getwsO(){return ROHStart;}
    std::string getweO(){return ROHEnd;}
    std::string getdist_MBO(){return Distance;}
};

using namespace std;
int main(int argc, char* argv[])
{
    time_t full_begin_time = time(0);
    cout<<"\n#############################################################\n";
    cout<<"###################################################   #######\n";
    cout<<"###############################################   /~|   #####\n";
    cout<<"############################################   _- `~~~', ####\n";
    cout<<"##########################################  _-~       )  ####\n";
    cout<<"#######################################  _-~          |  ####\n";
    cout<<"####################################  _-~            ;  #####\n";
    cout<<"##########################  __---___-~              |   #####\n";
    cout<<"#######################   _~   ,,                  ;  `,,  ##\n";
    cout<<"#####################  _-~    ;'                  |  ,'  ; ##\n";
    cout<<"###################  _~      '                    `~'   ; ###\n";
    cout<<"############   __---;                                 ,' ####\n";
    cout<<"########   __~~  ___                                ,' ######\n";
    cout<<"#####  _-~~   -~~ _          N                    ,' ########\n";
    cout<<"##### `-_         _           C                  ; ##########\n";
    cout<<"#######  ~~----~~~   ;         S                ; ###########\n";
    cout<<"#########  /          ;         U               ; ###########\n";
    cout<<"#######  /             ;                      ; #############\n";
    cout<<"#####  /                `                    ; ##############\n";
    cout<<"###  /                                      ; ###############\n";
    cout<<"#                                            ################\n";
    cout<<"-------------------------------------------------------------\n";
    cout<<"- ROH_CAL A PROGRAM TO OBTAIN AUTOZYGOSITY AND ROH STATUS   -\n";
    cout<<"- FROM SNP DATA                                             -\n";
    cout<<"- Author: Jeremy T. Howard and Ben Lemmer                   -\n";
    cout<<"- Institution: NCSU                                         -\n";
    cout<<"- Date: 8/12/2015                                           -\n";
    cout<<"-------------------------------------------------------------\n";
    ///////////////////////////////////////////////////
    ///       Read in Parameters from a file        ///
    ///////////////////////////////////////////////////
    if(argc != 2){cout << "Program ended due to a parameter file not given!" << endl; exit (EXIT_FAILURE);}
    string paramterfile = argv[1];
    /* Deletes Old Log File */
    fstream checklog; checklog.open("log_file.txt", std::fstream::out | std::fstream::trunc); checklog.close();
    std::ofstream log_file("log_file.txt", std::ios_base::out); /* open log file*/
    log_file << "===============================================" << endl;
    log_file << "==   Parameters passed to the program        ==" << endl;
    log_file << "===============================================" << endl;
    /* Parameters that are specified */
    string mapfile;                     /* map file name; Chr then position with a space delimeter */
    string genofile;                    /* genotype file; ID then string of genotypes with a space delimeter */
    int roh_cutoff;                     /* Mb cutoff for ROH */
    int roh_threshold;                  /* Number times the mean to cutoff if too low of SNP in ROH */
    string RemoveSNP;                   /* Do you want to remove SNP that are not contained within an ROH that passed threshold */
    string Outfilebasename;             /* base name for output file */
    /* read parameters file */
    vector <string> parm;
    string parline;
    ifstream parfile;
    parfile.open(paramterfile.c_str());
    if(parfile.fail()){cout << "Parameter file not found. Check log file." << endl; exit (EXIT_FAILURE);}
    while (getline(parfile,parline)){parm.push_back(parline);} /* Stores in vector and each new line push back to next space */
    int search = 0;
    while(1)
    {
        size_t fnd = parm[search].find("GENOTYPE_FILE:");
        if(fnd!=std::string::npos)
        {
            size_t pos = parm[search].find(": ", 0); parm[search].erase(0, pos + 2);
            parm[search].erase(remove(parm[search].begin(), parm[search].end(), ' '), parm[search].end()); genofile = parm[search];
            log_file << "Genotypes File: \t\t" << "'" << genofile << "'" << endl; break;
        }
        search++;
        if(search >= parm.size()){cout << endl << "Couldn't find 'GENOTYPE_FILE:' variable in parameter file!" << endl; exit (EXIT_FAILURE);}
    }
    search = 0;
    while(1)
    {
        size_t fnd = parm[search].find("MAP_FILE:");
        if(fnd!=std::string::npos)
        {
            size_t pos = parm[search].find(": ", 0); parm[search].erase(0, pos + 2);
            parm[search].erase(remove(parm[search].begin(), parm[search].end(), ' '), parm[search].end()); mapfile = parm[search];
            log_file << "Map File: \t\t\t"<< "'" << mapfile << "'" << endl; break;
        }
        search++;
        if(search >= parm.size()){cout << endl << "Couldn't find 'MAP_FILE:' variable in parameter file!" << endl; exit (EXIT_FAILURE);}
    }
    search = 0;
    while(1)
    {
        size_t fnd = parm[search].find("ROH_CUTOFF:");
        if(fnd!=std::string::npos)
        {
            size_t pos = parm[search].find(":", 0); parm[search].erase(0, pos + 2);
            parm[search].erase(remove(parm[search].begin(), parm[search].end(), ' '), parm[search].end()); roh_cutoff = atoi(parm[search].c_str());
            roh_cutoff = roh_cutoff * 1000000;
            log_file << "ROH threshold: \t\t\t" << "'" << roh_cutoff << "'" << endl; break;
        }
        search++; if(search >= parm.size()){cout << endl << "Couldn't find 'ROH_CUTOFF:' variable in parameter file!" << endl; exit (EXIT_FAILURE);}
    }
    search = 0;
    while(1)
    {
        size_t fnd = parm[search].find("ROH_THRESHOLD:");
        if(fnd!=std::string::npos)
        {
            size_t pos = parm[search].find(":", 0); parm[search].erase(0, pos + 1);
            parm[search].erase(remove(parm[search].begin(), parm[search].end(), ' '), parm[search].end()); roh_threshold = atoi(parm[search].c_str());
            log_file << "ROH Threshold: \t\t\t" << "'" << roh_threshold << "'" << endl; break;
        }
        search++;
        if(search >= parm.size())
        {
            roh_threshold = 2;
            log_file << "ROH Threshold: \t\t\t" << "'" << roh_threshold << "'" << " (Default)" << endl; break;
        }
    }
    search = 0;
    while(1)
    {
        size_t fnd = parm[search].find("REMOVE_SNP:");
        if(fnd!=std::string::npos)
        {
            size_t pos = parm[search].find(":", 0); parm[search].erase(0, pos + 1);
            parm[search].erase(remove(parm[search].begin(), parm[search].end(), ' '), parm[search].end()); RemoveSNP = parm[search];
            log_file << "Remove SNP: \t\t\t" << "'" << RemoveSNP << "'" << endl; break;
        }
        search++;
        if(search >= parm.size())
        {
            RemoveSNP = "yes";
            log_file << "Remove SNP: \t\t\t" << "'" << RemoveSNP << "'" << " (Default)" << endl; break;
        }
    }
    search = 0;
    while(1)
    {
        size_t fnd = parm[search].find("OUT_FILE:");
        if(fnd!=std::string::npos)
        {
            size_t pos = parm[search].find(": ", 0); parm[search].erase(0, pos + 1);
            parm[search].erase(remove(parm[search].begin(), parm[search].end(), ' '), parm[search].end()); Outfilebasename = parm[search]; break;
        }
        search++; if(search >= parm.size()){cout << endl << "Couldn't find 'OUTFILE:' variable in parameter file!" << endl; exit (EXIT_FAILURE);}
    }
    string OutputFileAuto = Outfilebasename + "_Autozygosity";       /* Name of output file Auto */
    string OutputFileROH = Outfilebasename + "_ROH";                 /* name of output file ROH */
    string OutputFileSummary = Outfilebasename + "_Summary";
    log_file << "Autozygosity output file: \t" << "'" << OutputFileAuto << "'" << endl;
    log_file << "ROH output file: \t\t" << "'" << OutputFileROH << "'" << endl;
    log_file << "Summary output file: \t\t" << "'" << OutputFileSummary << "'" << endl;
    log_file << "===============================================\n\n";
    log_file << "====================================\n";
    log_file << "===\tIndexing Positions \t====\n";
    log_file << "====================================" << endl;
    /* Read in map file don't need to know how many SNP are in the file */
    vector <string> numbers;
    /* Import file and put each row into a vector */
    string line;
    ifstream infile;
    infile.open(mapfile.c_str());
    if(infile.fail()){cout << "Error Opening map File\n"; exit (EXIT_FAILURE);}
    while (getline(infile,line))
    {
        numbers.push_back(line);        /* Stores in vector and each new line push back to next space */
    }
    int rows = numbers.size();          /* Determine number of SNP */
    log_file << "   - Total Number of SNP in Map file is " << rows << endl;
    vector < int > chr(rows,0);         /* stores chromosome in vector */
    vector < int > position(rows,0);    /* stores position */
    int index[rows];                    /* used for when grabs size of 4Mb */
    for(int i = 0; i < numbers.size(); i++)
    {
        string temp = numbers[i];                       // grab a line
        size_t pos = temp.find(" ", 0);                 // Find position where delimiter is at
        string tempa = temp.substr(0,pos);              // grab part
        chr[i] = atoi(tempa.c_str());                   // Convert it to an integer
        temp.erase(0, pos+1);                           // Erase phenotype from temp and delimeter
        position[i] = atoi(temp.c_str());               // Convert it to an integer
        index[i] = i;
    }
    numbers.clear();                                    // clear vector that holds each row
    vector<ROH_Index> roh_index;
    /* Create index to grab correct columns from genotype file when constructing ROH and Autozygosity*/
    for(int i = 0; i < rows; i++)
    {
        int j = i;
        while(1)
        {
            if(chr[i] == chr[j])
            {
                int temp = position[j] - position[i];
                if(temp > roh_cutoff)
                {
                    int numsnp = index[j] - index[i] + 1;
                    ROH_Index roh_temp(chr[i],position[i],position[j],index[i],index[j],numsnp);
                    roh_index.push_back(roh_temp);                  /* store in vector of roh_index objects */
                    break;
                }
                if(temp <= roh_cutoff){j++;}
            }
            if(chr[i] != chr[j]){break;}
        }
    }
    log_file << "   - Total number of ROH windows prior to removal " << roh_index.size() << endl;
    /* Figure out mean and Standard Deviation */
    double sum = 0.0;
    for(int i = 0; i < roh_index.size(); i++){sum += roh_index[i].getNumSNP();}
    double mean = sum / roh_index.size();
    double sq_sum;
    for(int i = 0; i < roh_index.size(); i++){sq_sum += roh_index[i].getNumSNP() * roh_index[i].getNumSNP();}
    double stdev = sqrt(sq_sum / roh_index.size() - mean * mean);
    log_file << "   - Mean +/- S.D. SNP size for an ROH: " << mean << " " << stdev << endl;
    int SNPSizeCutoff = mean - roh_threshold * stdev + 0.5;     /* Ensures it rounds up */
    log_file << "   - Any ROH window with SNP size below: " << SNPSizeCutoff << " removed." << endl;
    /* Remove ROH windows that fall below threshold */
    int row = roh_index.size();                             /* Figure out number of rows */
    int i = 0;                                      /* counter to determine where you are at */
    /* Keep ROH that are greater than a given threshold */
    while(i < row)
    {
       while(1)
        {
            if(roh_index[i].getNumSNP() < SNPSizeCutoff)
            {
                roh_index.erase(roh_index.begin()+i);   /* Remove from allfreq */
                row = row - 1;                          /* One less column in G and row in allfreq */
                break;
            }
            else{i++; break;}
        }
    }
    log_file << "   - Total number of ROH windows after removal " << roh_index.size() << endl;
    log_file << "============================================\n";
    log_file << "===\tCalculating ROH & Autozygosity \t====\n";
    log_file << "============================================" << endl;
    /* Output file with a header with index numbers */
    /* Make sure you know if overwriting file       */
    fstream file;
    file.open(OutputFileROH.c_str(), ios_base::out | ios_base::in);  // will not create file
    if (file.is_open())
    {
        log_file << "\nWARNING, output file" <<OutputFileROH  << " already exists and will be overwritten"<< endl;
        fstream check; check.open(OutputFileROH.c_str(), std::fstream::out | std::fstream::trunc); check.close();
    }
    file.close();
    fstream file1;
    file1.open(OutputFileAuto.c_str(), ios_base::out | ios_base::in);
    if (file1.is_open())
    {
        log_file << "WARNING, output file" <<OutputFileAuto  << " already exists and will be overwritten"<< endl;
        fstream check; check.open(OutputFileAuto.c_str(), std::fstream::out | std::fstream::trunc); check.close();
    }
    file1.close();
    fstream file2;
    file2.open(OutputFileSummary.c_str(), ios_base::out | ios_base::in);
    if (file2.is_open())
    {
        log_file << "WARNING, output file" << OutputFileSummary  << " already exists and will be overwritten"<< endl;
        fstream check; check.open(OutputFileSummary.c_str(), std::fstream::out | std::fstream::trunc); check.close();
    }
    file2.close();
    /* Save as a continuous string and then output */
    stringstream outputstring(stringstream::out);
    for(int i = 0; i < roh_index.size(); i++)               /* loop across and fill rohrange array */
    {
        if(i == 0){outputstring << roh_index[i].getChr() << "_" << roh_index[i].getStPos() << "_" << roh_index[i].getEnPos();}
        if(i > 0){outputstring << " " << roh_index[i].getChr() << "_" << roh_index[i].getStPos() << "_" << roh_index[i].getEnPos();}
    }
    outputstring << endl;
    std::ofstream output1(OutputFileROH.c_str(), std::ios_base::app | std::ios_base::out);
    output1 << outputstring.str(); outputstring.str(""); outputstring.clear();
    /* These are used to calculate proportion of genome in an ROH
    /* it will be the same for all individuals so for first animal remove ones that are 5 */
    vector < int > temp_chr; vector < int > temp_pos;
    vector < int > clb; vector < int > cub; double genome_length = 0.0; /* upper and lower bound for each chromosome */
    temp_chr.assign(chr.begin(),chr.end()); temp_pos.assign(position.begin(),position.end());
    int linenum = 1;
    /* Import file */
    ifstream infile1;
    infile1.open(genofile.c_str());
    if(infile1.fail()){cout << "Error Opening File\n";exit (EXIT_FAILURE);}
    while (getline(infile1,line))
    {
        if(linenum == 1){cout << "- Number of lines read: " << endl;}
        string ID;                                              /* Animal ID */
        vector <double> roh_status(roh_index.size(),0);         /* status of roh 0 is not and 1 is in roh */
        size_t pos = line.find(" ", 0);                         /* Find position where delimiter is at */
        ID = line.substr(0,pos);                                /* Place in phenotype vector need to make string to double */
        line.erase(0, pos + 1);                                 /* Erase phenotype from temp and delimeter */
        string geno = line;                                     /* Only thing that is left in temp is genotype string */
        vector < int > genotypes(geno.size(),0);                /* Array to store genotypes */
        vector < int > auto_status(geno.size(),5);              /* Vector to store autozygosity genotypes */
        double homozygosity = 0.0;
        /* Loop through and place in array */
        for(int i = 0; i < geno.size(); i++)
        {
            genotypes[i] = geno[i] - 48;                                    /* ASCI value is 48 for 0 */
            if(genotypes[i] == 3 || genotypes[i] == 4){genotypes[i] = 1;}   /* Convert 3 or 4 to 1 */
            if(genotypes[i] == 2){genotypes[i] = 0;}                        /* Convert 2 to 0 */
            if(genotypes[i] == 0){homozygosity += 1;}
        }
        homozygosity = homozygosity / double(geno.size());
        for(int i = 0; i < roh_index.size(); i++)                           /* loop across roh indexes */
        {
            double sumgeno = 0.0;                                                         /* sum up current roh index */
            /* sum array of genotypes; add one because need to include last snp */
            for(int j = roh_index[i].getStInd(); j < roh_index[i].getEnInd() + 1; j++)
            {
                sumgeno += genotypes[j];
            }
            sumgeno = sumgeno / roh_index[i].getNumSNP();                           /* calculate proportion */
            if(sumgeno > 0){roh_status[i] = 0;}
            if(sumgeno == 0){roh_status[i] = 1;}
            for(int j = roh_index[i].getStInd(); j < roh_index[i].getEnInd() + 1; j++)  /* fill autozygosity matrix */
            {
                if(auto_status[j] == 0 && roh_status[i] == 1)                       /* has been changed to a not in an ROH but is now in one */
                {
                    auto_status[j] = 1;                                             /* SNP is within an ROH */
                }
                if(auto_status[j] == 5 && roh_status[i] == 0)                       /* has not been changed yet and is not in ROH */
                {
                    auto_status[j] = 0;                                             /* SNP is not within an ROH */
                }
                if(auto_status[j] == 5 && roh_status[i] == 1)                       /* has not been changed yet and is in an ROH */
                {
                    auto_status[j] = 1;                                             /* SNP is within an ROH */
                }
                /* once been tagged as being in an ROH can't go back to not being in an ROH */
            }
        }
        /* Now begin to calculate the proportion of genome in an ROH; but first need to remove ROH */
        vector < int > temp_autostatus; temp_autostatus.assign(auto_status.begin(),auto_status.end());
        if(linenum == 1) /* need to set up tempchr and temppos */
        {
            int row = chr.size();                                           /* Figure out number of rows */
            int i = 0;                                                      /* counter to determine where you are at */
            while(i < row)
            {
                while(1)
                {
                    if(temp_autostatus[i] == 5)
                    {
                        temp_chr.erase(temp_chr.begin()+i);
                        temp_pos.erase(temp_pos.begin()+i);
                        temp_autostatus.erase(temp_autostatus.begin()+i);
                        row = row - 1; break;
                    }
                    else{i++; break;}
                }
            }
        }
        if(linenum > 1) /* just need to remove 5s */
        {
            int row = temp_autostatus.size();                               /* Figure out number of rows */
            int i = 0;                                                      /* counter to determine where you are at */
            while(i < row)
            {
                while(1)
                {
                    if(temp_autostatus[i] == 5){temp_autostatus.erase(temp_autostatus.begin()+i); row = row - 1; break;}
                    else{i++; break;}
                }
            }
        }
        if(RemoveSNP == "yes" && linenum == 1)
        {
            log_file << "   - Removing SNP that weren't in ROH." << endl;
            stringstream outputstring(stringstream::out);
            outputstring << "ID";
            for(int i = 0; i < temp_chr.size(); i++)
            {
                outputstring << " " << temp_chr[i] << "_" << temp_pos[i];
            }
            outputstring << endl;
            std::ofstream output2(OutputFileAuto.c_str(), std::ios_base::app | std::ios_base::out);
            output2 << outputstring.str(); outputstring.str(""); outputstring.clear();
            log_file << "   - SNP that remain: " << temp_chr.size() << "." << endl;
            log_file << "   - Number of lines read:";
        }
        if(RemoveSNP == "no" && linenum == 1)
        {
            log_file << "   - Did Not Remove SNP that weren't in ROH." << endl;
            log_file << "   - Therefore should see 5 in autozygosity file (need to remove)." << endl;
            log_file << "   - Number of SNP in autozygosity file: " << auto_status.size() << "." << endl;
            log_file << "   - Number of lines read: " << endl;
            stringstream outputstring(stringstream::out);
            outputstring << "ID";
            for(int i = 0; i < auto_status.size(); i++)
            {
                outputstring << " " << chr[i] << "_" << position[i];
            }
            outputstring << endl;
            std::ofstream output2(OutputFileAuto.c_str(), std::ios_base::app | std::ios_base::out);
            output2 << outputstring.str(); outputstring.str(""); outputstring.clear();
        }
        /* Now Generate Proportion of genome in an ROH */
        /*Store upper and lower bounds of chromosomes*/
        if(linenum == 1)
        {
            clb.push_back(0);
            for (int c=1; c< temp_chr.size(); c++)
            {
                if (temp_chr[c] != temp_chr[(c-1)])
                {
                    clb.push_back(c);
                    cub.push_back(c-1);
                }
            }
            //add end of last chromosome
            cub.push_back(temp_chr.size()-1);
            /* get total length */
            for(int i = 0; i < clb.size(); i++){genome_length += temp_pos[cub[i]] / double(1000000);}
        }
        /* Within each chromosome Figure out where start and end */
        vector < int > ROHsw;
        for (int c = 0; c < clb.size(); c++)  //for each chromosome
        {
            if(temp_autostatus[clb[c]] == 0)        /* Start is not an ROH */
            {
                for (int s = clb[c]+1; s<(cub[c]+1); s++ ) //from start to stop of chr
                {
                    if (temp_autostatus[s] != temp_autostatus[(s-1)]){ROHsw.push_back(s);}  //switch between 0/1 or 1/0
                }
            }
            else//same as above
            {
                ROHsw.push_back(clb[c]);    //first SNP is start of ROH
                for (int s =clb[c]+1;s<(cub[c]+1);s++ )
                {
                    if (temp_autostatus[s] != temp_autostatus[(s-1)]){ROHsw.push_back(s);}  //switch between 0/1 or 1/0
                }
            }
            if (ROHsw.size() % 2 != 0){ROHsw.push_back(cub[c]);}    //ROh goes to end of chromosome
        }
        /* Figure out start and end in megabases */
        vector < int > ws; vector < int > we; vector < double > dist_MB;
        for (int w=0; w < ROHsw.size(); w++)
        {
            if (w % 2==0) //odds are start, evens are ends
            {
                ws.push_back(temp_pos[ROHsw[w]]);
            }
            else
            {
                we.push_back(temp_pos[ROHsw[w]]-1);
            }
        }
        for (int i=0; i < we.size();i++){dist_MB.push_back((we[i] - ws[i]));} //length of ROh in megabases
        double proportion = 0.0;
        for (int i=0; i < we.size();i++){proportion += dist_MB[i] / double(1000000);}
        proportion = proportion / double(genome_length);
        /* If don't want it to spit out where at then just comment this out or change the number after % to something bigger */
        /* output ROH matrix */
        stringstream outputstring(stringstream::out);
        for(int i = 0; i < roh_index.size(); i++)               /* loop across and fill rohrange array */
        {
            if(i < roh_index.size() -1){outputstring << roh_status[i] << " ";}
            if(i == roh_index.size() -1){outputstring << roh_status[i];}
        }
        outputstring << endl;
        std::ofstream output1(OutputFileROH.c_str(), std::ios_base::app | std::ios_base::out);
        output1 << outputstring.str(); outputstring.str(""); outputstring.clear();
        /* output autozygosity matrix */
        if(RemoveSNP == "no")
        {
            outputstring << ID;
            for(int i = 0; i < auto_status.size(); i++)               /* loop across and fill rohrange array */
            {
                outputstring << " " << auto_status[i];
            }
            outputstring << endl;
            std::ofstream output2(OutputFileAuto.c_str(), std::ios_base::app | std::ios_base::out);
            output2 << outputstring.str(); outputstring.str(""); outputstring.clear();
        }
        if(RemoveSNP == "yes")
        {
            outputstring << ID;
            for(int i = 0; i < temp_autostatus.size(); i++)               /* loop across and fill rohrange array */
            {
                outputstring << " " << temp_autostatus[i];
            }
            outputstring << endl;
            std::ofstream output2(OutputFileAuto.c_str(), std::ios_base::app | std::ios_base::out);
            output2 << outputstring.str(); outputstring.str(""); outputstring.clear();
        }
        std::ofstream output3(OutputFileSummary.c_str(), std::ios_base::app | std::ios_base::out);
        output3 << ID << " " << homozygosity << " " << proportion << endl;
        if(linenum % 5 == 0)
        {
            cout << "  - " << linenum << endl;
            log_file << "  - " << linenum << endl;
        }
        linenum++;
    }
    cout << "- Finished calculating ROH and Autozygosity" << endl;
    cout << "- Program ended Normally" << endl;
    time_t full_end_time = time(0);
    cout << "- Time elapsed sec: " << difftime(full_end_time,full_begin_time) << endl;
}
