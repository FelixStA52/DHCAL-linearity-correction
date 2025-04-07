/*  2024.08.19: This program will create a hits distribution and fit a gaussian to the peak. The fit parameters are then saved in an file. 

            Instructions on how to run the program:
                1. Select either to run through the "good runs" list or a specific run and follow instructions below.
                2. Run with particle indentification or not. Explaination below
                3. Edit INPUTS section (lines 91-97) to match desired effect. An explanation of each parameter can be found below. 
                4. Set the path to the appropriate directory, line 106
                5. Run through root. 
                6. Sometimes have to quit root to re-run the program


            INPUTS Parameters;
                - 3 parameters first are for the start of shower criteria.
                    numinteractionlayersRequirement: defines the minimum number of consecutive layers in which we need to see 4+ hits to define the start of a shower
                    minshowerstartRequirement: (inclusive) minimum z layer in which the shower can start
                    maxshowerstartRequirement: (inclusive) maximum z layer in which the shower can start, expressed as a percentage of the total number of layers
                    To run with the start of shower criteria, a condition should be added on line .
                - fitThreshold: gaussian fitting only applied to peak, this sets what percentage of the to peak'S value is kept. If peak value is 1000, all bins with value above fitThreshold*1000 will be used in the fitting.


            Good run vs. Sprecific Run:
            1. Using the "good runs" list which itterates through the selected runs for a period. Selected runs are identified by their run numbers.
            2. Specifiying a run and a period for which one wants to run this code.

            To run with the "good runs" list:
                - line 102: run variabe should be an index of the selected good runs list, eg. string run = Jan2011GoodRuns[i];
                - line 91: period variable should match with the selected list.
                - line 100-102: itteration through runs should not be commented out and the max number of itterations should match the size of the selected run list.
                - line 540: end of itteration through the runs should not be commnented out. 

            To run with a specific run:
                - line 92; input the run to be analysed
                - line 91: period variable should match with the selected list.
                - line 100-100: itteration through runs should be commented out
                - line 540: end of itteration through the runs should not be commnented out. 


            Running with particle selection:
                When running with the particle selection, the program will itterate through all the 3 particle types and produce 3 seperate plots
                To activate: 
                - line 93: particleItteration = 1;
                - lines 160-174 and 373-376: make sure lines are un-commented, line 153 should be commented
                - line 387: in the criterias for filling the histograms and scatter plots, make sure there is an the [if (criteriasMet) condition] is there to check if event fits with the pariticle indentification criteria
                - line 390: make sure the amount of } fit the number of conditions for filling the histogram 
                - line 539: make sure line is un-commented
*/


/* Libraries*/
#include <cmath>
#include <string>     
#include <iostream>
#include <sstream>
#include <string>
#include <ctime>            // to get the date -- in getCurrentDate()
#include <iomanip>          // to convert the date to the right format -- in getCurrentDate()
#include <cstdio>           // used in the getUniqueCherekov() function to create the pipe, fgets, and interact with the OS
#include <vector>           
#include <sstream>          
//#include <omp.h>            // openMP used to speed up getSurfaceArea() function 
#include <unordered_map>    // used  disc surface analysis
#include <filesystem>

using namespace std;

// Helper functions
string getCurrentDate();
void setupStyle(int style);
int floorDoubleToInt(double value);
double getSurfaceArea(int array[96][96][52]);
double calculateSurfaceArea(int array[96][96][52]);
int get_energy(const string &period, int run_number);
vector<int> getUniqueCherekov(const string &filename);
void resetArray(int spacialRepresentation[96][96][52]);
void resetPositionArray(int postion[4000][3], int size);
string formatNumberWithLeadingZeros(int number, int width);
int setScale(const string &type, const string &period, int energy);
std::pair<int, int> SetFitRange(TH1F *histogram, double thresholdPrecentage);
tuple<double, double, double, int, int> getCriteria(int energy, int particle, string period);
void fitAndSaveParameters(TH1F* histrogram, double fitThreshold, string fileToSaveTo, int energy, string period, string run, string particle);
int selectionConditions(int particle, double surfaceArea, double ratio10, double hits, double electronCriteriaMin, double electronCriteriaMax, double muonSurfaceCriteria, int muonHitsCriteria, int pionHitsCriteria);
void setupHistogram(TH1* hist, const char* xTitle, const char* yTitle, double xTitleSize, double yTitleSize, double xLabelSize, double yLabelSize, double xTitleOffset, double yTitleOffset, int markerColor = -1, 
    double markerSize = -1, int markerStyle = -1, int surfacePlot = -1);


void modGaussianTrueFit(string run, int size = 3, int dimension = 2) { //string run

    vector<string> Jan2011GoodRuns = {"630175", "630129", "630138", "630148", "630170"};
    vector<string> Oct2010GoodRuns = {"600095", "600091", "600082", "600076", "600071", "600058", "600053"};

    /* INPUTS */
    string period = "Jun2011"; //period when the data was taken
    //string run = "630148"; 
    std::cout << "Run number: " << run << std::endl; // FSTA addition

    string outputDir = Form("correction_%d_%d", size, dimension);
    filesystem::create_directory(outputDir);
    
    int particleItteration = 1; // 1 for true, 0 for falsse
    int numinteractionlayersRequirement = 5;
    int minshowerstartRequirement = 1;  //inclusive
    double maxshowerstartRequirement = 0.75;  //inclusive -- refers to precentage of the detector depth
    double fitThreshold = 0.3; //when calculating the gaussian fit


    /* Run itterator if using the "good runs" list */
    //for (int runItterator = 0; runItterator < 5; runItterator++) {
    //    run = Jan2011runList[runItterator];  
    

    // Path to the data
    string filepath = Form("/Users/felix/Documents/McGill_Bsc/Winter_2025/Particle/ANL/%s/run%s.beam.txt", period.c_str(), run.c_str());
    if (!filesystem::exists(filepath)) {
        cerr << "File does not exist: " << filepath << endl;
        return;
    } cout << "File exists: " << filepath << endl;


    /* UTILISING DATE AND ENERGY FUNCTIONS */
    string date = getCurrentDate(); //date the program was run (to identify output files)
    int energy = get_energy(period, atoi(run.c_str())); //beam energy in GeV -- get_energy function need an integer for run_number
    if (energy == -1) {
        cerr << "get_energy() error: period " << period << ", run " << run << endl;
        return;
    } cout << "period: " << period << ", run: " << run << ", energy: " << energy  << endl;

    // Adding leading zeros for file names
    string energy_char = formatNumberWithLeadingZeros(energy, 3);


    /* DETECTOR INFORMATION */
    string detector_name = "DHCAL";
    // number of layers by direction
    int numX = 96;  
    int numY = 96;
    int numZ = 52;
    int detectorLength;
    if (period == "Oct2010") detectorLength = 37;
    else detectorLength = 45.;


    /* SET ENERGY SCALES (COMMON FOR ALL PLOTS) */
   int xScale = setScale("xscale", period, energy);
    int scaleRMS = setScale("rms", period, energy);
    int scaleSurface = setScale("surface", period, energy);
    int scaleFragment = setScale("fragment", period, energy); 
    int scaleRatioFragment = setScale("fragment hits", period, energy);

    cout << "xscale: " << xScale << ", rmsScale: " << scaleRMS << ", surface scale: " << scaleSurface << ", fragment scale " << scaleFragment << ", fragment/hits scale: " << scaleRatioFragment << endl;


    /* INTRODUCTION OF PARTICLE VARIALBES */
    double electronCriteriaMin;
    double electronCriteriaMax;
    double muonsSurfaceCriteriaMax;
    int muonsHitsCriteriaMin;
    int pionsHitsCriteriaMax;
    string particle = "Control";  // if running without particle selection, value will be for all data, hence the control
    //int particleType; // to be commented out if using the particle itteration


    /* ITTERATING THROUGH THE 3 PARTICLE TYPES 
    1-electtron
    2-muons
    3-pions */
    for (int particleType = 1; particleType < 4; particleType++) {
        cout << " Entring loop for particle " << particleType << endl;
        if (particleType==1) particle="Electron";
        else if (particleType==2) particle="Muons";
        else if (particleType==3) particle="Pions";

        // GET SELECTION CRITERIA 
        tuple<double, double, double, int, int> particleSelectionCriteria = getCriteria(energy, particleType, period);
        electronCriteriaMin = get<0>(particleSelectionCriteria);
        electronCriteriaMax = get<1>(particleSelectionCriteria); 
        muonsSurfaceCriteriaMax = get<2>(particleSelectionCriteria); 
        muonsHitsCriteriaMin = get<3>(particleSelectionCriteria);  
        pionsHitsCriteriaMax = get<4>(particleSelectionCriteria);   
    
        cout << "min electron criteria: " << electronCriteriaMin << ", max electron criteria: " << electronCriteriaMax << ", muon surface criteria: " << muonsSurfaceCriteriaMax << ", muons min hits criteria: " << muonsHitsCriteriaMin<< ", pions max hits criteria: " << pionsHitsCriteriaMax << endl;
    

        /* READING THE FILE */
        fstream file;
        file.open(filepath, ios::in);
        if (!file.is_open()) {
            cerr << "Error: Unable to open file\n";
            return; 
        }

        /* SETTING UP FIGURE NAME/TITLE AND HISTOGRAMS */
        string FileName;
        string FigureTitle;
        string outputFileName;
        if (particleItteration == 0) {
            FileName = Form("All_plots_%s_%sGeV_%s_run%s_%s.png", detector_name.c_str(), energy_char.c_str(), period.c_str(), run.c_str(),  date.c_str());
            FigureTitle = Form("All PLots %s %s %dGeV run%s", detector_name.c_str(), period.c_str(), energy, run.c_str());  // if using any cuts, should add them in the title
            outputFileName = Form("Fit_Params_%s_%s_%s_%s.txt", particle.c_str(), detector_name.c_str(), period.c_str(), date.c_str());
        }
        else if (particleItteration == 1) {
            FileName = Form("%s/All_plots_%s_selection_%s_%sGeV_%s_run%s_%s.png", 
                outputDir.c_str(), particle.c_str(), detector_name.c_str(), energy_char.c_str(), period.c_str(), run.c_str(), date.c_str());
            outputFileName = Form("%s/Fit_Params_%s_%s_%s_%s.txt", 
                outputDir.c_str(), particle.c_str(), detector_name.c_str(), period.c_str(), date.c_str());
            
            // Updated FigureTitle with parameters
            FigureTitle = Form("%s: %d GeV, %s, Run %s, Size: %d, Dimension: %d", 
                particle.c_str(), energy, period.c_str(), run.c_str(), size, dimension);
        }
        else {cout << "Error: particleItteration should have value 0 or 1" << endl;}

        cout << FileName << endl;
        cout << FigureTitle << endl; 


        /* SETTTING UP THE HISTOGRAM */
        TH1F *histHits = new TH1F(Form("histHits_%d", particleType), "", 100, -0.5, xScale);  // syntax: name, title, # of bins, lower bound, upper bound  
        TH1F *histCorrectedHits = new TH1F(Form("histCorrectedHits_%d", particleType), "", 100, -0.5, xScale);

        /* VARIABLES */
        int x,y,z;
        double t = 0;
        int hits = 0;
        double correctedHits = 0;
        int Event = 0;              //Index to find quantities that are inside vectors?
        int maxZ = 0;               //Greatest Z-position detected
        int TriggerTime = 0;        //Used to find duration of each event
        int spacialRepresentation[96][96][52];  //Used to calculate surface area
        double ratio5;              //ratio using the first 5 X channels
        double ratio10;             //etc.
        double ratio15;
        double HitsFirst5 = 0;      //Used to get ratios
        double HitsFirst10 = 0;
        double HitsFirst15 = 0;
        double Dispersion = 0;
        double MaxDisp=0;           //Maximal dispersion
        double maxDispZ[60];        //Keep track of maximum dispersion value in each Z-layer
        double cDisp=5;             //Use to determines the depth of the event based on the maximum dispersion values in each Z-layer
        double noHitsZ[60];         //Keep track of # hits in each Z-layer
        double RMS=0;
        double Sumx = 0; 
        double Sumy = 0; 
        double RMSx = 0;
        double RMSy = 0;
        double surfaceArea = 0;
        vector<int> X;               //Positions of each hit in event (kept in a vector so that the RMS can be calculated at the end of the event since the RMS formula requires both 
        vector<int> Y;               //     the individual position of each hit and the mean position calculated at the end, it is arduous to calculate otherwise).
        vector<int> Z;
        double meanxEvent;           //mean x-position of all hits in each event (vector size = # of events)
        double meanyEvent;
        bool isHit = false;
        resetArray(spacialRepresentation);  //Set all spacialRepresentation values to zero

        //Variables for start of shower analysis
        int countnot19or20 = 0;
        int count19 = 0;
        int count20 = 0;
        bool condition1meet = false;
        bool condition2meet = false;
        int totalevents = 0;
        int totalevents1stlayer = 0;
        std::vector<int> Zlayers;
        std::vector<int> Zvector;   //vector to store Z values (layers) for each event
        bool new_event = true;
        bool new_events = true;
        double averageOccurrence = 0;
        bool datacondition = false;
        bool eventheader = false;
        bool hitcondition = false;
        std::map<int, int> countMap;    //create a map to store the count of each number in the Zvector
        std::set<int> Zset;             //create a set to store unique numbers that appear four or more times in the Zvector
        std::vector<int> Zvec;          //create a vector to store the numbers that appear four or more times in the Zvector
        bool consecutive = false;       //initialize a boolean variable to track if there are consecutive numbers in Zvec
        int startNumber;
        int endNumber;
        int numConsecutive = 0;

        //Variables for intraction analysis
        double currentInteraction = 0;
        double totalInteraction = 0;
        const int positionArraySize = 4000;
        int position[4000][3];
        resetPositionArray(position, positionArraySize);

        //Variables for Fragment Analysis
        std::unordered_map<int, std::pair<double, int>> Zmap;

        while(file >> t >> x >> y >> z) {

            if  (x==-1 && z==-1) {  // the present line is not a hit but start of a new event
                if(hits>0){   //i.e. the event that just ended has at least one hit. 
                    //cout << "Events: " << Event << endl;
                    //cout << "hits: " << hits << endl;
                    Event++; 
                    ratio5=HitsFirst5/hits;
                    ratio10=HitsFirst10/hits;
                    ratio15=HitsFirst15/hits;
                
                    //Save mean X,Y positions in the corresponding vectors
                    meanxEvent = Sumx/hits;
                    meanyEvent = Sumy/hits;

                    correctedHits = 0.0; // Reset for the current event
                    for (int n = 0; n < hits; ++n) {
                        int x = X[n];
                        int y = Y[n];
                        int z = Zvector[n];
                        int m = size;
                        int halfM = m / 2;
                        int D_cube = 0;
                        int n_cube = 0;
                        
                        // Count hits in mxm slice around (x, y, z)
                        if (dimension == 3) {
                            for (int dx = -halfM; dx <= halfM; ++dx) {
                                for (int dy = -halfM; dy <= halfM; ++dy) {
                                    for (int dz = -halfM; dz <= halfM; ++dz) {
                                        // 3D cube processing
                                        int xi = x + dx, yi = y + dy, zi = z + dz;
                                        if (xi >= 0 && xi < 96 && yi >= 0 && yi < 96 && zi >= 0 && zi < 52) {
                                            if (spacialRepresentation[xi][yi][zi] == 1) {
                                                D_cube++;
                                            }
                                        }
                                    }
                                }
                            }
                            n_cube = m * m * m;
                        } else { // 2D
                            for (int dx = -halfM; dx <= halfM; ++dx) {
                                for (int dy = -halfM; dy <= halfM; ++dy) {
                                    int xi = x + dx, yi = y + dy, zi = z;
                                    if (xi >= 0 && xi < 96 && yi >= 0 && yi < 96 && zi >= 0 && zi < 52) {
                                        if (spacialRepresentation[xi][yi][zi] == 1) {
                                            D_cube++;
                                        }
                                    }
                                }
                            }
                            n_cube = m * m;
                        }
                        
                        if (D_cube > 0) {
                            if (D_cube >= n_cube){
                                double h_cube = log(1.0 - (double)(n_cube-0.37)/n_cube) / log(1.0 - 1.0/n_cube);
                                correctedHits += h_cube/D_cube;
                            }else{
                                double h_cube = log(1.0 - (double)D_cube/n_cube) / log(1.0 - 1.0/n_cube);
                                correctedHits += h_cube/D_cube;
                            }
                        }
                    }

                    //Calculate RMSx and RMSy for each hit of the event
                    for(int n = 0; n < X.size(); n++){
                        RMSx=RMSx + pow((X[n]-meanxEvent),2.0);
                        RMSy=RMSy + pow((Y[n]-meanyEvent),2.0);
                    }
                    RMS = sqrt((RMSx+RMSy)/(2*hits));
                    
                    if(MaxDisp<3){
                        cDisp = 1;
                    } 
                    else if(MaxDisp<10){
                        cDisp = 2; }

                    //Calculate Surface Area
                    surfaceArea = calculateSurfaceArea(spacialRepresentation);

                    //Calculating Start of Shower
                    std::sort(Zvector.begin(), Zvector.end()); // Sort the Zvector in ascending order
                    for (int Zs : Zvector) { // Iterate through the layers vector for the current event
                        countMap[Zs]++; // Increment count for the current number in the map
                    }
                    for (const auto& pair : countMap) { 
                        if (pair.second >= 4) { // If the layer appears four or more times, which means there are four or more hits appearing in that layer
                            Zset.insert(pair.first); // Add the number to the new vector
                            Zvec.push_back(pair.first);
                        }
                    }
                    int setSize = Zset.size(); // Get the size of the Zset

                    // Check if Zvec has any numbers
                    if (Zvec.size() > 0) {
                        size_t j = 0;
                        while (j < Zvec.size() - 1) {
                            startNumber = Zvec[j];
                            endNumber = startNumber;
                            // Find the end of the consecutive sequence
                            while (j < Zvec.size() - 1 && Zvec[j] + 1 == Zvec[j + 1]) {
                                ++j;
                                endNumber = Zvec[j];
                            }
                            // Check if the consecutive sequence has a length greater than 1
                            if (endNumber - startNumber > 1) {
                                //cout << "It does" << endl;
                                numConsecutive += endNumber - startNumber; // Increment the numConsecutive variable
                            }
                            if (numConsecutive > numinteractionlayersRequirement) { // checking here so that as soon as there is a consecutive series of layers that fit the criteria we stop imidiately. 
                                consecutive = true; 
                                break;
                            }
                            ++j;
                        }
                    }
                    
                    if (numConsecutive > numinteractionlayersRequirement) {
                        consecutive = true;

                        /* CENTER OF INTERACTION ANALYSIS -- very time consuming so it was commented out
                        for (int i=0; i<hits-1; i++) {
                            if (position[i][2] < startNumber) continue;

                            for (int j=i+1; j<hits; j++) {

                                if (position[j][2] < startNumber) continue;
                                currentInteraction += sqrt(abs( pow(position[i][0]-position[j][0],2) + pow(position[i][1]-position[j][1],2) + pow(position[i][2]-position[j][2],2) ));
                            }
                        }
                        totalInteraction += currentInteraction;*/
                    }

                    // Fragment Analysis (or Disc Volume)
                    double totalFragmentSum = 0;
                    for (const auto& entry : Zmap) {
                        double maxValue = entry.second.first;
                        int Zoccurance = entry.second.second;
                        totalFragmentSum += Zoccurance*maxValue;
                    }

                    /* Selection criteria per particle (if itterating through particles) */
                    int criteriaMet;
                    if (particleItteration) {
                        criteriaMet = selectionConditions(particleType, surfaceArea, ratio10, hits, electronCriteriaMin, electronCriteriaMax, muonsSurfaceCriteriaMax, muonsHitsCriteriaMin, pionsHitsCriteriaMax);
                    }


                    /* FILLING THE PLOTS 
                    Example of conditions you could have:
                        - if (consecutive == true) 
                        - if (startNumber >= minshowerstartRequirement && startNumber <= floorDoubleToInt(detectorLength*maxshowerstartRequirement)) 
                        - if (totalInteraction != 0) 
                        - if (criteriaMet==1)  -- Used when doing particle itteration
                        - if (hits >200) */

                    if (criteriaMet == 1 && hits > 30) {
                        //Save value in the histogram:
                        histHits->Fill(hits);
                        histCorrectedHits->Fill(correctedHits);
                    }

                    //Reset all appropriate variables so that the next event can be processed:
                    for(int w = 0; w <= 59; w++){
                        noHitsZ[w]=0;
                        maxDispZ[w]=0;
                    }
                    cDisp=5;
                    hits=0;
                    correctedHits = 0;
                    HitsFirst5=0;
                    HitsFirst10=0;
                    HitsFirst15=0;
                    maxZ=0;
                    Dispersion=0;
                    MaxDisp=0;
                    Sumx=0;
                    Sumy=0;
                    RMSx=0;
                    RMSy=0;
                    RMS=0;
                    X.clear();
                    Y.clear();
                    surfaceArea = 0;
                    resetArray(spacialRepresentation);
                    Zvector.clear();
                    countMap.clear();  
                    averageOccurrence=0;
                    Zvec.clear();
                    numConsecutive=0;

                    totalInteraction=0;
                    currentInteraction=0;
                    resetPositionArray(position, hits);

                    totalFragmentSum = 0; 
                    Zmap.clear();
                }
                else{ //i.e. the event that just ended had 0 hits
                    hits=0; 
                    correctedHits = 0;
                    RMS=0; 
                }
            }

            else {  //line corresponds to a hit
                hits++;
                spacialRepresentation[x][y][z] = 1;

                // int m = 3;  // Cube size (m x m x m)
                // int halfM = m/2;
                // int D_cube = 0;
                
                // // Count hits in cube
                // for(int dx = -halfM; dx <= halfM; dx++) {
                //     for(int dy = -halfM; dy <= halfM; dy++) {
                //         for(int dz = -halfM; dz <= halfM; dz++) {
                //             int xi = x + dx;
                //             int yi = y + dy;
                //             int zi = z + dz;
                //             if(xi >= 0 && xi < 96 && yi >= 0 && yi < 96 && zi >= 0 && zi < 52) {
                //                 if(spacialRepresentation[xi][yi][zi] == 1) D_cube++;
                //             }
                //         }
                //     }
                // }
                
                // // Calculate correction factor
                // int n_cube = m*m*m;
                // if(D_cube >= n_cube) D_cube = n_cube - 1;  // Avoid log(0)
                // if(D_cube > 0) {
                //     double h_cube = log(1.0 - (double)D_cube/n_cube) / log(1.0 - 1.0/n_cube);
                //     correctedHits += h_cube/D_cube;
                // }

                Zvector.push_back(z);
                /* center of interaction analysis
                if (hits<positionArraySize-1) {
                    position[hits-1][0]=x;
                    position[hits-1][1]=y;
                    position[hits-1][2]=z; 
                }*/

                // Updating Zmap with right value
                double value = std::pow(abs(x-48), 2) + std::pow(abs(y-28), 2);
                if (Zmap.find(z) == Zmap.end()) {  // update map
                    Zmap[z] = {value, 1};
                } else {
                    Zmap[z].first = std::max(Zmap[z].first, value);
                    Zmap[z].second += 1;
                }

                if(t >= 0){
                    //Save x and y temporarily, to later calculate RMS when the event ends:
                    X.push_back(x);
                    Y.push_back(y);
                
                    //Update highest Z attained in the event:
                    if(z > maxZ) maxZ=z; 
                    
                    //Keep track of how many hits are in the first N layers to later calculate ratioN
                    if(z < 16) HitsFirst15++;
                    if(z < 11) HitsFirst10++;
                    if(z < 6) HitsFirst5++;
                
                    //Update maximal dispersion attained in the event
                    //Dispersion is calculated as the distance between the position of hit and center of detector
                    Dispersion = sqrt(pow((x-48),2.0)+pow((y-48),2.0)); 
                    if(Dispersion > MaxDisp && x!=0 && x!=1 && x!=2 && x!=3){
                        MaxDisp=Dispersion;
                    }
                    noHitsZ[z]++;
                    if(Dispersion > maxDispZ[z]){
                        maxDispZ[z]=Dispersion;
                    }
                    
                    //Sums needed to find means in the end
                    Sumx = Sumx + x;
                    Sumy = Sumy + y;
                }
            }
            if(file.eof()) break;
        } 
        file.close();

        /* GAUSSIAN FITTING */
        fitAndSaveParameters(histHits, fitThreshold, outputFileName, energy, period, run, particle);
        fitAndSaveParameters(histCorrectedHits, fitThreshold, outputFileName, energy, period, run, particle+"_Corrected");

        /* PLOTTING */

        /* DEFINITION OF CANVAS SETTINGS */
        delete gROOT->FindObject(Form("0%sGeV%s_run%s_%s_Gaus_Fit", energy_char.c_str(),period.c_str(), run.c_str(), date.c_str()));
        TCanvas *cCanvas = new TCanvas(Form("0%sGeV%s_run%s_%s_Gaus_Fit", energy_char.c_str(),period.c_str(), run.c_str(), date.c_str()),"ComparisonPlots",2400,1200);
        cCanvas->Divide(2, 1);  // Split into 2 horizontal pads

        // Divide the canvas into 6x4 pads with custom margins
        //cCanvas->Divide(1, 1, 0.03, 0.02); 

        TPad *TCanvasTitle = new TPad("all","all",0,0.93,1,1);  // Place title at top 7% of canvas
        TCanvasTitle->SetFillStyle(4000);
        TCanvasTitle->Draw();
        TCanvasTitle->cd();

        TLatex *FigTitle = new TLatex();
        FigTitle->SetTextSize(0.5); 
        FigTitle->DrawLatexNDC(0.01, 0.5, FigureTitle.c_str());  // Center text vertically in title pad

        //1x1 Settings
        double ysize = 0.05;
        double yoffset1 = 0.8;
        double yoffset2 = 0.8;
        double yoffset3 = 1;
        double xsize = 0.05;
        double xoffset = 0.8;
        double labelsize1 = 0.04;
        double labelsize2 = 0.04;
        double labelsize3 = 0.04;
        double markersize1 = 0.03;
        double markersize2 = 0.05;
        int surfaceAxis = 505;

        setupHistogram(histHits, "Number of hits", "Number of events", 0.05, 0.05, 0.04, 0.04, 1.1, 1.2);
        setupHistogram(histCorrectedHits, "Corrected number of hits", "Number of events", 0.05, 0.05, 0.04, 0.04, 1.1, 1.2);

        // Setting up the style of the plot
        setupStyle(3);

        /* Drawing 1x1 plot */
        cCanvas->cd(1);  
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.14);
        gPad->SetBottomMargin(0.15);
        histHits->Draw();
        cCanvas->cd(2);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.14);
        gPad->SetBottomMargin(0.15);
        histCorrectedHits->Draw();

        //Saving the plots as a .png file
        cCanvas->Print(FileName.c_str());

    cout << "period: " << period << ", run: " << run << ", energy: " << energy  << endl;
    cout << "xscale: " << xScale << ", rmsScale: " << scaleRMS << ", surface scale: " << scaleSurface << ", fragment scale " << scaleFragment << ", fragment/hits scale: " << scaleRatioFragment << endl;

    }  // end of the particle looping
    //} // end of run looping

}

        

/* HELPER FUNCTIONS */      

int get_energy(const string &period, int run_number) {
    /*  
        Gets the energy for a particular run so that it does not have to be done manually 
        Inputs: period: as a constant string, period in which the run we are looking for it (Jan2011, Jun2011 or Oct2010 only)
                run_number: as int, the run number for which we want the energy
        Returns: as int, the energy of the run (-1 if run number was not found or wrong period was wrong)
    */
    string Jan2011 = "Jan2011", Jun2011 = "Jun2011", Nov2011 = "Nov2011", Oct2010 = "Oct2010", Apr2011 = "Apr2011";  // since can't compare directly to a string

    if (period == Jan2011) {
        if ((run_number >= 600111 && run_number <= 600135) || (run_number >= 600213 && run_number <= 600217)) return 32; 
        else if (run_number >= 600136 && run_number <= 600172) return 2;
        else if (run_number >= 600173 && run_number <= 600186) return 4;
        else if (run_number >= 600187 && run_number <= 600196) return 6;
        else if (run_number >= 600197 && run_number <= 600204) return 8;
        else if (run_number >= 600205 && run_number <= 600212) return 10;
        else if (run_number >= 600218 && run_number <= 600219) return 60;
        else return -1; // run_number not found
    }
    else if (period == Jun2011) { 
        if ((run_number >= 630103 && run_number <= 630123) || (run_number >= 630149 && run_number <= 630150) || (run_number >= 630157 && run_number <= 630160) || (run_number >= 630193 && run_number <= 630193)) return -1; // don't know as it is not in run_summary_jun_2011 
        else if (run_number >= 630124 && run_number <= 630136) return 40;
        else if (run_number >= 630137 && run_number <= 630145) return 50;
        else if ((run_number >= 630146 && run_number <= 630148) || (run_number >= 630151 && run_number <= 630156)) return 60;
        else if (run_number >= 630161 && run_number <= 630171) return 120;
        else if (run_number >= 630172 && run_number <= 630182 && run_number != 630174) return 32;
        else if (run_number >= 630183 && run_number <= 630186) return 16;
        else if (run_number >= 630187 && run_number <= 630192) return 8;
        else return -1; // run_number not found
    }
    else if (period == Nov2011) { 
        if (run_number == 650229 || run_number == 650233 || run_number == 650234 || run_number == 650236 || run_number == 650238 || run_number == 650242 || run_number == 650244 || run_number == 650247 || run_number == 650249 || run_number == 650251 || run_number == 650253 || run_number == 650260 || run_number == 650265 || run_number == 650266 || run_number == 650267 || run_number == 650270) return 32; 
        else if (run_number == 650418 || run_number == 650421 || run_number == 650424 || run_number == 650432 || run_number == 650435 || run_number == 650438 || run_number == 650441) return 10;
        else if (run_number == 650272 || run_number == 650273 || run_number == 650276 || run_number == 650404 || run_number == 650407 || run_number == 650409 || run_number == 650412 || run_number == 650415) return 8;
        else if (run_number == 650390 || run_number == 650393) return 6;
        else if (run_number == 650284 || run_number == 650285 || run_number == 650293 || run_number == 650297 || run_number == 650299 || run_number == 650303) return 4;
        else if (run_number == 650359 || run_number == 650360 || run_number == 650361 || run_number == 650370 || run_number == 650371) return 3;
        else if (run_number == 650318 || run_number == 650320 || run_number == 650321 || run_number == 650322 || run_number == 650355 || run_number == 650357 || run_number == 650374 || run_number == 650382 || run_number == 650385) return 2;
        else if (run_number == 650325 || run_number == 650328 || run_number == 650389 || run_number == 650395 || run_number == 650399) return 1;
        else if (run_number == 650291 || run_number == 650305 || run_number == 650306 || run_number == 650308 || run_number == 650312 || run_number == 650315 || run_number == 650316 || run_number == 650340 || run_number == 650343) return 16;
        else if (run_number == 650344 || run_number == 650346 || run_number == 650350 || run_number == 650351) return 12;
        else if (run_number == 650335 || run_number == 650336) return 32;
        else return -1; // run_number not found
    }
    else if (period == Oct2010) {
        if (run_number >= 600032 && run_number <= 600048 && run_number != 600041 && run_number != 600042 && run_number != 600046 && run_number != 600047) return 32;
        else if (run_number >= 600049 && run_number <= 600053 && run_number != 600051) return 25;
        else if (run_number >= 600054 && run_number <= 600062 && run_number != 600056 && run_number != 600060) return 20; 
        else if (run_number >= 600063 && run_number <= 600071 && run_number != 600066 && run_number != 600067 && run_number != 600068) return 16; 
        else if (run_number >= 600073 && run_number <= 600080 && run_number != 600074 && run_number != 600078) return 12; 
        else if (run_number >= 600082 && run_number <= 600084) return 8; 
        else if (run_number >= 600086 && run_number <= 600093 && run_number != 600090) return 4;  
        else if (run_number >= 600094 && run_number <= 600096) return 2; 
        else if (run_number >= 600097 && run_number <= 600098) return 10; 
        else return -1; // run_number not found  
    }
    else if (period == Apr2011) {
        if (run_number >= 630074 && run_number <= 630078) return 40; // don't know as it is not in run_summary_jun_2011 
        else return 50; // run_number not found  
    }
    else return -1; // wrong period
}


string getCurrentDate() {
    // Get the current time
    time_t t = time(nullptr);
    tm* now = localtime(&t);

    ostringstream oss;   //Create a stringstream to format the date
    oss << put_time(now, "%Y%m%d");

    return oss.str();
}



string formatNumberWithLeadingZeros(int number, int width) {
    /*  
        Function that formats a number with leading zeros
        Inputs: number: integer to extend with zeros
                width: how many char are wanted
        Returns: string with the leading zeros
    */
    string numberStr = to_string(number);
    int numZeros = width - numberStr.length();

    if (numZeros > 0) return string(numZeros, '0') + numberStr;
    else return numberStr;
}


double calculateSurfaceArea(int array[96][96][52]) {
    /*  
        Function that takes a 3D array of an event an input with DHCAL size and returns it's surface area.
        This is done by itterating through all cells corresponding to a hit and identifying for all 6 neighbour cells which are hits are which are empty.
        An empty neighbour cell is counted as a surface. Take into consideration the real life dimension of the cells (1cm in x and y, 1inch in z).
        If hit cell is at the edge of the detector, then that edge is considered as an empty cell and thus added in the surface area calculation.
        Input: 3D array of hits for a DHCAL event 
        Returns: surface area of the event
    */
    int size_i = 96;
    int size_j = 96;
    int size_k = 52;
    int surfaceAreaXY = 0;   // if iterating in x, surface will be in the zy plane (^x is normal direction of the plane)
    int surfaceAreaYZ = 0;
    int surfaceAreaXZ = 0;
    double totalArea = 0;
    int count = 0;

    //#pragma omp parallel for reduction(+:surfaceAreaXY, surfaceAreaYZ, surfaceAreaXZ, count) collapse(3)
    for (int i = 0; i < size_i; i++) {
        for (int j = 0; j < size_j; j++) {
           for (int k = 0; k < size_k; k++) { 
                if (array[i][j][k] == 1) {
                    count++;

                    // Check x-direction (need to be carful around edges of the detector)
                    if (i == size_i - 1 || array[i + 1][j][k] == 0) surfaceAreaYZ++;
                    if (i == 0 || array[i - 1][j][k] == 0) surfaceAreaYZ++;

                    // Check y-direction
                    if (j == size_j - 1 || array[i][j + 1][k] == 0) surfaceAreaXZ++;
                    if (j == 0 || array[i][j - 1][k] == 0) surfaceAreaXZ++;

                    // Check z-direction
                    if (k == size_k - 1 || array[i][j][k + 1] == 0) surfaceAreaXY++;
                    if (k == 0 || array[i][j][k - 1] == 0) surfaceAreaXY++;
    }   }   }   }
    totalArea = 1. * surfaceAreaXY + 2.54 * surfaceAreaYZ + 2.54 * surfaceAreaXZ;    // 2.54 is inches to cm conversion
    return totalArea;
}


void resetArray(int spacialRepresentation[96][96][52]) {
    /* 
        Function itterates through the array provided and sets all values to zero
        Input: spacialRepresentation: 3D array with DHCAL dimensions
        Returns: void, nothing is returned
    */
    for (int i = 0; i < 96; i++) {
        for (int j = 0; j < 96; j++) {
            for (int k = 0; k < 52; k++) {
                spacialRepresentation[i][j][k] = 0;
    }   }   }
}

void resetPositionArray(int postion[4000][3], int size) {
    /* 
        Function itterates through the array provided and sets all values to zero
        Input: position: 3D array with capacy of 4000 hits and x,y,z coords for each hit
        Returns: void, nothing is returned
    */
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < 3; j++) {
            postion[i][j] = 0;
    }   }
}


int setScale(const string& type, const string& period, int energy) {
    /* 
        Function used to streamline setting the various scales for the plots
        Inputs: type: determines what the value being scaled is (xscale, rms, surface)
                period: period in which run is (used in surface scale)
                energy: energy of the current run (used for xscale and rms)
        Returns: the scale for that particular type depending on the energy or period
    */
    int variableScale = -1;

    if (type=="xscale") {
        if (energy>=120) variableScale=2800; //FSTA mod
        else if (energy>=60) variableScale=2000;
        else if (energy>=50) variableScale=1500;
        else if (energy>=40) variableScale=1250;
        else if (energy>=30) variableScale=1100;
        else if (energy>=20) variableScale=800;
        else if (energy>=15) variableScale=600;
        else if (energy>=8) variableScale=500;
        else if (energy>=4) variableScale=300;
        else if (energy>=2) variableScale=250; 
        else variableScale=600;
    } 
    else if (type=="rms") {
        if (energy>=50) variableScale=30;
        else if (energy>=4) variableScale=20;
        else variableScale=16;
    } 
    else if (type=="surface") {
        if (energy>=120) variableScale=15000;
        else if (energy>=60) variableScale=13000;
        else if (energy>=50) variableScale=9000;
        else if (energy>=40) variableScale=8000;
        else if (energy>=30) variableScale=5500;
        else if (energy>=20) variableScale=4000;
        else if (energy>=15) variableScale=3000; 
        else if (energy>=8) variableScale=2000;
        else if (energy>=2) variableScale=1500;
        else variableScale=1200;
    } 
    else if (type=="fragment") {
        if (energy>=40) variableScale = 2e+6;
        else if (energy>=20) variableScale =1e+6;
        else variableScale = 3e+5; 
    } 
    else if (type=="fragment hits") {
        if (energy>=60) variableScale=6500;
        else if (energy>=50) variableScale=4000;
        else if (energy>=20) variableScale=3500;
        else if (energy>=15) variableScale=2000;
        else if (energy>=8) variableScale=1800;
        else variableScale = 1700; 
    } else {
        cerr << "Error: setScale() requires type: xscale, rms, surface, fragment or fragment hits" << endl;
    }
    return variableScale;
}


void setupStyle(int style) {
    /*
        Different ROOT Styles for the Cavas depending on how many subplots there are
        1: style for 4x4 grid (used until 20240702)
        2: style for 6X4 grid 
        3: style for 1x1 grid 
    */
    switch(style) {
        case 1:
            gStyle->SetStatFont(63);
            gStyle->SetStatFontSize(9);
            gStyle->SetStatW(0.25);
            gStyle->SetStatH(0.350);
            gStyle->SetTitleFontSize(0.1);
            break;
        case 2:
            gStyle->SetStatFont(63);
            gStyle->SetStatFontSize(8); 
            gStyle->SetStatW(0.27);   
            gStyle->SetStatH(0.350);
            gStyle->SetTitleFontSize(0.1);
            break;
        case 3:
            gStyle->SetStatFont(63);
            gStyle->SetStatFontSize(14);  // Smaller font
            gStyle->SetStatW(0.15);       // Narrower width (was 0.18)
            gStyle->SetStatH(0.15);       // Height unchanged
            gStyle->SetStatX(1.0);        // Align to right edge of pad
            gStyle->SetStatY(0.88);       // Vertical position
            break;
    }
}


void setupHistogram(TH1* hist, const char* xTitle, const char* yTitle, double xTitleSize, double yTitleSize, double xLabelSize, double yLabelSize, double xTitleOffset, 
    double yTitleOffset, int markerColor = -1, double markerSize = -1, int markerStyle = -1, int surfacePlot = -1) {
    /* TH2 object inherit from TH1 and thus can use for both the histogram and the scatter plots
    
    */
    hist->GetXaxis()->SetTitle(xTitle);
    hist->GetYaxis()->SetTitle(yTitle);
    hist->GetXaxis()->SetTitleSize(xTitleSize);
    hist->GetYaxis()->SetTitleSize(yTitleSize);
    hist->GetXaxis()->SetLabelSize(xLabelSize);
    hist->GetYaxis()->SetLabelSize(yLabelSize);
    hist->GetXaxis()->SetTitleOffset(xTitleOffset);
    hist->GetYaxis()->SetTitleOffset(yTitleOffset);
    //hist->SetTitleSize(0.05, "t"); 

    if (markerColor != -1) hist->SetMarkerColor(markerColor);
    if (markerSize != -1) hist->SetMarkerSize(markerSize);
    if (markerStyle != -1) hist->SetMarkerStyle(markerStyle);
    if (surfacePlot != -1) hist->GetXaxis()->SetNdivisions(surfacePlot);

    gStyle->SetOptStat("euo"); 
    gStyle->SetOptFit(1111);
    gStyle->SetStatFormat("6.3g");   //should mean that there is no scientific notation
}


int floorDoubleToInt(double value) {
    return static_cast<int>(std::floor(value));
}


int selectionConditions(int particle, double surfaceArea, double ratio10, double hits, double electronCriteriaMin, double electronCriteriaMax, double muonSurfaceCriteria, int muonHitsCriteria, int pionHitsCriteria) {
    /*
        Function that checks if the current event fits the selection criteria for the particle it's tested agaisnt
        Input: particle:  1 for electron, 2 for muon, 3 for pion
        Returns: 1 if conditions are met, -1 if conditions are not met

    */
    int conditionMet = -1;
    if (particle==1) {
        if (electronCriteriaMin < surfaceArea && surfaceArea < electronCriteriaMax && ratio10 > 0.65) {
            conditionMet=1;
        }
    }
    else if (particle==2) {
        if (surfaceArea < muonSurfaceCriteria && ratio10 < 0.5 && hits > muonHitsCriteria) {
            conditionMet=1;
        }
    }
    else if (particle==3) {
        if (hits<pionHitsCriteria && ((ratio10 < 0.65 && surfaceArea > muonSurfaceCriteria) || (ratio10 > 0.65 && surfaceArea > electronCriteriaMax))) {
            conditionMet=1;
        }
    }
    else cout << "Error selectionConditions(): No particle specified" << endl;
    return conditionMet; 
}


tuple<double, double, double, int, int> getCriteria(int energy, int particle, string period) {
    /* 
        Function that gets the selection criteria for a particle type for a specific energy
        Inputs: energy: energy of the run from which the event is from
                particle: type that we are trying to identify (1: electron, 2: muon, 3: pion)
                period: period from which the run is from
        Returns: tuple: min surface area for electron selection
                        max surface area for electron selection
                        max surface area for muons selection
                        min hits for muons selection
                        max hits for pions selection
    */
    tuple<double, double, double, int, int> returnTuple;

    if (particle == 1) { // Electron criteria
        if (energy >= 120)      returnTuple = std::make_tuple(1500, 6000, 1500, 30, 2200);
        else if (energy >= 60)  returnTuple = std::make_tuple(1500, 4000, 1250, 30, 1500);
        else if (energy >= 50)  returnTuple = std::make_tuple(1500, 3300, 1250, 30, 1100);
        else if (energy >= 40)  returnTuple = std::make_tuple(1250, 2200, 1250, 30, 900);
        else if (energy >= 32 && period == "Jun2011")  returnTuple = std::make_tuple(1250, 2180, 1250, 30, 900);
        else if (energy >= 32 && period == "Oct2010")  returnTuple = std::make_tuple(1000, 2100, 1000, 30, 700);
        else if (energy >= 25)  returnTuple = std::make_tuple(1000, 2000, 1400, 30, 650);
        else if (energy >= 20)  returnTuple = std::make_tuple(900, 1800, 1300, 30, 500);
        else if (energy >= 16)  returnTuple = std::make_tuple(700, 1500, 900, 30, 350);
        else if (energy >= 12)  returnTuple = std::make_tuple(650, 1200, 800, 30, 300);
        else if (energy >= 8)   returnTuple = std::make_tuple(400, 1100, 750, 30, 250);
        else if (energy >= 4)   returnTuple = std::make_tuple(250, 750, 700, 30, 150);
        else if (energy >= 2)   returnTuple = std::make_tuple(100, 500, 700, 30, 160);
        else                    returnTuple = std::make_tuple(100, 500, 700, 30, 160);
    }
    else if (particle == 2) { // Muon criteria
        if (energy >= 120)      returnTuple = std::make_tuple(500, 2000, 1000, 30, 1000);
        else if (energy >= 60)  returnTuple = std::make_tuple(500, 1800, 900, 30, 800);
        else if (energy >= 50)  returnTuple = std::make_tuple(500, 1500, 800, 30, 700);
        else if (energy >= 40)  returnTuple = std::make_tuple(400, 1200, 700, 30, 600);
        else if (energy >= 32)  returnTuple = std::make_tuple(300, 1000, 600, 30, 500);
        else if (energy >= 25)  returnTuple = std::make_tuple(200, 800, 500, 30, 400);
        else if (energy >= 20)  returnTuple = std::make_tuple(150, 700, 400, 30, 300);
        else if (energy >= 16)  returnTuple = std::make_tuple(100, 600, 300, 30, 200);
        else if (energy >= 12)  returnTuple = std::make_tuple(80, 500, 250, 30, 150);
        else if (energy >= 8)   returnTuple = std::make_tuple(50, 400, 200, 30, 100);
        else if (energy >= 4)   returnTuple = std::make_tuple(30, 300, 150, 30, 80);
        else if (energy >= 2)   returnTuple = std::make_tuple(20, 200, 100, 30, 50);
        else                    returnTuple = std::make_tuple(20, 200, 100, 30, 50);
    }
    else if (particle == 3) { // Pion criteria
        if (energy >= 120)      returnTuple = std::make_tuple(2000, 8000, 2000, 30, 3000);
        else if (energy >= 60)  returnTuple = std::make_tuple(1800, 7000, 1800, 30, 2500);
        else if (energy >= 50)  returnTuple = std::make_tuple(1500, 6000, 1500, 30, 2000);
        else if (energy >= 40)  returnTuple = std::make_tuple(1200, 5000, 1200, 30, 1500);
        else if (energy >= 32)  returnTuple = std::make_tuple(1000, 4000, 1000, 30, 1200);
        else if (energy >= 25)  returnTuple = std::make_tuple(800, 3000, 800, 30, 1000);
        else if (energy >= 20)  returnTuple = std::make_tuple(600, 2500, 600, 30, 800);
        else if (energy >= 16)  returnTuple = std::make_tuple(500, 2000, 500, 30, 600);
        else if (energy >= 12)  returnTuple = std::make_tuple(400, 1500, 400, 30, 500);
        else if (energy >= 8)   returnTuple = std::make_tuple(300, 1000, 300, 30, 400);
        else if (energy >= 4)   returnTuple = std::make_tuple(200, 800, 200, 30, 300);
        else if (energy >= 2)   returnTuple = std::make_tuple(100, 500, 100, 30, 200);
        else                    returnTuple = std::make_tuple(100, 500, 100, 30, 200);
    }
    else {
        cerr << "Error: Invalid particle type in getCriteria()" << endl;
        returnTuple = std::make_tuple(0, 0, 0, 0, 0);
    }

    return returnTuple;
}

std::pair<int, int> SetFitRange(TH1F *histogram, double thresholdPrecentage) {
        // Function to determine fit range based on max bin content
        int maxBin = histogram->GetMaximumBin();
        double maxContent = histogram->GetBinContent(maxBin);
        double threshold = thresholdPrecentage * maxContent;

        cout << "maxBin: " << maxBin << ", maxContent: " << maxContent << ", threshold" << threshold << endl;

        // Determine fitMin and fitMax bins
        int fitMinBin = maxBin, fitMaxBin = maxBin;
        while (fitMinBin > 0 && histogram->GetBinContent(fitMinBin) > threshold) fitMinBin--;
        while (fitMaxBin < histogram->GetNbinsX() && histogram->GetBinContent(fitMaxBin) > threshold) fitMaxBin++;

       cout << "fitMinBin: " << fitMinBin << ", fitMaxBin: " << fitMaxBin << endl; 

        return std::make_pair(fitMinBin, fitMaxBin);
}


void fitAndSaveParameters(TH1F* histrogram, double fitThreshold, string fileToSaveTo, int energy, string period, string run, string particle) {
    /*
        Function that adds a gaussian fit to a histogram for as long as the bins are above a certain threshold of the maximal value.
        Inputs: histrogram: histogram to add fit to
                fitThreshold: percentage of the maximal value needed to be taken into consideration in the histogram
                fileToSaveTo: name of the files the fit parameters should be saved to
        Returns: Fit parameters are saved in a text file. Fit parameter should be added to the histogram when plotting 
    
    */
    // Check if third highest bin has <15 events
    vector<double> binContents;
    for (int i = 1; i <= histrogram->GetNbinsX(); ++i) {
        binContents.push_back(histrogram->GetBinContent(i));
    }
    sort(binContents.rbegin(), binContents.rend());
    if (binContents.size() >= 3 && binContents[2] < 15 && binContents[0] < 25) {
        cout << "Skipping fit and save for " << particle << endl;
        return;  // Exit without saving
    }

    // Get the fit range in terms of bin numbers
    std::pair<int, int> fitRange = SetFitRange(histrogram, fitThreshold);
    int fitMinBin = fitRange.first;
    int fitMaxBin = fitRange.second;

    // Convert bin numbers to x-values
    double fitMin = histrogram->GetBinLowEdge(fitMinBin + 1);
    double fitMax = histrogram->GetBinLowEdge(fitMaxBin);

    // Perform the fit
    TF1* fitFunc = new TF1("fitFunc", "gaus", fitMin, fitMax);
    histrogram->Fit(fitFunc, "", "", fitMin, fitMax);

    // Retrieve the fit parameters
    double param0 = fitFunc->GetParameter(0); // constant term
    double param1 = fitFunc->GetParameter(1); // mean
    double param2 = fitFunc->GetParameter(2); // sigma

    double param0Error = fitFunc->GetParError(0);
    double param1Error = fitFunc->GetParError(1);
    double param2Error = fitFunc->GetParError(2);

    // Open a text file for writing
    std::ofstream outFile(fileToSaveTo, std::ios_base::app);  // writing in append mode

    if (param1 < 0) {
        cout << "Skipping save for " << particle << endl;
        return;  // Exit without saving
    }

    // Check if the file is open
    if (outFile.is_open()) {
        // Write the parameters to the file
        outFile << particle << ", " << period << ", " << energy << ", " << run << ", " << param0 << ", " << param0Error << ", " << param1 << ", " << param1Error << ", " << param2 << ", " << param2Error << std::endl;  //Hits, Energy, Parameter 0 (constant), Parameter 1 (mean), Parameter 2 (sigma)

        // Close the file
        outFile.close();
    } else {
        std::cerr << "Unable to open file for writing" << std::endl;
    }
}