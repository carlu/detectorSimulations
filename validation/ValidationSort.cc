// ValidationSort.cc - 5 Sept 2013 - C Unsworth
// To Compile:
// g++ ValidationSort.cc --std=c++0x -o ValidationSort -O2 -Wall `root-config --cflags --libs`
// To Run:
// ./ValidationSort [filename]
// Quick sort to process TIGRESS/GRIFFIN GEANT4 simulation output and produce spectra
// Which can be compared with experiment.

// C/C++ libraries:
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
using namespace std;

// ROOT libraries:
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TRandom.h>
#include <TRandom3.h>

// Custom header files for this sort
# include "ValidationSort.hh"

// Definitions
#define CRYSTALS 64   // Total number of HPGe crystals
#define CLOVERS  16   // Total number of clovers in array
#define COLOURS   4   // Number of crystals per clover (Blue,Green,Red,White)
#define MIN_EN   0.000001  // Minimum energy to be incremented to spectra

// File pointers:
static TFile *outfile = 0;

// Spectra pointers here:
static TH1F *CrystalEnergy[CRYSTALS];
static TH1F *CrystalEnergyReal[CRYSTALS];
static TH1F *ShieldEnergy[CRYSTALS]; // For now only a sum of shield energy, no individual shield elements.


int main(int argc, char **argv){
   
   // General stuff...
   int i, crystal;
   char name[256], title[256];  // For spectra 
   string str;
   // For holding the data items from each interaction....
   int DetN, SegN; 
   float Etrue, Tsec, Xmm, Ymm, Zmm; 
   char Process[24], Collection[64];
   // Derived quantities:
   float Ereal;
   float Sigma;
   // For holding data items for the whole event....   
   float CrysEn[CRYSTALS];
   float ShieldEn[CRYSTALS];
   float CrysEnReal[CRYSTALS];
   
   // Random Generator for ERes modelling
   TRandom *r1 = new TRandom3();  // TRandom3 is the closest to true random.
                                    // Faster options available if neccarsary.
   
   //  Open the input file
   ifstream infile;
   if(argc>1) {
      cout << "Opening input file: " << argv[1] << "...." << endl;
      infile.open(argv[1]);
      if(infile) {
         cout << "Success!" << endl;
      }
      else {
         cout << "Could not open file!" << endl;
         return -1;
      }   
   }
       
   // Open output file                                       
   outfile = new TFile("ValidationOut.root","RECREATE"); 
   
   // Initialise spectra...
   for(i=0;i<=CRYSTALS;i++) {
      sprintf(name,"Crys%d",i); 
      sprintf(title,"Crystal %d Energy (keV)",i); 
      CrystalEnergy[i] = new TH1F(name,title,8192,0,2048);
      sprintf(name,"Crys%d Real",i); 
      sprintf(title,"Realistic Crystal %d Energy (keV)",i); 
      CrystalEnergyReal[i] = new TH1F(name,title,8192,0,2048);
      sprintf(name,"Shield%d",i); 
      sprintf(title,"Shield %d Energy (keV)",i);      
      ShieldEnergy[i] = new TH1F(name,title,8192,0,2048);
   }
   
   // Start main loop....
   int line = 0;  // count lines to handle header
   while(!infile.eof()) {  
      getline(infile,str);
      line += 1;
      if(line<10) {continue;} // skip output header, should also skip first occurance of event header "Hits: x:"
      
      if(str.substr(0,4)=="Hits") {  // If this is a new event then fill spectra and reset event data
         for(i=0;i<CRYSTALS;i++) {
            if(CrysEn[i] > 0) {
               // Find sigma and apply Etrue -> Ereal
               Sigma = EResSigma[0] + (CrysEn[i] * (EResSigma[1] / (1332.5 - 59.5)  ) );
               CrysEnReal[i] = CrysEn[i]  + r1->Gaus(0,Sigma) ;
            } 
            if(CrysEn[i] > MIN_EN) {CrystalEnergy[i]->Fill(CrysEn[i]);}
            if(ShieldEn[i] > MIN_EN) {ShieldEnergy[i]->Fill(ShieldEn[i]);}  
            if(CrysEnReal[i] > MIN_EN) {CrystalEnergyReal[i]->Fill(CrysEnReal[i]);}
         }   
         memset(CrysEn,0.0,(CRYSTALS*sizeof(float)));
         memset(CrysEnReal,0.0,(CRYSTALS*sizeof(float)));
         memset(ShieldEn,0.0,(CRYSTALS*sizeof(float)));
         continue;
      }
      else {  // This is an interaction to be added to event record...
         // Scan input string for expected data items:
         sscanf(str.c_str(),"%i %i %f %f %f %f %f %s %s", &DetN, &SegN, &Etrue, &Tsec, &Xmm, &Ymm, &Zmm, Process, Collection); 
         
         crystal = ((DetN -1) * 4) + SegN -1; // Get crystal number from detector and segment
                     // note: "segment" in GEANT output means crystal in real world 
         //printf("%s\n%s\n",Process, Collection);
         if(strcmp(Collection,"CollectionGriffinForwardGe") || strcmp(Collection,"CollectionGriffinBackGe")) {
            CrysEn[crystal] += Etrue; // sum energy in this crystal  
         }
         else{ 
            ShieldEn[crystal] += Etrue;
         }
         
         
                     // note: At the moment this is summing HPGe and BGO energy for each colour, need to parse
                     // "Collection" string to fix this        
      } // close secttion for processing interactions   
   }  // Close main while loop
   
   // Now write the spectra....  (I think I can do these all at once but not sure how just yet.
   for(i=0;i<CRYSTALS;i++) {
      CrystalEnergy[i]->Write();
      CrystalEnergyReal[i]->Write();
      ShieldEnergy[i]->Write();
   }
   
   outfile->Close();

   return 0;
}

