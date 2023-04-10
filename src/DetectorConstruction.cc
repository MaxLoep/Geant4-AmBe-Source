/*
This file contains examples on how to define materials and how to create geometries and place them.
There are no external dependencies, you find everything in this file. -> You can copy code snippets and paste them in your simulation.
Remember to include the header-files in your simulation, e.g. if you want to place a Box, you have to put ' #include "G4Box.hh" ' in your file as well.
*/


#include "DetectorConstruction.hh"      //Header file where functions classes and variables may be defined (...)
#include "DetectorMessenger.hh"         //Header file for own macro commands
#include "G4RunManager.hh"              //Necessary. You need this.

#include "G4NistManager.hh"             //for getting material definitions from the NIST database
#include "G4Material.hh"

#include "G4Box.hh"                     //for cuboid
#include "G4Sphere.hh"                  //for sphere
#include "G4Tubs.hh"                    //for cylinder
#include "G4LogicalVolume.hh"           //Necessary. You need this.
#include "G4PVPlacement.hh"             //Necessary. You need this.
#include "G4SystemOfUnits.hh"           //for units
#include "G4UnitsTable.hh"              //for units
#include "G4PhysicalConstants.hh"       //for physial constants like pi

#include "G4VisAttributes.hh"           //for Visualization
#include "G4Color.hh"                   //for Visualization

//Primitive Scorer from example B4d
#include "G4SDManager.hh"
#include "G4SDChargedFilter.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackLength.hh"
#include "G4PSFlatSurfaceCurrent.hh"

#include "G4SDParticleWithEnergyFilter.hh"
#include "G4SDParticleFilter.hh"
#include "G4SDChargedFilter.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "SensitiveDetector.hh"                     //the SensitiveDetector
#include "CADMesh.hh"                   // for importing CAD-files (.stl, .obj, ...). Read all about it at: https://github.com/christopherpoole/CADMesh


DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
 fAbsorMaterial(nullptr), fLAbsor(nullptr), world_mat(nullptr), fDetectorMessenger(nullptr),
 fScoringVolume(0)
{
  // World Size
  world_sizeXYZ = 20.*cm;

  //set box parameters
  boxX     = 10. *cm;
  boxY     = 10. *cm;
  boxZ     = 10. *cm;

  // set dummy variables
  a = 20.*cm; // used for x- and y-width of Sensitive Detectors
  b = 10.*cm; 
  c = 1.*cm;
  d = 1.*cm;
  e = 1.*cm;

  // materials
  DefineMaterials(); // see below for this function
  // SetAbsorMaterial("G4_Co");
  //Print all defined materials to console
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  // create commands for interactive definition of the geometry
  fDetectorMessenger = new DetectorMessenger(this);
}

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//Define materials and compositions you want to use in the simulation
void DetectorConstruction::DefineMaterials()
{
  //MATERIALS:
  //
  //How to define Materials using the NIST database
  //see https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/Appendix/materialNames.html for a list of available materials
  //
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // define world material as vacuum (Galactic) and boxMaterial as Copper using the NIST database
  boxMaterial  = nist->FindOrBuildMaterial("G4_Cu");

  //Print all defined materials to console
  // G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  // 
  // Materials
  // 
  world_mat  = nist->FindOrBuildMaterial("G4_Galactic");
  Vacuum     = nist->FindOrBuildMaterial("G4_Galactic");
  Steel      = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  Be         = nist->FindOrBuildMaterial("G4_Be");  //1.848 g/cm3
  O          = nist->FindOrBuildMaterial("G4_O");   //13.67 g/cm3
  Am         = nist->FindOrBuildMaterial("G4_Am");  //0.001429 g/cm3  - 13.67g/cm3 according to geant4

  //Define AmBeO via number of atoms
  //AmO2-Be
  //AmO2: 0.37g
  //Be:   4.6g
  G4int ncomponents, natoms;
  BeO =  new G4Material("BeO", 3.05*g/cm3, ncomponents=3);
  BeO->AddMaterial(Am, natoms=1);
  BeO->AddMaterial(Be, natoms=1);
  BeO->AddMaterial(O, natoms=1);

  //Define AmBeO via percent
  AmBeO= new G4Material("AmBeO",                           //name
                                      3.05*g/cm3,          //density
                                      3);                  //number of elements

  //Add Elements to Material
  AmBeO->AddMaterial(Am, 95.*perCent);
  AmBeO->AddMaterial(Be, 1.6*perCent);
  AmBeO->AddMaterial(O, 3.4*perCent);
}


G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();


  //Start with creating a World-Volume to place things in
  // 
  // World
  // 
  // World box where the simulation takes place
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXYZ, 0.5*world_sizeXYZ, 0.5*world_sizeXYZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //boolean operation?
                      0,                     //copy number
                      true);                 //overlaps checking?

  //Make world-volume invisible
  auto logicWorldVisAtt = new G4VisAttributes(G4Color(1, 1, 1, 0.01)); //(r, g, b , transparency)
  logicWorldVisAtt->SetVisibility(true);
  logicWorld->SetVisAttributes(logicWorldVisAtt);


// 
//AmBe Neutron Source - Capsule type X21 - Emission: ~2.2x10^6 n/s/Ci (Amersham) -> ~6x10^-5 n/s/Bq (ISO 8529-1:2021) and 1Ci = 37 GBq
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#pragma region
  //
  // Cylinder - outer Steel Cylinder 
  //
  G4Tubs* sCylinder_out = 
    new G4Tubs("Cylinder_out",                   //name
              0., 7.8/2*mm,                      //inner radius, outer radius
              15./2*mm,                          //z half length
              0., twopi);                        //min phi, max phi

  G4LogicalVolume* lCylinder_out = 
    new G4LogicalVolume(sCylinder_out,           //shape
                        Steel,                   //material
                        "Cylinder_out");         //name

  new G4PVPlacement(0,                           //no rotation
              G4ThreeVector(0,0,0),              //position
              lCylinder_out,                     //logical volume
              "Cylinder_out",                    //name
              logicWorld,                        //mother  volume
              false,                             //boolean operation?
              0,                                 //copy number
              true);                             //overlaps checking?

  //Make (in-)visible and give it a color
  auto lCylinder_outVisAtt = new G4VisAttributes(G4Color(0, 0, 1, 0.9)); //(r, g, b , transparency)
  lCylinder_outVisAtt->SetVisibility(true);
  lCylinder_out->SetVisAttributes(lCylinder_outVisAtt);

  //
  // Cylinder - inner Beryllium-Oxide Cylinder 
  //
  G4Tubs* sCylinder_in = 
    new G4Tubs("Cylinder_in",                    //name
              0., 4.6/2*mm,                      //inner radius, outer radius
              11.8/2*mm,                         //z half length
              0., twopi);                        //min phi, max phi

  G4LogicalVolume* lCylinder_in = 
    new G4LogicalVolume(sCylinder_in,            //shape
                        BeO,                     //material
                        // AmBeO,                     //material
                        "Cylinder_in");          //name

  new G4PVPlacement(0,                           //no rotation
              G4ThreeVector(0,0,0),              //position
              lCylinder_in,                      //logical volume
              "Cylinder_in",                     //name
              lCylinder_out,                     //mother  volume
              false,                             //boolean operation?
              0,                                 //copy number
              true);                             //overlaps checking?

  //Make (in-)visible and give it a color
  auto lCylinder_inVisAtt = new G4VisAttributes(G4Color(1, 0, 0, 0.8)); //(r, g, b , transparency)
  lCylinder_inVisAtt->SetVisibility(true);
  lCylinder_in->SetVisAttributes(lCylinder_inVisAtt);

#pragma endregion
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------


  //
  // Sphere - SD to detect gammas
  //
  G4Sphere* sSphere =
    new G4Sphere("Sphere",                    //name
              2.*cm, 2.1*cm,                  //inner radius, outer radius
              0., twopi,                      //min phi, max phi
              0., pi);                        //min rho, max rho
            
  G4LogicalVolume* lSphere =
    new G4LogicalVolume(sSphere,              //shape
                        Vacuum,             //material
                        "Sphere");            //name

  new G4PVPlacement(0,                        //no rotation
              G4ThreeVector(0,0,0),           //position
              lSphere,                        //logical volume
              "Sphere",                       //name
              logicWorld,                     //mother volume
              false,                          //boolean operation?
              0,                              //copy number
              true);                          //overlaps checking?

  //Make (in-)visible and give it a color
  auto lSphereVisAtt = new G4VisAttributes(G4Color(0, 1, 0, 0.1)); //(r, g, b , transparency)
  lSphereVisAtt->SetVisibility(true);
  lSphere->SetVisAttributes(lSphereVisAtt);


  PrintParameters();
  
  //always return the root volume
  //
  return physWorld;
}

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n The Absorber is " << G4BestUnit(boxX,"Length")
         << " of " << AmBeO->GetName() 
         << "\n \n" << AmBeO << G4endl;
}


//
//Functions for custom GUI and macro commands - see DetectorConstruction.hh, DetectorMessenger.cc, DetectorMessenger.hh
//
void DetectorConstruction::SetAbsorMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  
  if (pttoMaterial) { 
    dummyMat = pttoMaterial;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
    G4cout << "\n The dummyMat is now "
           << dummyMat->GetName() 
           << "\n \n" << dummyMat << G4endl;

    if(fLAbsor) { fLAbsor->SetMaterial(fAbsorMaterial); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    G4cout << "\n Weird if-case happened..." << G4endl;
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}

// Change Parameters via Macro file with these
// Change a
void DetectorConstruction::change_a(G4double value)
{
  a = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
  G4cout  << "\n a is now " << G4BestUnit(a,"Length") << G4endl;
}

// Change b
void DetectorConstruction::change_b(G4double value)
{
  b = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
  G4cout  << "\n b is now " << G4BestUnit(b,"Length") << G4endl;
}

// Change c
void DetectorConstruction::change_c(G4double value)
{
  c = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
  G4cout  << "\n c is now " << G4BestUnit(c,"Length") << G4endl;
}

// Change d
void DetectorConstruction::change_d(G4double value)
{
  d = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
  G4cout  << "\n d is now " << G4BestUnit(d,"Length") << G4endl;
}

// Change e
void DetectorConstruction::change_e(G4double value)
{
  e = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
  G4cout  << "\n e is now " << G4BestUnit(e,"Length") << G4endl;
}


//
//Assign Detectors and Scorers to Volume
//
void DetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);


  auto sphereSD = new SphereSD("SphereSD");                   //create a new Sensitive Detector
  G4SDManager::GetSDMpointer()->AddNewDetector(sphereSD);     //add new SD to SDManager
  SetSensitiveDetector("Sphere", sphereSD);                   //Apply Sensitive Detector 'SphereSD' to Volume 'Box'


}