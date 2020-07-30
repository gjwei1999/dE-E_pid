//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// 
/// \file B4DetectorConstruction.cc
/// \brief Implementation of the B4DetectorConstruction class

#include "B4DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* B4DetectorConstruction::fMagFieldMessenger = nullptr; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::B4DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fSiliconPV(nullptr),
   fScintillatorPV(nullptr),
   fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::~B4DetectorConstruction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4DetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::DefineMaterials()
{ 
  // material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Air");
  nistManager->FindOrBuildMaterial("G4_Au");
  nistManager->FindOrBuildMaterial("G4_Si");
  nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4DetectorConstruction::DefineVolumes()
{

  //world size
  auto worldSizeXY = 60.0*cm;
  auto worldSizeZ  = 80.0*cm; 

  //taget size
  G4double tar_rmin =  0.*cm; 
  G4double tar_rmax = 1.*cm;
  G4double tar_h = 1.0*mm;
  G4double tar_sphi = 0.0*deg;
  G4double tar_deltaphi = 360.0*deg;
  
  //silicon size
  G4double sil_x = 10.*cm;
  G4double sil_y = 10.*cm;
  G4double sil_z = 1.534*mm;
  
  //scintillator size 
  G4double scintillator_x = 10.0*cm;
  G4double scintillator_y = 10.0*cm;
  G4double scintillator_z = 10.0*cm;
  
  //position and rotation
  G4double phi = 10.0*deg;
  G4double theta = 50.0*deg;
  G4double l = 16.185*cm;
  
  
  G4ThreeVector u2 = G4ThreeVector(std::cos(phi), std::sin(theta)*std::sin(phi), -std::cos(theta)*std::sin(phi));
  G4ThreeVector v2 = G4ThreeVector(0.0, std::cos(theta), std::sin(theta));
  G4ThreeVector w2 = G4ThreeVector( std::sin(phi), -std::sin(theta)*std::cos(phi),std::cos(theta)*std::cos(phi));
  
  G4ThreeVector u1 = G4ThreeVector(std::cos(phi), std::sin(theta)*std::sin(phi), -std::cos(theta)*std::sin(phi));
  G4ThreeVector v1 = G4ThreeVector(std::sin(phi), -std::sin(theta)*std::cos(phi), std::cos(theta)*std::cos(phi));
  G4ThreeVector w1 = G4ThreeVector(0.0, -std::cos(theta),-std::sin(theta));
  
  G4RotationMatrix rotm1  = G4RotationMatrix(u1, v1, w1);
  G4RotationMatrix rotm2  = G4RotationMatrix(u2, v2, w2);
  //for silicon
  G4ThreeVector pos1 = G4ThreeVector(l*std::cos(theta)*std::sin(phi), l*std::sin(theta), l*std::cos(theta)*std::cos(phi));
  
  G4Transform3D transform1 = G4Transform3D(rotm1,pos1);
  
  //for scintillator
  G4double l_scin = (l + 0.5*sil_z + 0.5*scintillator_z);
  
  G4ThreeVector pos2 = G4ThreeVector(l_scin*std::cos(theta)*std::sin(phi), l_scin*std::sin(theta), l_scin*std::cos(theta)*std::cos(phi));
  
  G4Transform3D transform2 = G4Transform3D(rotm1,pos2);
  
  
  // Get materials
  auto nist = G4NistManager::Instance();
  
  auto defaultMaterial = nist->FindOrBuildMaterial("G4_AIR");//G4Material::GetMaterial("G4_Air");
  auto targetMaterial = nist->FindOrBuildMaterial("G4_Au");//G4Material::GetMaterial("G4_Au");
  auto siliconMaterial = nist->FindOrBuildMaterial("G4_Si");//G4Material::GetMaterial("G4_Si");
  auto scintillatorMaterial = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");//G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  
  if ( ! defaultMaterial || ! siliconMaterial || ! scintillatorMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("B4DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }  
   
  //     
  // World
  //
  auto worldS 
    = new G4Box("World",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size
                         
  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name
                                   
  auto worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
    
    
  //
  //target 
  //
    auto targetS
        = new G4Tubs("Target",                    //its name
               tar_rmin, tar_rmax, 0.5*tar_h, tar_sphi, tar_deltaphi); //its size
    
    auto targetLV
        = new G4LogicalVolume(targetS,            //its solid
                        targetMaterial,             //its material
                        "Target");         //its name
    
    new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    targetLV,                //its logical volume
                    "Target",              //its name
                    worldLV,//logicEnv,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    fCheckOverlaps);          //overlaps checking
    
    //silicon
    auto siliconS
        = new G4Box("Silicon",                    //its name
               0.5*sil_x, 0.5*sil_y, 0.5*sil_z); //its size
        
    auto siliconLV
        = new G4LogicalVolume(siliconS,         //its solid
                        siliconMaterial,          //its material
                        "Silicon");           //its name
    
    fSiliconPV
        = new G4PVPlacement(transform1,                       //
                    siliconLV,             //its logical volume
                    "Silicon",                //its name
                    worldLV,//logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    fCheckOverlaps);          //overlaps checking
 
    //scintillator
    auto scintillatorS
        = new G4Box("Sintillator",                      //its name
              0.5*scintillator_x, 0.5*scintillator_y, 0.5*scintillator_z); //its size
        
    auto scintillatorLV
        = new G4LogicalVolume(scintillatorS,         //its solid
                        scintillatorMaterial,          //its material
                        "Scintillator");           //its name
    
    fScintillatorPV
        = new G4PVPlacement(transform2,                       //rotation
                    scintillatorLV,             //its logical volume
                    "Scintillator",                //its name
                    worldLV,//logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                          //copynumber
                   fCheckOverlaps);          //
        
    
    
  

  
  //
  // print parameters
  //
//  G4cout
//    << G4endl 
//    << "------------------------------------------------------------" << G4endl
//    << "---> The calorimeter is " << nofLayers << " layers of: [ "
// << absoThickness/mm << "mm of " << absorberMaterial->GetName() 
//    << " + "
//    << gapThickness/mm << "mm of " << gapMaterial->GetName() << " ] " << G4endl
//    << "------------------------------------------------------------" << G4endl;
  
  //                                        
  // Visualization attributes
  //
//  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

//  auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
//  simpleBoxVisAtt->SetVisibility(true);
//  calorLV->SetVisAttributes(simpleBoxVisAtt);

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//void B4DetectorConstruction::ConstructSDandField()
//{ 
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
//  G4ThreeVector fieldValue;
//  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
//  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
//  G4AutoDelete::Register(fMagFieldMessenger);
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
