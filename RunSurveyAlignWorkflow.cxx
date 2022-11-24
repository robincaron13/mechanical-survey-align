// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// Run the workflow to extract align parameters from survey.
///
//-----------------------------------------------------------------------------

#include <TFile.h>
#include <iostream>

#include "ModuleTransform.h"

#include "CCDB/CcdbApi.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "MFTBase/GeometryTGeo.h"
#include "MFTSimulation/GeometryMisAligner.h"
#include "TGeoManager.h"

namespace o2 {

namespace mft {

class SurveyAlignPar {
public:
  GeometryTGeo *mGeometryTGeo; //! access to geometry details

  //______________________________________________________________________
  void GetAlignParamsFromSurveyPositions(Bool_t verbose, Int_t level,
                                         const std::string &ccdbHost, long tmin,
                                         long tmax,
                                         const std::string &objectPath) {
    /// loop on detection elements

    mGeometryTGeo = GeometryTGeo::Instance();

    LOG(info) << " >>>> GetAlignParamsFromSurveyPositions ";

    double lPsi, lTheta, lPhi = 0.;
    Int_t nAlignID = 0;
    Int_t nChip = 0;
    Int_t nHalf = mGeometryTGeo->getNumberOfHalfs();

    // The function to get misalignment transformation matrices of DE from
    // surveyed positions
    o2::mft::ModuleTransform mTransform;

    bool emptyAlgPar = false;
    Double_t parErr[6] = {0.001, 0.001, 0.001, 0.001, 0.001, 0.001};

    TString sname = mGeometryTGeo->composeSymNameMFT();
    TGeoCombiTrans localDeltaTransform;

    mTransform.SetSensorId(0, 0, 0, 0);
    mTransform.SetSymName(sname);
    emptyAlgPar = true;
    Int_t resultTransform = mTransform.GetModuleMeanTransform(
        localDeltaTransform, parErr, 0, -1, emptyAlgPar);

    for (Int_t hf = 0; hf < nHalf; hf++) {
      Int_t nDisks = mGeometryTGeo->getNumberOfDisksPerHalf(hf);
      // New module transformation
      // localDeltaTransform = MisAlignHalf();

      sname = mGeometryTGeo->composeSymNameHalf(hf);

      mTransform.SetSensorId(hf, 0, 0, 0);
      mTransform.SetSymName(sname);
      emptyAlgPar = true;
      Int_t resultTransform = mTransform.GetModuleMeanTransform(
          localDeltaTransform, parErr, 0, -1, emptyAlgPar);

      for (Int_t dk = 0; dk < nDisks; dk++) {
        // localDeltaTransform = MisAlignDisk();
        sname = mGeometryTGeo->composeSymNameDisk(hf, dk);

        Int_t nLadders = 0;
        for (Int_t sensor = mGeometryTGeo->getMinSensorsPerLadder();
             sensor < mGeometryTGeo->getMaxSensorsPerLadder() + 1; sensor++) {
          nLadders += mGeometryTGeo->getNumberOfLaddersPerDisk(hf, dk, sensor);
        }
        if (level == 2) {

          mTransform.SetSensorId(hf, dk, 0, 0);
          mTransform.SetSymName(sname);
          emptyAlgPar = false;
          Int_t resultTransform = mTransform.GetModuleMeanTransform(
              localDeltaTransform, parErr, level, -1, emptyAlgPar);
        }
        if (level == 4) {

          mTransform.SetSensorId(hf, dk, 0, 0);
          mTransform.SetSymName(sname);
          emptyAlgPar = true;
          // Int_t resultTransform = mTransform.GetModuleMeanTransform(
          //     localDeltaTransform, parErr, level, -1, emptyAlgPar);
        }
        if (level == 10) {

          mTransform.SetSensorId(hf, dk, 0, 0);
          mTransform.SetSymName(sname);
          emptyAlgPar = false;
          Int_t resultTransform = mTransform.GetModuleMeanTransform(
              localDeltaTransform, parErr, 2, -1, emptyAlgPar);
        }

        for (Int_t lr = 0; lr < nLadders; lr++) { // nLadders
          // localDeltaTransform = MisAlignLadder();
          sname = mGeometryTGeo->composeSymNameLadder(hf, dk, lr);
          Int_t nSensorsPerLadder =
              mGeometryTGeo->getNumberOfSensorsPerLadder(hf, dk, lr);
          TString path = "/cave_1/barrel_1/" + sname;

          mTransform.SetSensorId(hf, dk, lr, 0);
          mTransform.SetSymName(sname);
          emptyAlgPar = true;
          Int_t resultTransform = mTransform.GetModuleMeanTransform(
              localDeltaTransform, parErr, 0, -1, emptyAlgPar);

          for (Int_t sr = 0; sr < nSensorsPerLadder; sr++) {
            // localDeltaTransform = MisAlignSensor();
            sname = mGeometryTGeo->composeSymNameChip(hf, dk, lr, sr);
            if (!matrixToAngles(localDeltaTransform.GetRotationMatrix(), lPsi,
                                lTheta, lPhi)) {
              LOG(error) << "Problem extracting angles from sensor";
            }
            if (level == 4) {

              mTransform.SetSensorId(hf, dk, lr, sr);
              mTransform.SetSymName(sname);
              emptyAlgPar = false;
              Int_t resultTransform = mTransform.GetModuleMeanTransform(
                  localDeltaTransform, parErr, level, nChip, emptyAlgPar);
            }
            if (level == 2) {

              mTransform.SetSensorId(hf, dk, lr, sr);
              mTransform.SetSymName(sname);
              emptyAlgPar = true;
              // Int_t resultTransform = mTransform.GetModuleMeanTransform(
              //     localDeltaTransform, parErr, level, nChip, emptyAlgPar);
            }
            if (level == 10) {

              mTransform.SetSensorId(hf, dk, lr, sr);
              mTransform.SetSymName(sname);
              emptyAlgPar = false;
              Int_t resultTransform = mTransform.GetModuleMeanTransform(
                  localDeltaTransform, parErr, 4, nChip, emptyAlgPar);
            }
            nChip++;

          } // semsors

        } // ladders

      } // disk

    } // half

    if (verbose) {
      mTransform.PrintAlignPars();
    }

    if (!ccdbHost.empty()) {
      std::string path = objectPath.empty()
                             ? o2::base::DetectorNameConf::getAlignmentPath(
                                   o2::detectors::DetID::MFT)
                             : objectPath;
      path = "MFT/test_CCDB/MFT"; // testing the ccdb
      LOG(info) << " --> Storing alignment object on ccdb, path:" << path
                << " host:" << ccdbHost;

      o2::ccdb::CcdbApi api;
      std::map<std::string, std::string> metadata; // can be empty
      metadata["test"] = fmt::format("Test survey align object for MFT:{}",
                                     o2::detectors::DetID::MFT);

      api.init(
          ccdbHost
              .c_str()); // or http://localhost:8080 for a local installation
      // store abitrary user object in strongly typed manner
      api.storeAsTFileAny(&mTransform.algVec, path, metadata, tmin, tmax);
    }
    if (ccdbHost.empty()) {
      std::string path = ""; // testing the ccdb
      LOG(info) << " --> Storing alignment object locally";
      /*
      // current date and time on the current system
         time_t now = time(0);
         // convert now to string form
         char* date_time = ctime(&now);
      LOG(info) << "The current date and time is: " << date_time ;
      */

      TFile *fout;
      if (level == 2)
        fout = new TFile("MFTSurveyAlignParDisk.root", "RECREATE");
      if (level == 4)
        fout = new TFile("MFTSurveyAlignParChip.root", "RECREATE");
      if (level == 10)
        fout = new TFile("MFTSurveyAlignParDiskChip.root", "RECREATE");

      fout->WriteObjectAny(&mTransform.algVec,
                           "std::vector<o2::detectors::AlignParam>",
                           "ccdb_object");
      // fout->WriteObject(this,"SurveyAlignPar");
      fout->Close();
    }
  }

  //______________________________________________________________________
  bool matrixToAngles(const double *rot, double &psi, double &theta,
                      double &phi) {
    /// Calculates the Euler angles in "x y z" notation
    /// using the rotation matrix
    /// Returns false in case the rotation angles can not be
    /// extracted from the matrix

    if (std::abs(rot[0]) < 1e-7 || std::abs(rot[8]) < 1e-7) {
      LOG(error) << "Failed to extract roll-pitch-yall angles!";
      return false;
    }
    psi = std::atan2(-rot[5], rot[8]);
    theta = std::asin(rot[2]);
    phi = std::atan2(-rot[1], rot[0]);
    return true;
  }

  //______________________________________________________________________
  void MisAlignFromFile(Bool_t verbose, const std::string &filenameDisk,
                        const std::string &filenameChip) {

    if (!filenameDisk.empty()) {
      TFile algDiskFile(filenameDisk.c_str(), "read");
      std::vector<o2::detectors::AlignParam> *lAPdisk;
      algDiskFile.GetObject("ccdb_object", lAPdisk);
      // loop on disks
      for (int indice = 0; indice < lAPdisk->size(); indice++) {
        o2::detectors::AlignParam &algParDisk = lAPdisk->at(indice);
        algParDisk.print();
        algParDisk.applyToGeometry();
      }
    }

    if (!filenameChip.empty()) {
      TFile algChipFile(filenameChip.c_str(), "read");

      std::vector<o2::detectors::AlignParam> *lAPchip;
      algChipFile.GetObject("ccdb_object", lAPchip);

      // loop on for sensors
      for (int indice = 0; indice < lAPchip->size(); indice++) {
        o2::detectors::AlignParam &algPar = lAPchip->at(indice);
        algPar.applyToGeometry();
      }
    }

    /*
  if (!filenameChip2.empty()) {
    TFile algChipFile2(filenameChip2.c_str(), "read");
    std::vector<o2::detectors::AlignParam>* lAPchip2;
    algChipFile2.GetObject("ccdb_object", lAPchip2);
    for (int indice = 0; indice < lAPchip2->size(); indice++) {
      o2::detectors::AlignParam& algPar = lAPchip2->at(indice);
      algPar.applyToGeometry();
    }
  }
  if (!filenameChip3.empty()) {
    TFile algChipFile3(filenameChip3.c_str(), "read");
    std::vector<o2::detectors::AlignParam>* lAPchip3;
    algChipFile3.GetObject("ccdb_object", lAPchip3);
    for (int indice = 0; indice < lAPchip3->size(); indice++) {
      o2::detectors::AlignParam& algPar = lAPchip3->at(indice);
      algPar.applyToGeometry();
    }
  }
     */
  }

  //_____________________________________________________________________________
  void GetFileAlignedSurvey(Int_t level = 4) {
    /// Obtain the o2sim_geometry.root (corrected or pre-aligned geometry) using
    /// the given alignment file

    // Give the name of the stored alignment parameters for chips or disks
    // e.g. "vectorAlignParamsChip.root" or vectorAlignParamsDiskafterP2.root
    // and for disk: filenameDisk = "vectorAlgPar.root"

    std::string filenameDisk = "";
    std::string filenameChip = "";

    if (level == 2)
      filenameDisk = "MFTSurveyAlignParDisk.root";
    if (level == 4)
      filenameChip = "MFTSurveyAlignParChip.root";
    if (level == 10) {
      filenameChip = "MFTSurveyAlignParChip.root";
      filenameDisk = "MFTSurveyAlignParDisk.root";
    }

    MisAlignFromFile(true, filenameDisk, filenameChip);
    // MisAlignFromFile(true, filenameDisk, filenameChip, filenameChip2,
    // filenameChip3);

    if (!gGeoManager) {
      LOG(fatal) << "TGeoManager doesn't exist !";
      return;
    }

    // lock the geometry, or unlock with "false"
    gGeoManager->RefreshPhysicalNodes();

    // store the misaligned geometry in a root file
    gGeoManager->Export("geofiles/o2sim_geometry-survey-aligned.root");
  }

  /*
  //_____________________________________________________________________________
  void GenerateMisAlignment()
  {
    /// Obtain the o2sim_geometry.root after giving ramdom misalignment using
  the GeometryMisAligner

      o2::base::GeometryManager::loadGeometry("", false);

      // The misaligner
      o2::mft::GeometryMisAligner aGMA0;

    const std::string& ccdbHost = ""; // http://ccdb-test.cern.ch:8080
    long tmin = 0;
    long tmax = -1;
    const std::string& objectPath = "";
    double convMuToCm = 0.0001;
    double xladder = 2*convMuToCm;
    double yladder = 2*convMuToCm;
    double zladder = 2*convMuToCm;
    double xchip = 10*convMuToCm;
    double ychip = 20*convMuToCm;
    double zchip = 80*convMuToCm;
    double xdisk = 0.05;
    double ydisk = 0.05;
    double zdisk = -0.26;
    double pitchXdisk = 0.0001;
    double rollYdisk = 0.0012;
    double yawZdisk = 0.0012;
    double anglechip = 0.0001;

   if (!gGeoManager) {
     LOG(fatal) << "TGeoManager doesn't exist !";
     return;
   }

    gGeoManager->RefreshPhysicalNodes(false);
    aGMA0.SetDiskCartMisAlig(0.00, xdisk, 0.0, ydisk, zdisk, 0.05);
    aGMA0.SetSensorCartMisAlig(0.0, xchip, 0.0, ychip, 0.0, zchip);
    aGMA0.SetLadderCartMisAlig(0.0, xladder, 0.0, yladder, 0.0, zladder);
    //aGMA0.SetDiskAngMisAlig(0.0, pitchXdisk, rollYdisk, 0.0003, yawZdisk,
  0.0003);
    //aGMA0.SetSensorAngMisAlig(0.0, anglechip, 0.0, anglechip, 0.0, anglechip);

    aGMA0.MisAlign(true,  ccdbHost, tmin, tmax, objectPath, filename);

    // lock the geometry, or unlock with "false"
    gGeoManager->RefreshPhysicalNodes();

    // store the misaligned geometry in a root file
    gGeoManager->Export("o2sim_geometry-surveyaligned.root");

  }

   */
};

} // namespace mft
} // namespace o2

void loadMFTGeometry(
    Bool_t prealigned = false,
    const std::string &geoname = "geofiles/o2sim_geometry.root") {
  o2::base::GeometryManager::loadGeometry(geoname, prealigned);
  LOG(info) << " >> MFT geometry loaded";
}

void RunSurveyAlignWorkflow(Bool_t verbose = true, Int_t level = 4,
                            const std::string &ccdbHost = "", long tmin = 0,
                            long tmax = -1,
                            const std::string &objectPath = "") {

  /// Standard host name for ccdb http://ccdb-test.cern.ch:8080
  /// compute alignparam for alignable volume from survey values
  /// level 4 = chip
  /// level 2 = half-disk
  /// level 10 = both: half-disk and chip

  LOG(info) << " >>> RunSurveyAlignWorkflow ";
    
  loadMFTGeometry();

  o2::mft::SurveyAlignPar surAli;

  surAli.GetAlignParamsFromSurveyPositions(verbose, level, ccdbHost, tmin, tmax,
                                           objectPath);

  surAli.GetFileAlignedSurvey(level);
}
