#ifndef DP_OSCILLATIONHELPER_HXX_SEEN
#define DP_OSCILLATIONHELPER_HXX_SEEN

#include "TH1D.h"

class BargerPropagator;

struct OscillationHelper {
  enum NuTypes {
    kNuebarType = -1,
    kNumubarType = -2,
    kNutaubarType = -3,
    kNueType = 1,
    kNumuType = 2,
    kNutauType = 3,
  };

  double DipAngle_degrees; // = 5.8;
  double OscParams[6];     // = {0.825, 0.10, 1.0, 7.9e-5, 2.5e-3, 0.0};
  double LengthParam;      // = 0xdeadbeef;

  bool IsSetUp;

  BargerPropagator* bp;
  Int_t FromPDG, ToPDG;
  NuTypes FromType, ToType;

  NuTypes GetNuType(int pdg);

  void Setup_dipangle(double OscParams[6], double DipAngle_degrees = 5.8);
  void Setup_baseline(double OscParams[6], double baseline_km = 295);

  OscillationHelper() : IsSetUp(false), bp(NULL){};
  OscillationHelper(OscillationHelper const &other);
  ~OscillationHelper();

  void SetOscillationChannel(int PDGFrom, int PDGTo);
  double GetWeight(double ENu_GeV);
  void OscillateHistogram(TH1 *h);
};

#endif
