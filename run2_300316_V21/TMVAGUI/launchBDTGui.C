#include "TMVAGui.C"

void launchBDTGui() {
// Launch the GUI for the root macros
//source /afs/cern.ch/sw/lcg/external/gcc/4.3.2/x86_64-slc5/setup.csh
//source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.32.02/x86_64-slc5-gcc43-opt/root/bin/thisroot.csh
// gui not working in root6

  gROOT->LoadMacro("TMVAGui.C");
  TMVAGui("../TMVA_ZnunuHighPt.root");
}
