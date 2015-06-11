rootcint -f linkdefs.cxx -c `root-config --cflags` defs/linkdefs.h

g++ -Wall linkdefs.cxx -I$ROOTINC `root-config --libs --glibs` class_def/JetAna.C class_def/Tracks.C class_def/HLT.C class_def/HiTree.C class_def/Skim.C class_def/GenParticles.C class_def/pfcand.C make_ntuples.C -o HT_make_ntuples.out
