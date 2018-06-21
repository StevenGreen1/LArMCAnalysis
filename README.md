Commands:
Build local larsoft as normal 
cd srcs
git clone https://github.com/StevenGreen1/LArMCAnalysis.git
mrb uc
cd $MRB_BUILDDIR
mrbsetenv
mrb i -j4
mrbslp
