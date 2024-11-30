echo "loading clas12root"
module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
module purge
module load sqlite/dev
module load clas12
echo "done loading clas12root."
echo 
echo "git pull && python macros/execute_skim"
git pull && python macros/execute_skim
