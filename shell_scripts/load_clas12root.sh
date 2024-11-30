echo "loading clas12root"
module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
module purge
module load sqlite/dev
module load clas12
echo "done loading clas12root."

echo "git pull && macros/execute_skim"
git pull && macros/execute_skim
