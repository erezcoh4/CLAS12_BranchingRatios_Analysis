echo "loading clas12root"
module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
module purge
module load sqlite/dev
module load clas12
echo "done loading clas12root."
echo 
echo "git pull && python ./macros/execute_skim --fdebug=5 --no-email --FirstEvent=0 --NeventsMax=69"
git pull && python ./macros/execute_skim --fdebug=5 --no-email --FirstEvent=0 --NeventsMax=69

echo "Done."
echo
echo "Examples:"
echo "git pull && python ./macros/execute_skim --fdebug=5 --no-email --FirstEvent=0 --NeventsMax=69"
echo "git pull && python ./macros/execute_skim --fdebug=1 --NeventsMax=-1"
echo "python ./macros/execute_skim --fdebug=1 --NeventsMax=-1 --Nruns=10"
echo
