
echo "Running basic script"

## script-specific info - setup .h files correctly
mkdir -p class_def/run2
cp *.h class_def/run2/
cp GenParticles.h class_def/

#untar the input filelists
tar xf input_files.tar.gz

ss="${2#*=}"

echo "Running script: $ss"

root -b -l -q $ss+\(1,$1\)

echo "getting dummy job report"
echo "================= Dumping Input files ===================="
python -c "import PSet; print '\n'.join(list(PSet.process.source.fileNames))"
#cmsRun -j FrameworkJobReport.xml -p PSet.py
echo "Done!"
