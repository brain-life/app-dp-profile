#!/bin/bash

module load matlab

log=compiled/commit_ids.txt
true > $log
echo "/N/u/brlife/encode-dp" >> $log
(cd /N/u/brlife/git/encode && git log -1) >> $log
echo "/N/u/brlife/git/vistasoft" >> $log
(cd /N/u/brlife/git/vistasoft && git log -1) >> $log
#echo "/N/dc2/projects/lifebid/code/mba " >> $log
#(cd /N/dc2/projects/lifebid/code/mba && git log -1) >> $log
echo "/N/u/brlife/git/jsonlab" >> $log
(cd /N/u/brlife/git/jsonlab && git log -1) >> $log

cat > build.m <<END
addpath(genpath('/N/u/brlife/git/encode-dp'));
addpath(genpath('/N/u/brlife/git/vistasoft'));
addpath(genpath('/N/dc2/projects/lifebid/code/mba'))
addpath(genpath('/N/u/brlife/git/jsonlab'));

mcc -m -R -nodisplay -a /N/u/brlife/git/vistasoft/mrAnatomy/Segment -d compiled compute_profiles
exit
END
matlab -nodisplay -nosplash -r build && rm build.m


