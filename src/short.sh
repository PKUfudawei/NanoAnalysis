#!/bin/bash 
# short.sh: a short discovery job 

printf "Start time: "; /bin/date 
printf "Job is running on node: "; /bin/hostname 
printf "Job running as user: "; /usr/bin/id 
printf "Job is running in directory: "; /bin/pwd 

echo
echo "Working hard..."
sleep ${1-15}
echo "Sleep for ${1} seconds"
echo "Science complete!"
