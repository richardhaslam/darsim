# runtests.sh
cd CITests

LOGFILE=log.txt

matlab -nodesktop -nosplash -logfile "$LOGFILE" -r 'RunTests';
CODE=$?

cat "$LOGFILE"

exit $CODE



