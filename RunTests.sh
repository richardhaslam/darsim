# runtests.sh
cd IntegratedTests

LOGFILE=log.txt

matlab -nodesktop -nosplash -minimize -wait -logfile "$LOGFILE" -r 'RunTests';
CODE=$?

cat "$LOGFILE"

exit $CODE



