# runtests.sh

LOGFILE=log.txt

matlab -nodesktop -nosplash -logfile "$LOGFILE" -r 'CITests/RunTests';
CODE=$?

cat "$LOGFILE"

exit $CODE



