for i in 0.001  0.01  0.1  1  10  100 ; do cp mean.py $i; cd $i;tail -n 20  EdgeInfo> nn; echo $i;python mean.py ; cd ..; done
