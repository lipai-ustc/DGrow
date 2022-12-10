../../teaks.x infile 
cp EdgeInfo nn
sed -i '1,3d' nn
python mean.py
