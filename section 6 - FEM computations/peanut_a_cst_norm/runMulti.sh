gmsh disk_T_DtN.geo;

python3 sample_k.py -1.1;
python3 sample_k.py -0.9;

make;

./exec-x86_64-linux-g++-9-Release -1.2 -1.1 est;
./exec-x86_64-linux-g++-9-Release -0.9 -0.8 est;
