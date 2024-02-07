# /opt/homebrew/Cellar/open-mpi/4.1.5/bin/mpicxx -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -I /opt/homebrew/Cellar/open-mpi/4.1.5/include/ -o SampleSort_radixSort ./SampleSort_radixSort.cpp
# /opt/homebrew/Cellar/open-mpi/4.1.5/bin/mpirun -np 3 ./SampleSort_radixSort randomString_50_alpha_26_set_1.txt

g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ rlapp_tokenBased.cpp -o rlapp_tokenBased
./rlapp_tokenBased ds6_800k

# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ rlapp_tokenBased.cpp -o rlapp_tokenBased
# ./rlapp_tokenBased ds7_1M_fl
# ./rlapp_tokenBased ds_test
# ds7_1M_fl

# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -I /usr/local/  -o kwayMergeSort kwayMergeSort.cpp
# ./kwayMergeSort ds11_5M_fl

# ds11_5M_fl

# /opt/homebrew/Cellar/open-mpi/5.0.0/bin/mpicxx -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -I /opt/homebrew/Cellar/open-mpi/4.1.5/include/ prla_dynamic_hybrid.cpp -o prla_dynamic_hybrid
# /opt/homebrew/Cellar/open-mpi/5.0.0/bin/mpirun -np 1 ./prla_dynamic_hybrid ds1_50k_fl
#  --map-by node:PE=3 --bind-to core
# ds11_5M_fl

# /opt/homebrew/Cellar/open-mpi/5.0.0/bin/mpicxx -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -I /opt/homebrew/Cellar/open-mpi/4.1.5/include/ testHybrid.cpp -o testHybrid
# /opt/homebrew/Cellar/open-mpi/5.0.0/bin/mpirun -np 4 ./testHybrid
