!lscpu |grep 'Model name'
!ulimit -s 8192
%cd /content
!git clone https://github.com/DeltaGroupNJUPT/Vina-GPU.git

%cd /usr/local
!wget https://boostorg.jfrog.io/artifactory/main/release/1.80.0/source/boost_1_80_0.tar.gz
!chmod 777 boost_1_80_0.tar.gz
!tar -xzvf boost_1_80_0.tar.gz


%cd /content/Vina-GPU
import os

# Read in the file
with open('Makefile', 'r') as file :
  filedata = file.read()

# Replace the target string
filedata = filedata.replace('../boost_1_77_0', '/usr/local/boost_1_80_0')
filedata = filedata.replace('-DOPENCL_3_0', '-DOPENCL_1_2')

# Write the file out again
with open('Makefile', 'w') as file:
  file.write(filedata)

%cd /content/Vina-GPU
!make clean
!make source
!./Vina-GPU --config ./input_file_example/2bm2_config.txt
# Remake to build without compiling kernel files
!make clean
!make