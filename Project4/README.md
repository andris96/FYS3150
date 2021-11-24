Welcome to Project 4. Below are the instructions for running each program related to this project.

To run the c++ program (remove -fopenmp if you want to run serially) :

To build:
g++ main.cpp src/IsingModel.cpp src/TestIsingModel.cpp -I include -o main.exe -larmadillo -fopenmp

To run:
./main.exe

To run the python program: 

python3 plot.py
