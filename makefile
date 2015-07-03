# makefile for ADMM decoder
ldpc.out: LDPC_Simulator_Example.cpp
	g++ -O3 -o ldpc.out  LDPC_Class.cpp LDPC_Simulator_Example.cpp