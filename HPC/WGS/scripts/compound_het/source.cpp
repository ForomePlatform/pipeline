#include<iostream>
#include<vector>
#include<list>
#include<map>
#include<sstream>
#include<string>
#include<fstream>
#include<stdexcept>
#include<cstdlib>
#include<ctime>
#include<cmath>
#include<limits>
#include<algorithm>
#include<boost/math/special_functions/gamma.hpp>
#include"decs_1.h"
#include"tests_1.h"


// monty



int main(int argc,char** argv){
        using namespace std;

        {
                //(infilename,trio_filename,outfilename,P_T,ExAC_AF_T,AF_unrel_T)
                string infilename,trio_filename,outfilename,P_T_str,ExAC_AF_T_str,AF_unrel_T_str;
                for(int i(0);i<argc;++i){
                        if(string(argv[i])=="-I")
                                infilename=argv[i+1];
                        if(string(argv[i])=="-TR")
                                trio_filename=argv[i+1];
                        if(string(argv[i])=="-O")
                                outfilename=argv[i+1];
                        if(string(argv[i])=="-T")
                                P_T_str=argv[i+1];
                        if(string(argv[i])=="-E")
                                ExAC_AF_T_str=argv[i+1];
                        if(string(argv[i])=="-A")
                                AF_unrel_T_str=argv[i+1];
                }
                ////
                const double P_T(atof(P_T_str.c_str()));
                const double ExAC_AF_T(atof(ExAC_AF_T_str.c_str()));
                const double AF_unrel_T(atof(AF_unrel_T_str.c_str()));

                cout<<"infilename="<<infilename<<"\n";
                cout<<"trio_filename="<<trio_filename<<"\n";
                cout<<"outfilename="<<outfilename<<"\n";
                cout<<"P_T="<<P_T<<"\n";
                cout<<"ExAC_AF_T="<<ExAC_AF_T<<"\n";
                cout<<"AF_unrel_T="<<AF_unrel_T<<"\n";

                worker(infilename,trio_filename,outfilename,P_T,ExAC_AF_T,AF_unrel_T);
                //tester_3(infilename,trio_filename);


        }


        //tester_4();

        return 0;
}
