# Parker-jet-code
 Parker's code for the jet collectivity study

Code to make trees -> local or CRAB, produces roots

Macros in src/
Main macro( list of files, job number 1-N)
    loads root file
    2PC, jet multiplicity, save histograms
    Line 136-137 comment out
        MC corrections
    Line 141-144 comment out
        HLT efficiency

    Line 79 
    Line 149 Main code starts



Header file
    include/
    3 header files: coordinate tools, constants, Tree details
    
Once changes are made to macro and header, back up one directory so you can see both src and include âmake cleanâ removes old executable. âMakeâ compiles new changes (very fast, not like scram b)

Make file
Change the name to reflect the macro you are working on, (change it in 5 locations)

Produces executable in âbinâ directory

Data file lists
    data_v2_list.txt, contains all paths to all files


Executable
Condor jobs in batch
    divides, gives args, etc
