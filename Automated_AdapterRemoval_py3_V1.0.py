

"""
Created on Mon Jan 18 17:41:10 2021

@author: minkj
"""
#________________ option parse _______________________________
import sys 

args = sys.argv[1:]

option_dict={'-cores':"32",'-paired':"1", '-exclude':""}

help="""

________________________________________________________________________________________________

AdapterRemoval is a tool to remove adapter and trim low-quality reads
This pipeline needs the following prerequisite software installed:
AdapterRemoval: conda install bioconda::adapterremoval (Anaconda or Miniconda is needed.)
________________________________________________________________________________________________

Run example: python Automated_AdapterRemoval_py3_V1.0.py -include .fq -exclude filtered -paired 1
________________________________________________________________________________________________

Usage;

-help       show the option list
-include    keywords for filtering files to include 
-exclude    keywords for filtering files to exclude (default is "")
-paired     1: paired-end, 2: single-end (default is 1)
-cores      number of cores for AdapterRemoval (default is 32)
________________________________________________________________________________________________
"""
if args==[]:
    print (help)
    quit()

for i in range(len(args)):
    if args[i].startswith("-"):
        try:
            option_dict[args[i]]=args[i+1]
        except:
            if args[0]=="-help":
                print (help)

                quit()


def file_list(infilter, exfilter): #inputs are string
    import os
    f_list=os.listdir('.')
    infiltered=list(filter(lambda x: infilter in x, f_list))
    if option_dict["-exclude"]=="":
        exfiltered=infiltered
    else:
        exfiltered=list(filter(lambda x: exfilter not in x, infiltered))
    exfiltered.sort()
    print ("\n\nFiltered file list:\n", exfiltered)
    return exfiltered

def seq_pairing(filtered_files): # input is list
    n=0
    paired_list=[]
    infile=filtered_files
    #infile.sort()
    #print (infile)
    infile_len=len(infile)
    for k in range(int(len(infile)/2)):
        #print (n)
        f1=infile[n]
        f2=infile[n+1]
        n=n+2
        
        f1_list=list(f1)
        f2_list=list(f2)
        if len(f1_list)==len(f2_list):
            #diff=0
            for k in range(len(f1_list)):
                if f1_list[k]!=f2_list[k]:
                    if f1_list[k] in ["1","2"] and f2_list[k] in ["1","2"]:
                        diff="paired"
                    else:
                        diff="error"
            if diff=="paired":
                paired_list.append((f1,f2))
            else:
                print ("\n\n", f1,"    and    ", f2,"   are compared.")
                print ("Error occurred in pairing.. \nF and R names must be same except for the distinguishing number, 1 and 2.")
                quit()

        else:
            print ("\n\n", f1,"    and    ", f2,"   are compared.")
            print ("Error occurred in pairing.. \nF and R names must be same except for the distinguishing number, 1 and 2.")
            quit()

    return paired_list

def main():
    infilter_cont=option_dict["-include"]
    exfilter_cont=option_dict["-exclude"]

    filtered=file_list(infilter_cont,exfilter_cont)
    import os
    if option_dict['-paired']=="1":
        paired=seq_pairing(filtered)
        n=0
        for tuple_file in paired:
            n=n+1
            print ("\n\n"+str(n)+"/"+str(len(paired)), tuple_file, "<=== AdapterRemoval is running..")
            run_fastqc=os.system("AdapterRemoval --threads %s --file1 %s --file2 %s --output1 %s --output2 %s " % (option_dict['-cores'], tuple_file[0], tuple_file[1],tuple_file[0]+"_1_filtered.fq", tuple_file[1]+"_2_filtered.fq"))

    elif option_dict['-paired']=="2":
        #filtered.sort()
        n=0
        for file in filtered:
            n=n+1
            print ("\n\n"+str(n)+"/"+str(len(filtered)), file, "<=== AdapterRemoval is running..") 
            run_fastqc=os.system("AdapterRemoval --threads %s --file1 %s --output1 %s" % (option_dict['-cores'], file, file+"_filtered.fq"))

    print ("Completed!!")

if __name__=="__main__":
    main()




