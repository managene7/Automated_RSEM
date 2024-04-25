

"""
Created on Mon Jan 18 17:41:10 2021

@author: minkj
"""
#________________ option parse _______________________________
import sys 

args = sys.argv[1:]

option_dict={'-cores':"32",'-paired':"1", '-exclude':"",'-ref':"", '-target':"1", '-out':"RSEM_results", '-build_index':"2", '-parsing_only':"2", "-skip_filtering":""}
help="""
________________________________________________________________________________________________

RSEM is a pipeline to calculate counts, TPM, and FPKM from genes or transcripts.
This pipeline fully automatizes the RSEM from the build index for mapping to parsing RSEM results.
This pipeline needs the following prerequisite software installed:

AdapterRemoval: conda install bioconda::adapterremoval

RSEM (pipeline): sudo apt-get update -y
                 sudo apt-get install -y rsem

Bowtie2 (mapping) : sudo apt-get install -y bowtie2

Samtools: sudo apt-get install -y samtools
________________________________________________________________________________________________

The pipeline consists of building an index, running RSEM, and parsing the RSEM result
You can run this pipeline by following methods:
1) run all steps: use all options
   example: python Automated_RSEM-Bowtie2_multi-processing_py3_v2.0.py -skip_filtering 2 \\
            -build_index 1 -include .fq -exclude .fa -paired 1 -parsing_only 2 -target 1

2) skip AdapterRemoval and run all the rest of the steps: use all options
   example: python Automated_RSEM-Bowtie2_multi-processing_py3_v2.0.py -skip_filtering 1 \\
            -build_index 1 -include filtered.fq -exclude .fa -paired 1 -parsing_only 2 -target 1
            
3) run RSEM and parse RSEM results: set '-build_index' option as '2' and use all options
   example: python Automated_RSEM-Bowtie2_multi-processing_py3_v2.0.py -skip_filtering 1 \\
            -build_index 2 -ref genome.fa -include .fq -exclude .fa -paired 1 -parsing_only 2 \\
            -target 1

4) only parse RSEM results: set '-parsing_only' option as '1' and use '-target' option.
   example: python Automated_RSEM-Bowtie2_multi-processing_py3_v2.0.py -skip_filtering 1 \\
            -parsing_only 1 -target 1

________________________________________________________________________________________________

Usage;

-help                       show options
-skip_filtering (required)  1: yes (skip filtering), 2: no (run AdapterRemoval)
-build_index    (option)    1: yes, 2: no (default is 2)
-ref            (required)  ref seq file name (index file name must be same with ref seq name)
-include        (required)  key letters to select RNA-seqs among files (example: '.fastq')
-exclude        (option)    key letters to exclude the rest of files (default is "")
-paired         (option)    1: paired-end, 2: single-end (default is 1)
-parsing_only   (option)    1: parsing the RSEM results only, 2: run all steps (default is 2)
-target         (option)    1: representative genes, 2: isoforms (default is 1)
-out            (option)    tag for the output file name (default is RSEM_result)                       
-cores          (option)    number of cores for RSEM (default is 32)
______________________________________________________________________________________________
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
    for k in range(int(len(infile)/2)):
        f1=infile[n]
        f2=infile[n+1]
        n=n+2

        f1_list=list(f1)
        f2_list=list(f2)
        if len(f1_list)==len(f2_list):
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

def rsem_pair(thread, tuple_file, ref):
    import os
    run_RSEM=os.system("rsem-calculate-expression -p %s --paired-end --bowtie2 --estimate-rspd --append-names %s %s %s %s" % (thread, tuple_file[0], tuple_file[1], ref, tuple_file[0]+"--"+tuple_file[1]+"_rsem_output.txt"))    
    if run_RSEM !=0:
        print ("\n\nRunning RSEM raised error!\nCheck input files and tryagain.\n\n")
        quit()

def rsem_single(thread, file, ref):
    import os
    run_RSEM=os.system("rsem-calculate-expression -p %s --bowtie2 --estimate-rspd --append-names %s %s %s" % (thread, file, ref, file+"_rsem_output.txt"))
    if run_RSEM !=0:
        print ("\n\nRunning RSEM raised error!\nCheck input files and tryagain.\n\n")
        quit()

def build_index(ref_genome, type, gene_info, gene_info_type):
    import os
    print ("\n\nStart to build index file for Bowtie2 mapping..\n\n")
    if type=="1":
        if gene_info_type=="1":
            build_index=os.system(f"rsem-prepare-reference --gtf {gene_info} --bowtie2 -p 32 {ref_genome} {ref_genome}")
            if build_index !=0:
                print ("\n\nBuilding index failed.\nCheck the input sequence file.")
                quit()
        elif gene_info_type=="2":
            build_index=os.system(f"rsem-prepare-reference --gff3 {gene_info} --bowtie2 -p 32 {ref_genome} {ref_genome}")
            if build_index !=0:
                print ("\n\nBuilding index failed.\nCheck the input sequence file.")
                quit()
        
    elif type=="2":
        transcript_to_gene_map=ref_genome+"_transcript_to_gene_map.txt"
        build_index=os.system(f"rsem-prepare-reference --bowtie2 -p 32 --transcript-to-gene-map {transcript_to_gene_map} {ref_genome} {ref_genome}")
        if build_index !=0:
                print ("\n\nBuilding index failed.\nCheck the input sequence file.")
                quit()

    print ("\n\nBuilding index file has completed!")


def fasta_reformat(fasta):
    print ("Converting CDS sequence..\n\n")
    infile=open(fasta,'r')
    infile_lines=infile.readlines()
    out_file=open(fasta,'w')
    init=0
    for line in infile_lines:
        line=line.strip()

        if ">" in line:
            if init==0:
                name=line.split()[0]
                seq_list=[]
                init=1
            else:
                seqs="".join(seq_list)
                out_file.write(name+"\n")
                out_file.write(seqs+"\n")
                name=line.split()[0]
                seq_list=[]
        else:
            seq_list.append(line)

    seqs="".join(seq_list)
    out_file.write(name+"\n")
    out_file.write(seqs+"\n")


def main():
    import os

    if option_dict['-skip_filtering']=="":
        print ("\n\nChoose '-skip_filtering' option and try again. (1: skip filtering, 2: run AdapterRemoval)\n\n" )
        quit()


    if option_dict['-parsing_only'] !="1":

        #___ build index file for bowtie2_____________________________
        if option_dict['-build_index']=="1":
            target_seq=option_dict['-ref']
            
            while 1:
                seq_type=input("Fasta file contents | 1: genome sequences, 2: gene sequences | : ")
                if seq_type in ["1","2"]:
                    if seq_type =="2":
                        reformat=fasta_reformat(target_seq)

                        gene_map_name=target_seq+"_transcript_to_gene_map.txt"
                        transcript_to_gene_map=os.system(f"extract-transcript-to-gene-map-from-trinity {target_seq} {gene_map_name}")
                    break
                else:
                    print ("\n\nChoose the type by entering 1 or 2.\n\n")

            if seq_type=="1":
                while 1:
                    gene_info=input("Enter the name of gtf or gff3 file: ")
                    if gene_info!="":
                        
                        break
                    else:
                        print ("\n\nThe name is missing.. Try it again.\n\n")
                while 1: 
                    gene_info_type=input("Choose the type | 1: gtf, 2: gff3 |: ")
                    if gene_info_type in ["1","2"]:
                        break
                    else:
                        print ("\n\nChoose the type by entering 1 or 2.\n\n")
            else:
                gene_info=""
                gene_info_type=""
            build_index(target_seq, seq_type, gene_info, gene_info_type)
            
        #__________________________________________________________________
            


        #___ run AdapterRemoval____________________________________
        if option_dict['-skip_filtering']!="1":

            infilter_cont=option_dict["-include"]
            exfilter_cont=option_dict["-exclude"]

            filtered=file_list(infilter_cont,exfilter_cont)
            if option_dict['-paired']=="1":
                paired=seq_pairing(filtered)
                n=0
                for tuple_file in paired:
                    n=n+1
                    print ("\n\n"+str(n)+"/"+str(len(paired)), tuple_file, "<=== AdapterRemoval is running..\n\n")
                    run_AdapterRemoval=os.system("AdapterRemoval --threads %s --file1 %s --file2 %s --output1 %s --output2 %s " % (option_dict['-cores'], tuple_file[0], tuple_file[1],tuple_file[0]+"_1_filtered.fq", tuple_file[1]+"_2_filtered.fq"))
                    if run_AdapterRemoval!=0:
                        print ("\n\nRunning AdapterRemoval Raised an error.\nCheck input files and try again.\n\n")
                        quit()
                    else:
                        print ("\n\nAdaptRemoval has completed!\n\n")

            elif option_dict['-paired']=="2":
                #filtered.sort()
                n=0
                for file in filtered:
                    n=n+1
                    print ("\n\n"+str(n)+"/"+str(len(filtered)), file, "<=== AdapterRemoval is running..\n\n") 
                    run_AdapterRemoval=os.system("AdapterRemoval --threads %s --file1 %s --output1 %s" % (option_dict['-cores'], file, file+"_filtered.fq"))
                    if run_AdapterRemoval!=0:
                        print ("\n\nRunning AdapterRemoval Raised an error.\nCheck input files and try again.\n\n")
                        quit()
                    else:
                        print ("\n\nAdaptRemoval has completed!\n\n")
            infilter_cont="_filtered.fq"
        #_________________________________________________________


        #___ RSEM mapping _________________________

        filtered=file_list(infilter_cont,exfilter_cont)

        import multiprocessing
        if option_dict['-paired']=="1":
            paired=seq_pairing(filtered)
            len_paired=len(paired)
            num_cores=int(option_dict['-cores'])

            num_iter=int(len_paired/num_cores)+1
            posi=0
            multiprocessing.freeze_support()
            for i in range(num_iter):
                sub_list=paired[posi:posi+num_cores]
                posi+=num_cores
                len_sublist=len(sub_list)

                #___allocate threads per process when number of seqs are less than threads___
                thread=int(num_cores/len_sublist)
                if thread <1:
                    thread=1
                #__________________________________________

                sub_args=[]
                for id in sub_list:
                    sub_args.append((thread, id, option_dict['-ref']))

                print ("\n\nThis subsets are in processing..\n\n",sub_list,"\n\n")
                pool=multiprocessing.Pool(processes=num_cores)
                pool.starmap(rsem_pair, sub_args)
                pool.close()
                pool.join()

        elif option_dict['-paired']=="2":
            len_filtered=len(filtered)
            #___allocate threads per process when number of seqs are less than threads___
            thread=int(int(option_dict['-cores'])/len_filtered)
            if thread <1:
                thread=1
            #__________________________________________

            num_cores=int(option_dict['-cores'])

            num_iter=int(len_filtered/num_cores)+1
            posi=0

            multiprocessing.freeze_support()
            for i in range(num_iter):
                sub_list=filtered[posi:posi+num_cores]
                posi+=num_cores
                len_sublist=len(sub_list)

                #___allocate threads per process when number of seqs are less than threads___
                thread=int(num_cores/len_sublist)
                if thread <1:
                    thread=1
                #__________________________________________

                sub_args=[]
                for id in sub_list:
                    sub_args.append((thread, id, option_dict['-ref']))

                print ("\n\nThis subsets are in processing..\n\n",sub_list,"\n\n")
                pool=multiprocessing.Pool(processes=num_cores)
                pool.starmap(rsem_single, sub_args)
                pool.close()
                pool.join()

        print ("RSEM running has completed. \n\nParsing of the RSEM results starts.")
    
    
    #_________RSEM result parsing___________

    if option_dict['-target']=="1":
        infilter_cont="genes.results"
        exfilter_cont="isoforms.results"
        option_dict['-out']=option_dict['-out']+"_genes"
    elif option_dict['-target']=="2":
        infilter_cont="isoforms.results"
        exfilter_cont="genes.results"
        option_dict['-out']=option_dict['-out']+"_isoforms"
    else:
        print (f"\n\nYou entered -target option as {option_dict['-target']}. \nSo, it will be ignored, and isoforms will be parsed.")

    filtered=file_list(infilter_cont,exfilter_cont)

    import csv
    n=0

    seq_list=[]
    init=0
    cont_dic={} # {transcript:[count,tpm, fpkm]}
    for f in filtered:
        n=n+1
        f_open=open(f,'r')
        sub_cont_dic={}
        transcript_list=[]
        while 1:
            line=f_open.readline().strip()
            if line=="":
                break
            else:
                line=line.split()
                if line[0]!="transcript_id":
                    transcript=line[0].split("_")[0]
                    # if transcript not in ['gene', 'transcript']:
                        
                    
                    if transcript not in ['gene', 'transcript']:
                        transcript_list.append(transcript)
                        sub_cont_dic[transcript]=[line[4],line[5],line[6]]
        init=1
        #seq_name=f.split("_")[0]+"_"+f.split("_")[1]+"_"+f.split("_")[2]
        seq_name=f.split(".f")[0]#+"_"+f.split("_")[1]+"_"+f.split("_")[2]
        seq_list.append(seq_name)
        cont_dic[seq_name]=sub_cont_dic
    seq_list.sort()


    count_dic={}
    tpm_dic={}
    fpkm_dic={}

    for t_id in transcript_list:
        count_dic[t_id]=[]
        tpm_dic[t_id]=[]
        fpkm_dic[t_id]=[]
        for seq in seq_list:
            #print (t_id)
            # print (cont_dic[seq][t_id])
            count_dic[t_id].append(str(int(float(cont_dic[seq][t_id][0]))))
            tpm_dic[t_id].append(str(int(float(cont_dic[seq][t_id][1]))))
            fpkm_dic[t_id].append(str(int(float(cont_dic[seq][t_id][2]))))
    count_csv=csv.writer(open(option_dict['-out']+"_Count_all.csv",'w', newline=""))
    TPM_csv=csv.writer(open(option_dict['-out']+"_TPM_all.csv",'w', newline=""))
    FPKM_csv=csv.writer(open(option_dict['-out']+"_FPKM_all.csv",'w', newline=""))

    #print (count_dic)
    #print (list(["Transcript ID"])+seq_list)
    count_csv.writerow(list(["Gene_ID"])+seq_list)
    TPM_csv.writerow(list(["Gene_ID"])+seq_list)
    FPKM_csv.writerow(list(["Gene_ID"])+seq_list)

    for transcript_id in transcript_list:
        count_csv.writerow(list([transcript_id])+count_dic[transcript_id])
        TPM_csv.writerow(list([transcript_id])+tpm_dic[transcript_id])
        FPKM_csv.writerow(list([transcript_id])+fpkm_dic[transcript_id])

    print ("\n\nRSEM result parsing has completed!!")

if __name__ == '__main__':
    main()    
