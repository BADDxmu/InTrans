# -*- coding: utf-8 -*-
import configparser
import os,sys,shutil
sys.path.append('..')
import logging
import time
import subprocess

logging.basicConfig(level=logging.DEBUG,
                format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                datefmt='%a, %d %b %Y %H:%M:%S',
                filename='./log',
                filemode='w')
cfg=sys.argv[1]
__config__ = configparser.ConfigParser()
__config__.read(cfg,encoding="utf-8")

def loggingtime(funcname,starttime):
    logging.info(funcname + ' Time Running %f s' %(time.time()-starttime))
    global start
    start = time.time()

# Check if the folder exists
def CheckFileDirs():
    try:
        if not os.path.exists('./output'):
            os.mkdir('./output')
        if not os.path.exists('./temp_output'):
            os.mkdir('./temp_output')
    except Exception as e:
        logging.error("CheckFileDirs Function Wrong",e)
        print("CheckFileDirs Wrong",e)

# find the first different word's index between str1 and str2
def FindStrIndex(str1,str2):
    length = len(str1)
    if len(str1) > len(str2): length = len(str2)
    for i in range(length):
        if str1[i] != str2[i]:
            break
    return i

# IDBA find the best file to use
def IDBAFindBestFile(filedir):
    if not os.path.exists(filedir):
        print("Didn't find the file, Check '" + filedir + "/' has the log file")
        return ''
    try:
        BestKmer = {}; TempKmer = ''
        logfile = filedir + '/log'
        fp = open(logfile,'r')
        line = fp.readline()
        i = 0
        while line:
            if 'kmer' in line and "kmers" not in line:
                TempKmer = line.replace('kmer','').strip()
                BestKmer[TempKmer] = ''
            if 'contigs' in line and 'length' in line:
                i += 1
            if i == 3 and 'contigs' in line:
                if line.find('n50') != -1:
                    tempN50 = line[line.find('n50')+4:line.find('max')-1]
                    BestKmer[TempKmer] = int(tempN50.strip())
                i = 0
            line = fp.readline()
        BestKmerNum = sorted(BestKmer.items(), key=lambda item:item[1], reverse=True)[0][0]
        fp.close()
        # print(BestKmerNum)
        logging.info("for " + filedir +" we seletced the de novo assemble result when kmer=" + str(BestKmerNum) + " for the biggest N50 we got with this kmer parameter")
        return 'transcript-' + BestKmerNum + '.fa'
    except Exception as e:
        logging.error("IDBAFindBestFile Function Wrong", e)
        print("IDBAFindBestFile Wrong",e)
        exit(0)

# IDBA Function #
def IDBA():
    logging.info("Start Running IDBA:")
    try:
        path = __config__.get("seq data","fastq_dir")
        if path is None or not os.path.exists(path):
            logging.error("Run.cfg didn't set right fastq_dir, can't find your path")
            print("Run.cfg didn't set right fastq_dir")
            exit(0)
        if path[-1] != '/':
            path = path + '/'
        try:
            Pairsamplelist = eval(__config__.get('seq data', 'fastq_files'))
        except Exception as e:
            logging.error("fastq_files format error!",e)
            print("fastq_files format error!")

        mink = __config__.get("IDBA","mink")
        maxk = __config__.get("IDBA","maxk")
        step = __config__.get("IDBA","step")
        num_threads = __config__.get("IDBA","num_threads")

        IDBABestFiles = []
        for sample in Pairsamplelist:

            IDBACombineStr = "fq2fa --merge "

            samplelist = sample.strip().split(" ")
            # print(samplelist)
            for sam in samplelist:

                if not os.path.exists(path + sam):
                    logging.error("Didn't find your file. Check your fastq_files's format is right.")
                    print("Didn't find your file. Check your fastq_files's format is right.")
                    exit(0)

                IDBACombineStr = IDBACombineStr + str(path + sam) + ' '

            if len(samplelist) > 1:
                sampledir = './temp_output/' + samplelist[0][0:FindStrIndex(samplelist[0],samplelist[1])-1] + '.fa'
            else:
                sampledir = './temp_output/' + samplelist[0] + '.fa'

            IDBACombineStr = IDBACombineStr + str(sampledir) + ' '
            # print(IDBACombineStr)
            try:
                subprocess.call(IDBACombineStr,shell=True)
                # pp.communicate()
                # os.system(IDBACombineStr)
            except subprocess.CalledProcessError as e:
                logging.error("IDBA first Command Line is wrong.",e.output,e.returncode)
                exit(0)

            tempdir = sampledir.replace('.fa','')
            IDBAExecute = "idba_tran -r &sample -o &output --mink &mink --maxk &maxk --step &step --num_threads &num_threads".replace("&sample",sampledir).replace("&output",tempdir).replace("&mink",mink).replace("&maxk",maxk).replace("&step",step).replace("&num_threads",num_threads)
            # print(IDBAExecute)
            try:
                os.system(IDBAExecute)
            except Exception as e:
                logging.error("IDBA second Command line is wrong.")
                exit(0)

            # print(tempdir)
            IDBABestFiles.append(tempdir + '/' + IDBAFindBestFile(tempdir))
        # print(IDBABestFiles)

        logging.info("IDBA Completed")
        return IDBABestFiles
    except Exception as e:
        logging.error("IDBA Function Wrong:" + str(e))
        print("IDBA Wrong",e)
        exit(0)

# change name of HeteroFiles
def HeteroNewFiles(heterodir):
    if not os.path.exists(heterodir):
        logging.error("Didn't find the heterodir file")
        exit(0)
    try:
        newHeteroFiles = './temp_output/HeteroFile/' + heterodir.split('/')[-1]
        if not os.path.exists("./temp_output/HeteroFile/"):
            os.makedirs('./temp_output/HeteroFile/')

        rb = open(heterodir,'r')
        wb = open(newHeteroFiles, 'w')
        rline = rb.readline()
        while rline:
            if '>' in rline:
                rline = rline[:rline.find('>') + 1] + " Hetero " + rline [rline.find('>') + 1:]
            wb.write(rline)
            rline = rb.readline()
        wb.close()
        rb.close()
        return newHeteroFiles
    except Exception as e:
        logging.error("HeteroNewFiles Function Wrong", e)
        print("HeteroNewFiles Wrong:" + str(e))
        exit(0)

# CAT Function
def CATCombine(IDBABestList):
    logging.info("Start Running CAT:")
    if IDBABestFiles is None or (type(IDBABestList[0]) is str and '.fa' not in IDBABestList[0]):
        print("Can't find the IDBA files")
        logging.error("Can't find the IDBA files")
        exit(0)
    try:
        heter_prefix_path = __config__.get("seq data","heter_prefix")

        # define heter_file
        is_heterfile = True
        if heter_prefix_path is None or len(heter_prefix_path.strip(' ')) == 0:
            logging.warning("Run.cfg didn't set heter_prefix")
            is_heterfile = False
        if len(heter_prefix_path) !=0 and heter_prefix_path[-1] != '/':
            heter_prefix_path = heter_prefix_path + '/'
        try:
            Heterogeneous_data_files = eval(__config__.get('seq data', 'heter_files'))
            #Heterogeneous_data_files = __config__.get('seq data', 'heter_files')
            if Heterogeneous_data_files is None:
                logging.warning("Run.cfg didn't set heter_files list")
                is_heterfile = False
            elif len(Heterogeneous_data_files) == 0:
                logging.warning("Run.cfg didn't set heter_files list")
                is_heterfile = False
        except Exception as e:
            logging.error("Check your heter_files format.")
            print("heter_files format Error!")
            exit(0)

        if is_heterfile:
            HeteroNewFileList = []
            for heterofile in Heterogeneous_data_files:
                heterofiledir = heter_prefix_path + heterofile
                temp = HeteroNewFiles(heterofiledir)
                if temp is None:
                    return
                HeteroNewFileList.append(temp)

        # filename = IDBABestList[0].split('/')[-2]

        cmdline = 'cat '
        for IDBAFile in IDBABestList:
            cmdline = cmdline + IDBAFile + ' '

        if is_heterfile:
            for newfile in HeteroNewFileList:
                cmdline = cmdline + newfile + ' '

        cmdline = cmdline + '> ./temp_output/total.fa'

        try:
            subprocess.call(cmdline,shell=True)
        except Exception as e:
            logging.error("CAT Command line is wrong")
            exit(0)

        if not os.path.exists('./temp_output/total.fa'):
            logging.error("Must be something wrong here!")
            exit(0)

        logging.info("CAT Completed")
        return './temp_output/total.fa'
    except Exception as e:
        logging.error("CAT Function Wrong:" + str(e))
        print("CAT Combine Wrong",e)
        exit(0)

# CD-HIT-EST Function
def CDHITEST(CATFile):
    logging.info("Start Runing CD-HIT-EST")
    if CATFile is None or ('.fa' not in CATFile):
        logging.error("Can't find the CAT files")
        exit(0)
    try:
        c = __config__.get("CD-hit-est","c")
        aS = __config__.get("CD-hit-est","aS")
        T = __config__.get("CD-hit-est","T")

        if not os.path.exists(CATFile):
            logging.error(str(CATFile) + " not exit.")
            exit(0)

        chefunction = "cd-hit-est -i &inputfile -o &outputfile -c &c -aS &aS -T &T -M 0".replace("&inputfile",CATFile).replace("&outputfile","./temp_output/CD-HIT").replace("&c",c).replace("&aS",aS).replace("&T",T)
        # print(chefunction)

        try:
            subprocess.call(chefunction,shell=True)
        except Exception as e:
            logging.error("CD-HIT-EST Command line is wrong.")
            exit(0)

        if not os.path.exists('./temp_output/CD-HIT'):
            logging.error("Didn't generate the CD-HIT file")
            exit(0)

        logging.info("CD-HIT-EST Completed")
        return './temp_output/CD-HIT'

    except Exception as e:
        logging.error("CD-HIT-EST Function Wrong:" + str(e))
        print("CD-HIT-EST Wrong",e)
        exit(0)

# CAP3 Funtion
def CAP3(CHEFile):
    logging.info("Start Running CAP3:")
    if CHEFile is None:
        logging.error("Can't find the CD-HIT-EST files")
        exit(0)
    try:
        if not os.path.exists("./temp_output/CAP"):
            os.makedirs("./temp_output/CAP")

        filename = CHEFile.split('/')[-1]
        newFiledir = './temp_output/CAP/'+ filename
        shutil.copy(CHEFile,newFiledir)

        cmdline = 'cap3 ' + newFiledir
        subprocess.call(cmdline,shell=True)

        if not os.path.exists(newFiledir+'.cap.contigs') or not os.path.exists(newFiledir+'.cap.singlets'):
            logging.error("cap files are not exits.")
            exit(0)

        cmdline = 'cat ' + newFiledir+'.cap.contigs' + ' ' + newFiledir+'.cap.singlets' + ' > ./output/Final.fa'
        try:
            subprocess.call(cmdline,shell=True)
        except Exception as e:
            logging.error("CAP3 Command line is wrong.")
            exit(0)

        # print(cmdline)

        logging.info("CAP3 Completed")
        return './output/Final.fa'

    except Exception as e:
        logging.error("CAP3 Function Wrong:" + str(e))
        print("CAP3 Wrong",e)
        exit(0)

# Refine Function
def Refine(FinalFile):
    try:
        removelen = __config__.get("transcript","removelen")
        try:
            if 'bp' not in removelen:
                logging.error("removelen must be interger + bg. like 300bp.")
                print("removelen must be interger + bg. like 300bp")
                exit(0)
            removelen = int(removelen[:-2])
            print(removelen)
        except Exception as e:
            logging.error("removelen must be interger + bg. like 300bp.")
            print("removelen must be interger + bg. like 300bp")
            exit(0)
        newFileName = './output/transcript.fa'
        fp = open(FinalFile,'r')
        wr = open(newFileName,'w')
        line = fp.readline()
        tempname = 'ThisisByAnt'
        tempGene = []
        tempGeneNum = 0

        if os.path.exists('./temp_output/HeteroFile/'):
            HeteroFiles = [ './temp_output/HeteroFile/' + x for x in os.listdir('./temp_output/HeteroFile/')]
            is_HeteroFiles = True
        else:
            is_HeteroFiles = False
        # print(HeteroFiles)
        while line:
            if '>' in line:
                # print(tempGeneNum)
                if tempGeneNum < removelen:
                    # print("a:")
                    tempname = line
                    line = fp.readline()
                    continue
                if is_HeteroFiles:
                    for file in HeteroFiles:
                        # print("b:")
                        fr = open(file,'r')
                        if line in fr.read():
                            tempname = line
                            fr.close()
                            break
                        fr.close()
                # print(tempname,line)
                if tempname != line:
                    # print("c:")
                    wr.write(tempname)
                    wr.writelines(tempGene)
                    tempname = line
                tempGene = []
                tempGeneNum = 0
            else:
                tempGeneNum = tempGeneNum + len(line) - 1
                tempGene.append(line)
            line = fp.readline()
        wr.close()
        fp.close()
        return True

    except Exception as e:
        logging.error("Refine Function Wrong:" + str(e))
        print("Refine Wrong",e)
        exit(0)

if __name__=='__main__':

    global start
    start = time.time()

    CheckFileDirs()
    loggingtime('ChenkFileDirs', start)

    IDBABestFiles = IDBA()
    loggingtime('IDBA', start)

    Totalfile = CATCombine(IDBABestFiles)
    loggingtime('CAT', start)

    CHEFile = CDHITEST(Totalfile)
    loggingtime('CD-HIT-EST', start)

    FinalFile = CAP3(CHEFile)
    loggingtime('CAP3', start)

    defaultFunction = __config__.get("transcript","refinement")
    if defaultFunction.lower() != 'yes' and FinalFile is not None:
        if os.path.exists(FinalFile):
            os.rename(FinalFile,"./output/transcript.fa")
        print("The file is integrated in the current output directory.")
        #shutil.rmtree('./temp_output')

    elif defaultFunction.lower() == 'yes' and FinalFile is not None:
        if Refine(FinalFile):
            if os.path.exists(FinalFile):
                os.remove(FinalFile)
            print("The file is integrated in the current output directory.")
            #shutil.rmtree('./temp_output')
    else:
        logging.error("Output Wrong")
        print("Output Wrong")
    loggingtime('Final', start)
    logging.info("Success!")


