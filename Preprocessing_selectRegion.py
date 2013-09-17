#!/usr/bin/python

import sys,os,argparse
import pp
import TableIO

sys.path.append("./")
def ParseArg():
    p=argparse.ArgumentParser( description='Data preprocessing and region selction for GATE')
    p.add_argument('-i','--info_bed',type=str,help='information of aligned bed files to be inputed,must be stored in the folder "bed_file"')
    p.add_argument('-g','--genome',type=str,help='genome information for input data, can be mm9/hg19/hg18/')
    p.add_argument('-p','--parameters',type=int,nargs='+',default=[5,5],help="parameters [m,n] to select regions. those regions with normalized signal no less than 'm' in at least one time point for at least 'n' epi-marks were selected, default:[5,5]")
    p.add_argument('-n','--normalize',type=int,default=1E8,help='normalize total read counts into this number, default:1E8')
    p.add_argument('-o','--output',type=str,help='output txt file ')
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

def chr_info(genome):
    chr_file=open("./chromSizes/"+genome+".txt",'r')
    chr_list={}
    for c in chr_file:
        c=c.strip().split('\t') # chr_name\tlength
        chr_list[c[0]]=int(c[1])
    chr_file.close()
    return chr_list

def BEDtoCount(bed,chr_list,bins,binsize):
    bedF=TableIO.parse("./bed_file/"+bed,'bed')
    total_bin_num=0
    count={}
    for c,length in chr_list.items():
        count[c]=[0]*bins[c]

    read_num=0
    for i in bedF:
        if i.chr not in chr_list.keys(): continue
        start_bin=i.start/binsize
        end_bin=i.stop/binsize
        for n in range(start_bin,(end_bin+1)):
            try:
                count[i.chr][n]+=1
            except IndexError:
                continue
        read_num+=1
    return count,read_num



def main():
    args=ParseArg()
    global cutoff, mark_num, binsize,normalize,time_list,mark_list
    if len(args.parameters)!=2:
        print >> sys.stderr,"length of parameters doesn't hit, choose default: 5,5"
        parameters=[5,5]
    else:
        parameters=args.parameters
    cutoff=parameters[0]
    mark_num=parameters[1]

    binsize=200
    normalize=args.normalize
    genome=args.genome
    chr_list=chr_info(genome)

    # get bed list
    files=open(args.info_bed)
    bed_list={}
    time_list=[]
    mark_list=[]
    for f in files:
        f=f.strip().split('\t')
        if f[0] not in time_list: time_list.append(f[0])
        if f[1] not in mark_list: mark_list.append(f[1])
        bed_list[(f[0],f[1])]=f[2]
    files.close()

    # get bin_nums
    bins={}
    for c,length in chr_list.items():
        bins[c]=length/binsize
        if length%binsize!=0: bins[c]+=1

    # parallel counting the reads
    print >> sys.stderr,"# Start counting reads for binsize: %dbp"%binsize
    ppservers=()
    job_server=pp.Server(ppservers=ppservers)
    counts={}
    read_num={}
    for key,bedfile in bed_list.items():
        print "  # counting for file: "+key[0]+"_"+key[1]
        counts[key],read_num[key]=BEDtoCount(bedfile,chr_list,bins,binsize)
        #job1=job_server.submit(BEDtoCount,(bedfile,chr_list,bins,binsize),(),('TableIO',))
        #counts[key],read_num[key]=job1()

    output=open(args.output,"w")
    print "# Select regions:"
    #select regions
        #print header
    header="chr\tstart\tend"
    for time in time_list:
        for mark in mark_list:
            header+=('\t'+mark+"_"+time)
    print >>output,header
    for c in chr_list.keys():
        print >>sys.stderr,"  # Select regions for "+c
        for i in range(bins[c]):
            count_list=[0]*(len(mark_list)*len(time_list))
            exist_mark_num=0
            for m,mark in enumerate(mark_list):
                exist=0
                for t,time in enumerate(time_list):
                    key=(time,mark)
                    norm_count=float(counts[key][c][i])*normalize/read_num[key]
                    count_list[t*len(mark_list)+m]=norm_count
                    if norm_count>=cutoff: exist=1
                exist_mark_num+=exist
                if exist_mark_num-(m+1)<mark_num-len(mark_list): break
            if (exist_mark_num>=mark_num) and (sum(count_list)/len(count_list)<600): 
                # 600 to avoid extremely high value for all marks
                print >>output,c+"\t%d\t%d\t"%(i*binsize+1,(i+1)*binsize)+"\t".join('%.4f'%f for f in count_list)
    output.close()

if __name__=="__main__":
    main()






