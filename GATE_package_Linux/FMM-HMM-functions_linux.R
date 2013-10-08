Inputdata<-function(name,nmark,ntime)
{
  data=read.table(name,sep='\t',header=T)
  location=data[,1:3]
  colnames(location)=c("chr","start","end")
  average = mean(as.matrix(log(data[,4:dim(data)[2]]+1)))
  constant=floor(10/average)
  cat ("\n Average log-value for all data is:",average)
  cat ("\n Transform data using 'constant * log (x+1)', here constant is: ",constant,"\n\n")
  data=constant*log(data[,4:dim(data)[2]]+1)
  data=floor(data)
  observation=array(0,c(dim(data)[1],ntime,nmark))
  for (i in 1:nmark) { for (j in 1:ntime) {observation[,j,i]=data[,(nmark*j-nmark+i)] } }

  return (list(location=location,observation=observation))
}


FMM.HMM.program <-function (observation,ncluster,maxiteration,nstep,ndistance,initial) 
{
set.seed(20110617)
nregion=dim(observation)[1]
#ncluster=8
nmarker=dim(observation)[3]
ntime=dim(observation)[2]
nstate=2
#maxiteration=6000
#nstep=30
#ndistance=0.005
largenum=1E50

# transition probability: b
b=array(1, dim=c(ncluster,2,2))
for (k in 1:ncluster) 
{
  b[k,1,1]=0.5
  b[k,2,1]=0
  b[k,1,2]=1-b[k,1,1]
  b[k,2,2]=1
}

if (initial==1)
   { cluster<-kmeans(observation[,1,],ncluster)$cluster }
else
   { cluster<-sample(gl(ncluster,4,nregion)) }

lambda=array(1, dim=c(ncluster,2,nmarker))
for (k in 1:ncluster) 
{
  for (m in 1:nmarker)
  {
     lambda[k,1,m]=mean(observation[cluster==k,1,m])
     lambda[k,2,m]=mean(observation[cluster==k,ntime,m])
  }
}

## Basic R functions
initModel <-function (b,lambda,cluster)
{
  pie=NULL
  for (k in 1:ncluster)
  {
    pie=cbind(pie,length(cluster[cluster==k])/nregion)
  }
  return (list(b=b,lambda=lambda,cluster=cluster,pie=pie))
}


###########################################################################################################


## Initial values
FMHMM = initModel(b, lambda, cluster)
temp.FMHMM = FMHMM

#for.ward=r_forward(temp.FMHMM,observation)
#back.ward=r_backward(temp.FMHMM,observation)
#temp.FMHMM=r_cluster_matrix(temp.FMHMM,for.ward)
#probability=temp.FMHMM$probability

## C variable input
observation_c<-as.vector(observation)
b_c<-as.vector(temp.FMHMM$b)
lambda_c<-as.vector(temp.FMHMM$lambda)
cluster_c<-as.numeric(temp.FMHMM$cluster)
pie_c<-as.vector(temp.FMHMM$pie)
probability_c<-rep(0,nregion*ncluster)
forwar_c<-rep(0,nregion*ntime*2*ncluster)
backwar_c<-rep(0,nregion*ntime*2*ncluster)
diff_r<-rep(0,maxiteration+10)
pie_r<-rep(0,(ncluster+1)*maxiteration)
H_r<-rep(0, nregion*ntime+10)
log_likeli_r<-rep(0, maxiteration+10)
stop<-1

dyn.load("FMM-HMM.so")
results<-NULL
results<-.C("f_HMM",observation=as.double(observation_c),b=as.double(b_c),lambda=as.double(lambda_c),cluster=as.integer(cluster_c), 
piee=as.double(pie_c),forwar=as.double(forwar_c),backwar=as.double(backwar_c), probability = as.double(probability_c), diff_r=as.double(diff_r),
log_likeli_r=as.double(log_likeli_r),pie_r=as.double(pie_r),as.integer(nstate),as.integer(ncluster),as.integer(nmarker),as.integer(ntime),
as.integer(nregion),as.integer(maxiteration),stop=as.integer(stop),as.integer(nstep),as.double(ndistance),hidden=as.integer(H_r))

dyn.unload("FMM-HMM.so")

results$lambda=matrix(results$lambda,nrow=ncluster,ncol=2*nmarker,byrow=F)
colnames(results$lambda)=as.vector(rbind(paste("Mark",1:nmarker,"_HS",0,sep=''),paste("Mark",1:nmarker,"_HS",1,sep='')))
rownames(results$lambda)=1:ncluster
results$probability=matrix(results$probability,nrow=nregion,ncol=ncluster,byrow=F)
pdf(file=paste("detection of process_cluster",ncluster,".pdf",sep=''))
i=results$stop
log_likeli=results$log_likeli_r[1:i]
pie=matrix(results$pie_r[1:(i*ncluster)],nrow=ncluster,ncol=i,byrow=F)
par(mar=c(5, 4, 4, 8) + 0.1)
plot(1:i,log_likeli/min(log_likeli),type='l',pch=1,col="black",ylim=c(0,1),yaxt='n',xlab='# iteration',ylab='-log_likelihood')
for (j in 1:ncluster)
{
  lines(1:i,pie[j,],type='l',pch=2,col=j+1)   
}
axis(2,at=(1:20)/20,labels=(1:20)/20,col.axis="black",las=2)
axis(4,at=(1:20)/20,labels=(1:20)/20,col.axis="red",las=2)
mtext("ratio of size in each cluster",side=4,las=0,line=3,col='red')
dev.off()
results$hidden=matrix(results$hidden[1:(nregion*ntime)],nrow=nregion,ncol=ntime,byrow=F)
colnames(results$hidden)=paste("HS_t",1:ncol(results$hidden),sep='')
return (results)

}
 
grouping.clusters <-function(location,result,group_num)
{
  hclust_result =hclust(d = dist(result$lambda), method = "ave")
  groups=cutree(hclust_result, k=group_num)
  output=cbind(location,results$cluster,groups[results$cluster])
  colnames(output)=c("chr","start","end","cluster","group")
  result_cl.gr.hid=cbind(output,result$hidden)
  colnames(result_cl.gr.hid)=c("chr","start","end","cluster","group",paste("HS_t",1:ncol(result$hidden),sep=''))
  write.table(result_cl.gr.hid,"result_cl-gr-hid.txt",quote=F,sep='\t',row.names=F)
  return (output)
}

color.BED <-function(data, hidden)
{
  palette(rainbow(length(unique(data$group))))
  color=data$group;
  for (i in 1:length(unique(data$group)))
  {
    color[data$group==i]=paste(col2rgb(i)[1],col2rgb(i)[2],col2rgb(i)[3],sep=',')
  }
  thick.start=array(0,dim=dim(data)[1])
  ntime=ncol(hidden)
  frag_len=mean(data$end-data$start)
  for (i in 1:ntime) {
  thick.start[rowSums(hidden)==(i-1)]=data$end[rowSums(hidden)==(i-1)]-floor((i-1)*frag_len/ntime)-1
  }
  data_col_BED=cbind(as.character(data$chr),data$start,data$end,as.character(data$cluster),0,'+',thick.start,data$end,color)
  write.table(data_col_BED,paste("data_col_BED","bed",sep='.'),quote=F,sep='\t',row.names=F,col.names=F)
  
}

